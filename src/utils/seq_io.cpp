// 该文件负责 FASTA/FASTQ 的读写（基于 kseq），并包含若干针对大文件/高吞吐量的优化说明。
// 关键点（性能说明，中文）：
// - kseq 是一个轻量的基于 stdio/gzread 的读取器，常用于高性能序列读取。它本身在单次读取上已经很高效，
//   但在面对非常大的文件（GB 级或更大）时，IO 缓冲和系统调用次数会成为瓶颈。
// - 对策：对于未压缩的文件，我们通过为 std::FILE* 设置一个较大的用户空间缓冲（setvbuf + malloc）来
//   减少 fread 的系统调用次数；对于压缩文件（zlib），我们可以通过 gzbuffer 提升 zlib 内部缓冲。
// - 这些优化都有权衡：更大的缓冲减少系统调用但占用更多内存；如果内存紧张，应降低缓冲大小。
// - 此外，写出端（FastaWriter）也尽量减少小的写调用，构造好一次性的 header/序列缓冲并批量写出，
//   避免 per-character 的写入，显著提升写吞吐。

#include "utils.h"
#include <cstdio>
#include <cstdlib> // for malloc/free

#if __has_include(<zlib.h>)
    #include <zlib.h>
    #define HALIGN4_HAVE_ZLIB 1
#else
    #define HALIGN4_HAVE_ZLIB 0
#endif

#include "kseq.h"

// kseq 初始化：根据是否有 zlib 选择不同的底层读函数。
// 说明：kseq 的 KSEQ_INIT 宏会生成 kseq_read / kseq_init / kseq_destroy 等函数。
// - 当有 zlib 时使用 gzFile/gzread；此时可以通过 gzbuffer 调整 zlib 的内部缓冲大小以提升吞吐。
// - 当没有 zlib 时使用 std::FILE* 与自定义的 fileRead（调用 fread），此时我们通过 setvbuf 来替换 stdio 的缓冲。

#if HALIGN4_HAVE_ZLIB
    KSEQ_INIT(gzFile, gzread)
#else
    static int fileRead(std::FILE* fp, void* buf, int len)
    {
        const std::size_t n = std::fread(buf, 1, static_cast<std::size_t>(len), fp);
        return static_cast<int>(n);
    }
    KSEQ_INIT(std::FILE*, fileRead)
#endif

namespace seq_io
{
    // Impl 结构体：封装底层文件描述符、kseq 对象以及用于加大缓冲的持久化缓冲区指针（仅在非 zlib 时使用）。
    // 设计说明（性能相关）：
    // - io_buf 与 io_buf_size 用于在 fread 场景中提升 stdio 缓冲到用户自定义大小（例如 1 MiB），
    //   这样可以显著减少内核态系统调用次数（read/fread），在顺序读取大文件时非常有效。
    // - 对于压缩文件（zlib），我们不能替换 zlib 内部的缓冲，但可以使用 gzbuffer 来告诉 zlib 增大内部缓冲区，
    //   以减少对底层 read 的调用频率与解压函数调用开销。
    // - 这些缓冲策略对 SSD/HDD/网络文件系统的效果不同：SSD 更能从更大缓冲中受益，网络文件系统需谨慎以避免占用大量内存。
    struct KseqReader::Impl
    {

#if HALIGN4_HAVE_ZLIB
        gzFile fp{nullptr};
#else
        std::FILE* fp{nullptr};
#endif
        kseq_t* seq{nullptr};
        FilePath file_path; // store the input path for diagnostics and error messages

        // 新增：用于加大 stdio / zlib 缓冲的持久化缓冲区（非 zlib 情况）
        // 这样可以减少 fread/gzread 的系统调用次数，提高大文件读取吞吐量
        char* io_buf{nullptr};
        std::size_t io_buf_size{1 << 20}; // 默认 1 MiB，可根据需要调整
    };

    static std::runtime_error makeIoError(const std::string& msg, const FilePath& p)
    {
        return std::runtime_error(msg + ": " + p.string());
    }

    static void assignKstring(std::string& dst, const kstring_t& ks)
    {
        if (ks.s && ks.l > 0) dst.assign(ks.s, ks.l);
        else dst.clear();
    }

    // KseqReader 构造/析构：在这里对 IO 缓冲做了优化
    // 关键点：
    // - 对于 gzip 文件（HALIGN4_HAVE_ZLIB），使用 gzopen/gzbuffer/gzread；gzbuffer 的 size 参数是 int，
    //   因此不要将 io_buf_size 设得过大以避免接口问题。建议 256 KiB ~ 4 MiB 根据内存与硬件调优。
    // - 对于未压缩文件，使用 fopen + setvbuf，分配用户空间缓冲区并设为 _IOFBF（全缓冲）。
    // - 为什么不总是使用更大缓冲？因为内存占用和系统资源是稀缺的：对多线程/多流同时读取场景，过大缓冲会导致内存压力。
    KseqReader::KseqReader(const FilePath& file_path)
        : impl_(std::make_unique<Impl>())
    {
        impl_->file_path = file_path;

#if HALIGN4_HAVE_ZLIB
        // gzopen/gzread 是 kseq 常用组合。
        impl_->fp = gzopen(file_path.string().c_str(), "rb");
        if (!impl_->fp) {
            throw makeIoError("failed to open input", file_path);
        }

        // 提升 zlib 的内部缓冲区大小以加速大文件读取（单位：字节，接口为 int）
        // 如果 io_buf_size 超过 int 上限，gzbuffer 的行为未定义；这里选用 1MiB，安全且有效。
        gzbuffer(impl_->fp, static_cast<int>(impl_->io_buf_size));

        impl_->seq = kseq_init(impl_->fp);
#else
        impl_->fp = std::fopen(file_path.string().c_str(), "rb");
        if (!impl_->fp) {
            throw makeIoError("failed to open input", file_path);
        }

        // 为 stdio 分配一个较大的缓冲区并设置为全缓冲，减少 fread 的系统调用次数
        // 注意：setvbuf 要在首次 I/O 之前调用，且缓冲区要一直保持直到 fclose
        impl_->io_buf = static_cast<char*>(std::malloc(impl_->io_buf_size));
        if (impl_->io_buf) {
            if (setvbuf(impl_->fp, impl_->io_buf, _IOFBF, static_cast<size_t>(impl_->io_buf_size)) != 0) {
                // 如果设置失败，释放并继续（仍能工作，只是没有额外加速）
                std::free(impl_->io_buf);
                impl_->io_buf = nullptr;
            }
        }

        impl_->seq = kseq_init(impl_->fp);
#endif

        if (!impl_->seq) {
            throw makeIoError("failed to init kseq", file_path);
        }
    }

    KseqReader::~KseqReader()
    {
        if (!impl_) return;

        if (impl_->seq) {
            kseq_destroy(impl_->seq);
            impl_->seq = nullptr;
        }

#if HALIGN4_HAVE_ZLIB
        if (impl_->fp) {
            gzclose(impl_->fp);
            impl_->fp = nullptr;
        }
#else
        if (impl_->fp) {
            std::fclose(impl_->fp);
            impl_->fp = nullptr;
        }
        // 释放为 stdio 分配的缓冲区（必须在 fclose 之后）
        if (impl_->io_buf) {
            std::free(impl_->io_buf);
            impl_->io_buf = nullptr;
        }
#endif
    }

    KseqReader::KseqReader(KseqReader&& other) noexcept = default;
    KseqReader& KseqReader::operator=(KseqReader&& other) noexcept = default;

    bool KseqReader::next(SeqRecord& rec)
    {
        if (!impl_ || !impl_->seq) {
            throw std::runtime_error("KseqReader is not initialized");
        }

        const int ret = kseq_read(impl_->seq);

        if (ret >= 0) {
            assignKstring(rec.id,   impl_->seq->name);
            assignKstring(rec.desc, impl_->seq->comment);
            assignKstring(rec.seq,  impl_->seq->seq);
            return true;
        }

        if (ret == -1) {
            return false; // EOF
        }

        // 其它负值视为解析错误（如 FASTQ 质量串截断等场景）
        throw std::runtime_error("kseq_read() failed with code " + std::to_string(ret) +
                                 " for file: " + impl_->file_path.string());
    }

    // ------------------------- FastaWriter -------------------------

    FastaWriter::FastaWriter(const FilePath& file_path, std::size_t line_width)
        : out_(file_path, std::ios::binary), line_width_(line_width)
    {
        if (!out_) {
            throw makeIoError("failed to open output", file_path);
        }
        if (line_width_ == 0) line_width_ = 80;
    }

    void FastaWriter::writeWrapped(std::ofstream& out, std::string_view s, std::size_t width)
    {
        for (std::size_t i = 0; i < s.size(); i += width) {
            const std::size_t n = (i + width <= s.size()) ? width : (s.size() - i);
            out.write(s.data() + static_cast<std::streamoff>(i), static_cast<std::streamsize>(n));
            out.put('\n');
        }
    }

    // FastaWriter::write 的实现包含了若干写优化：
    // - 先把 header 构造为一个小字符串并一次写出，避免多次写调用；
    // - 将序列按行宽组装到一个临时缓冲 seqbuf 中，然后一次性写出，避免 per-character/line 的多次系统调用；
    // - 这种策略在写入巨型 fasta 文件时能显著减少 I/O 开销，但会占用额外的内存用于 seqbuf（其大小 ≈ L + L/width）。
    // - 如果需要进一步提升写性能，可以考虑：
    //     * 使用 unbuffered write（低层 write）并管理自己的较大缓冲区；
    //     * 使用多线程生产者-消费者模型：写线程负责把已经准备好的缓冲写入磁盘，计算线程负责生成序列字符串。
    //   这些会增加实现复杂度（并发、同步与错误处理）。
    // - 小提示：通用做法是把写缓冲大小设置为几百 KB 到几 MB，根据目标磁盘与系统内存调整。

    // cpp
    void FastaWriter::write(const SeqRecord& rec)
    {
        if (!out_) throw std::runtime_error("FastaWriter output stream is not ready");

        // 构造并一次性写入 header
        std::string header;
        header.reserve(1 + rec.id.size() + (rec.desc.empty() ? 0 : 1 + rec.desc.size()) + 1);
        header.push_back('>');
        header.append(rec.id);
        if (!rec.desc.empty()) {
            header.push_back(' ');
            header.append(rec.desc);
        }
        header.push_back('\n');
        out_.write(header.data(), static_cast<std::streamsize>(header.size()));

        // 构造带换行的序列缓冲区并一次性写出，减少多次 write/put 调用
        const std::size_t L = rec.seq.size();
        const std::size_t width = (line_width_ == 0) ? 80 : line_width_;

        if (L == 0) {
            out_.put('\n');
            return;
        }

        std::string seqbuf;
        seqbuf.reserve(L + (L / width) + 1);
        for (std::size_t i = 0; i < L; i += width) {
            const std::size_t n = (i + width <= L) ? width : (L - i);
            seqbuf.append(rec.seq.data() + static_cast<std::size_t>(i), n);
            seqbuf.push_back('\n');
        }

        out_.write(seqbuf.data(), static_cast<std::streamsize>(seqbuf.size()));
    }

    void FastaWriter::flush()
    {
        out_.flush();
    }

} // namespace seq_io
