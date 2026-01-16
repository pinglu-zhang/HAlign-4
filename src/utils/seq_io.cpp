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

// kseq 初始化：根据是否有 zlib 选择不同的底层句柄类型与 read 函数。
//
// 这里最容易“看不懂”的点是 KSEQ_INIT：
// - KSEQ_INIT(T, readfn) 会基于句柄类型 T 和读取函数 readfn，生成一组静态函数/类型：
//   * kseq_t         : 解析器状态（内部包含 name/comment/seq/qual 等 kstring）
//   * kseq_init(T)   : 用句柄创建解析器
//   * kseq_read(kseq_t*) : 读取下一条 FASTA/FASTQ 记录（返回码见下方 next() 契约）
//   * kseq_destroy(kseq_t*) : 释放解析器
//
// 我们的策略：
// - 有 zlib：用 gzFile + gzread，能读 .gz；并通过 gzbuffer() 增大 zlib 内部 buffer。
// - 无 zlib：用 FILE* + 自定义 fileRead(fread)，并通过 setvbuf() 增大 stdio buffer。

#if HALIGN4_HAVE_ZLIB
    KSEQ_INIT(gzFile, gzread)
#else
    static int fileRead(std::FILE* fp, void* buf, int len)
    {
        // 约定：kseq 期待 readfn 返回“实际读取到的字节数”。
        // - 返回 0 表示 EOF；
        // - 返回正数表示读取成功；
        // - 这里不返回负数（stdio 的错误通过 ferror 暗含），kseq 会按读取结果推进。
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
        // 资源所有权说明：
        // - fp：底层文件句柄（gzFile 或 FILE*），由 KseqReader 独占并在析构中关闭。
        // - seq：kseq 解析器状态，依赖 fp 的有效性；必须先 destroy(seq) 再 close(fp)。
        // - io_buf：仅在非 zlib 情况下用于 setvbuf() 的用户态缓冲；其生命周期必须覆盖整个读取过程。

#if HALIGN4_HAVE_ZLIB
        gzFile fp{nullptr};
#else
        std::FILE* fp{nullptr};
#endif
        kseq_t* seq{nullptr};
        FilePath file_path; // store the input path for diagnostics and error messages

        // io_buf/io_buf_size：用于提升底层顺序读取吞吐的“大块缓冲”。
        // 性能动机：
        // - 顺序扫描巨大 FASTA/FASTQ 时，小 buffer 会导致 read()/fread() 调用次数激增，系统调用开销明显；
        // - 通过更大的 buffer，可以让底层一次性把更多数据搬到用户态，kseq 再在内存中解析。
        // 代价与权衡：
        // - 每个 Reader 实例都会分配一块 buffer；多线程并行读取多个文件时，内存占用会线性放大。
        // - 因此这里用 8MiB 作为吞吐与内存的折中默认值。
        char* io_buf{nullptr};
        std::size_t io_buf_size{8 << 20}; // 默认 8 MiB，可根据需要调整
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

    // KseqReader 构造/析构：在这里对 IO 缓冲做 best-effort 优化。
    //
    // 异常安全说明：
    // - 任何“无法继续工作”的错误（open 失败、kseq_init 失败）都会抛异常；
    // - buffer 调整（gzbuffer/setvbuf）失败时不会抛异常：因为这只影响性能，不影响正确性。
    //
    // 性能说明：
    // - 对 gzip：gzbuffer() 能减少 zlib 内部 refill 次数（以及潜在的系统调用/解压调度开销）。
    // - 对纯文本：setvbuf() 能显著减少 fread() 触发的内核读次数。
    KseqReader::KseqReader(const FilePath& file_path)
        : impl_(std::make_unique<Impl>())
    {
        impl_->file_path = file_path;

#if HALIGN4_HAVE_ZLIB
        // gzopen/gzread 是 kseq 的常用组合：
        // - 读 .gz 场景无需上层做任何区分；
        // - kseq 始终通过 gzread 拉取原始字节流并解析成记录。
        impl_->fp = gzopen(file_path.string().c_str(), "rb");
        if (!impl_->fp) {
            throw makeIoError("failed to open input", file_path);
        }

        // 提升 zlib 的内部缓冲区大小以加速大文件读取（接口参数为 int）。
        // 这里用 io_buf_size(默认 1MiB) 并 static_cast<int>。
        // 说明：即使 gzbuffer 调用失败，gzread 仍然可用，所以这里不检查返回值，保持 best-effort。
        gzbuffer(impl_->fp, static_cast<int>(impl_->io_buf_size));

        impl_->seq = kseq_init(impl_->fp);
#else
        impl_->fp = std::fopen(file_path.string().c_str(), "rb");
        if (!impl_->fp) {
            throw makeIoError("failed to open input", file_path);
        }

        // 为 stdio 分配一个较大的用户态缓冲并设置为全缓冲，减少 fread 的系统调用次数。
        // 注意：
        // - setvbuf 必须在首次 I/O 前调用；
        // - 缓冲区内存在 fclose 前必须保持有效；
        // - 如果 setvbuf 失败，这个 reader 仍然可以正常工作，只是吞吐可能下降。
        impl_->io_buf = static_cast<char*>(std::malloc(impl_->io_buf_size));
        if (impl_->io_buf) {
            if (setvbuf(impl_->fp, impl_->io_buf, _IOFBF, static_cast<size_t>(impl_->io_buf_size)) != 0) {
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
        // 析构必须不抛异常：I/O 清理失败最多影响资源回收，不应在栈展开期间触发 terminate。
        if (!impl_) return;

        // 释放顺序很关键：kseq_t 内部可能持有指向底层缓冲/句柄的状态，必须先 destroy 再 close。
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
        // 释放为 stdio 分配的用户态缓冲（必须在 fclose 之后）。
        if (impl_->io_buf) {
            std::free(impl_->io_buf);
            impl_->io_buf = nullptr;
        }
#endif
    }

    // 移动构造/赋值标注 noexcept：
    // - 允许 KseqReader 放入 std::vector 等容器时走 move 而非 copy；
    // - move 期间不应抛异常，符合“资源句柄转移”的语义。
    KseqReader::KseqReader(KseqReader&& other) noexcept = default;
    KseqReader& KseqReader::operator=(KseqReader&& other) noexcept = default;

    bool KseqReader::next(SeqRecord& rec)
    {
        if (!impl_ || !impl_->seq) {
            throw std::runtime_error("KseqReader is not initialized");
        }

        // 返回契约（与 kseq 保持一致）：
        // - ret >= 0：成功读到一条记录；本函数会覆盖 rec.id/desc/seq 并返回 true。
        // - ret == -1：EOF；返回 false。注意：EOF 时 rec 保持调用者上一次的内容（本实现不清空）。
        // - ret < -1：解析错误（例如 FASTQ 质量串长度不匹配/截断等）；抛异常并带文件路径便于定位。
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

        throw std::runtime_error("kseq_read() failed with code " + std::to_string(ret) +
                                 " for file: " + impl_->file_path.string());
    }

    // ------------------------- SeqWriter -------------------------

    SeqWriter::SeqWriter(const FilePath& file_path, std::size_t line_width)
        : SeqWriter(file_path, Format::fasta, line_width, 8ULL * 1024ULL * 1024ULL)
    {}

    SeqWriter::SeqWriter(const FilePath& file_path, std::size_t line_width, std::size_t buffer_threshold_bytes)
        : SeqWriter(file_path, Format::fasta, line_width, buffer_threshold_bytes)
    {}

    SeqWriter::SeqWriter(const FilePath& file_path, Format fmt, std::size_t line_width, std::size_t buffer_threshold_bytes)
        : out_(file_path, std::ios::binary),
          format_(fmt),
          line_width_(line_width == 0 ? 80 : line_width),
          buffer_threshold_bytes_(buffer_threshold_bytes)
    {
        if (!out_) {
            throw makeIoError("failed to open output", file_path);
        }

        // 自建 buffer 的意义：把大量小 append 聚合成“大块 write”，减少系统调用次数。
        // - buffer_threshold_bytes_ == 0：禁用自建 buffer，直接写入 out_（仍有 ofstream 内部缓冲）。
        // - buffer_threshold_bytes_  > 0：先积累到 buffer_，到阈值再一次性 write。对吞吐更友好。
        if (buffer_threshold_bytes_ > 0) {
            // reserve：减少 buffer_ 在增长过程中的多次 realloc/memcpy；上限 8MiB 作为折中。
            buffer_.reserve(std::min<std::size_t>(buffer_threshold_bytes_, 8ULL * 1024ULL * 1024ULL));
        }
    }

    SeqWriter SeqWriter::Sam(const FilePath& file_path, std::size_t buffer_threshold_bytes)
    {
        return SeqWriter(file_path, Format::sam, /*line_width=*/0, buffer_threshold_bytes);
    }

    SeqWriter::~SeqWriter()
    {
        try {
            flush();
        } catch (...) {
            // best-effort
        }
    }

    void SeqWriter::flushBuffer_()
    {
        // flushBuffer_ 是性能关键路径：尽量一次 write 出去，避免频繁的小写。
        if (buffer_.empty()) return;
        out_.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
        buffer_.clear();
    }

    void SeqWriter::appendOrFlush_(std::string_view s)
    {
        if (!out_) throw std::runtime_error("SeqWriter output stream is not ready");

        // buffer_threshold_bytes_ == 0：不做额外聚合，直接写。
        // 说明：这里仍然是“成块写入”的（一次 write s.size()），不是逐字符写。
        if (buffer_threshold_bytes_ == 0) {
            out_.write(s.data(), static_cast<std::streamsize>(s.size()));
            return;
        }

        buffer_.append(s.data(), s.size());
        if (buffer_.size() >= buffer_threshold_bytes_) {
            flushBuffer_();
        }
    }

    void SeqWriter::appendWrapped_(std::string& dst, std::string_view s, std::size_t width)
    {
        // 说明：FASTA 常见的“每行固定宽度”只是可读性要求，不影响序列语义。
        // width==0 时表示不换行（直接原样追加）。
        if (width == 0) {
            dst.append(s.data(), s.size());
            return;
        }

        // 复杂度：O(|s|)，每 width 字符插入一个 '\n'。
        for (std::size_t i = 0; i < s.size(); i += width) {
            const std::size_t n = (i + width <= s.size()) ? width : (s.size() - i);
            dst.append(s.data() + i, n);
            dst.push_back('\n');
        }
    }

    void SeqWriter::writeFasta(const SeqRecord& rec)
    {
        if (format_ != Format::fasta) {
            throw std::runtime_error("SeqWriter::writeFasta called but writer is not in FASTA mode");
        }

        // line_width_ 在构造时已把 0 归一到 80；这里仍保留防御式写法以确保行为稳定。
        const std::size_t width = (line_width_ == 0) ? 80 : line_width_;

        // 这里用 recordbuf 拼出“完整的一条 FASTA 记录”再一次 appendOrFlush_，目的：
        // - 减少多次 appendOrFlush_ 带来的分支/阈值判断开销；
        // - 使底层 write 更集中、系统调用更少。
        // FASTA 格式约定：
        // >{id}[ {desc}]\n
        // {seq (wrap)}\n
        // 边界：seq 为空时，仍写出一个空行作为序列行（保持与许多工具一致，且保持历史行为）。
        std::string recordbuf;
        recordbuf.reserve(1 + rec.id.size() + (rec.desc.empty() ? 0 : 1 + rec.desc.size()) + 1 +
                          rec.seq.size() + (rec.seq.size() / width) + 2);

        // header
        recordbuf.push_back('>');
        recordbuf.append(rec.id);
        if (!rec.desc.empty()) {
            recordbuf.push_back(' ');
            recordbuf.append(rec.desc);
        }
        recordbuf.push_back('\n');

        // sequence (wrapped)
        if (rec.seq.empty()) {
            recordbuf.push_back('\n');
        } else {
            for (std::size_t i = 0; i < rec.seq.size(); i += width) {
                const std::size_t n = (i + width <= rec.seq.size()) ? width : (rec.seq.size() - i);
                recordbuf.append(rec.seq.data() + i, n);
                recordbuf.push_back('\n');
            }
        }

        appendOrFlush_(recordbuf);
    }

    void SeqWriter::writeSamHeader(std::string_view header_text)
    {
        if (format_ != Format::sam) {
            throw std::runtime_error("SeqWriter::writeSamHeader called but writer is not in SAM mode");
        }

        if (header_text.empty()) return;

        // SAM header 允许用户传入多行文本；这里只做一个“确保以换行结尾”的修正，保证后续记录从新行开始。
        appendOrFlush_(header_text);
        if (header_text.back() != '\n') {
            appendOrFlush_("\n");
        }

        // 仅记录状态：本实现不强制要求“先写 header 再写 record”。（保持现有逻辑，仅释义）
        sam_header_written_ = true;
    }

    void SeqWriter::writeSam(const SamRecord& r)
    {
        if (format_ != Format::sam) {
            throw std::runtime_error("SeqWriter::writeSam called but writer is not in SAM mode");
        }

        // 组装一整行 SAM record 再一次性写入：减少多次 write/append 的开销。
        // SAM 共有 11 列必选字段，TAB 分隔：
        // QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
        std::string line;
        line.reserve(r.qname.size() + r.rname.size() + r.cigar.size() + r.rnext.size() + r.seq.size() + r.qual.size() + 64);

        // 注意：这里用 to_string 拼接整数，属于“可读性优先”的实现；
        // 若未来这里成为热路径瓶颈，可考虑 fmt::format_to 或自实现 itoa（但那属于行为等价的性能重构，非本次改动范围）。
        line.append(r.qname.data(), r.qname.size());
        line.push_back('\t');
        line.append(std::to_string(r.flag));
        line.push_back('\t');
        line.append(r.rname.data(), r.rname.size());
        line.push_back('\t');
        line.append(std::to_string(r.pos));
        line.push_back('\t');
        line.append(std::to_string(static_cast<unsigned int>(r.mapq)));
        line.push_back('\t');
        line.append(r.cigar.data(), r.cigar.size());
        line.push_back('\t');
        line.append(r.rnext.data(), r.rnext.size());
        line.push_back('\t');
        line.append(std::to_string(r.pnext));
        line.push_back('\t');
        line.append(std::to_string(r.tlen));
        line.push_back('\t');
        line.append(r.seq.data(), r.seq.size());
        line.push_back('\t');
        line.append(r.qual.data(), r.qual.size());

        // OPT: 可选字段区。为兼容调用方传入的 opt 是否自带前置 '\t'，这里做一次轻量修正。
        if (!r.opt.empty()) {
            if (r.opt.front() != '\t') line.push_back('\t');
            line.append(r.opt.data(), r.opt.size());
        }

        // 每条 record 必须以 '\n' 结束。
        line.push_back('\n');
        appendOrFlush_(line);
    }

    void SeqWriter::flush()
    {
        if (!out_) return;
        flushBuffer_();
        out_.flush();
    }

    // ------------------------------------------------------------------
    // 函数：makeSamRecord
    // 功能：从 SeqRecord 和比对信息构建 SAM 记录
    //
    // 实现说明：
    // 1. 这是一个便捷函数，封装了从 SeqRecord 到 SamRecord 的转换逻辑
    // 2. 质量值处理：若 query.qual 为空，则使用 SAM 默认值 "*"（已在 SamRecord 定义中初始化）
    // 3. 所有 string_view 字段直接引用输入参数，生命周期由调用者保证
    // 4. 性能：O(1)，仅赋值操作，无内存分配或拷贝
    //
    // 参数：
    //   - query: 查询序列的 SeqRecord（必须在 SamRecord 使用期间保持有效）
    //   - ref_name: 参考序列名称（必须在 SamRecord 使用期间保持有效）
    //   - cigar_str: CIGAR 字符串（必须在 SamRecord 使用期间保持有效）
    //   - pos: 1-based 比对起始位置（默认 1）
    //   - mapq: mapping quality（默认 60）
    //   - flag: SAM flag（默认 0）
    //
    // 返回：SeqWriter::SamRecord（所有 string_view 字段引用输入参数）
    // ------------------------------------------------------------------
    SeqWriter::SamRecord makeSamRecord(
        const SeqRecord& query,
        std::string_view ref_name,
        std::string_view cigar_str,
        std::uint32_t pos,
        std::uint8_t mapq,
        std::uint16_t flag)
    {
        // 构建 SAM 记录（所有字段均为 string_view，零拷贝）
        SeqWriter::SamRecord sam_rec;
        sam_rec.qname = query.id;      // query 名称（使用 SeqRecord 的 id 字段）
        sam_rec.flag  = flag;          // SAM flag（0 表示未比对/正向链）
        sam_rec.rname = ref_name;      // 参考序列名称
        sam_rec.pos   = pos;           // 1-based 比对起始位置
        sam_rec.mapq  = mapq;          // mapping quality
        sam_rec.cigar = cigar_str;     // CIGAR 字符串
        sam_rec.seq   = query.seq;     // query 序列

        // 质量值处理：若 query.qual 不为空，则使用；否则保持默认值 "*"
        if (!query.qual.empty()) {
            sam_rec.qual = query.qual;
        }
        // 否则使用默认值 "*"（已在 SamRecord 定义中初始化）

        return sam_rec;
    }

} // namespace seq_io


