#include "utils.h"
#include <cstdio>

#if __has_include(<zlib.h>)
    #include <zlib.h>
    #define HALIGN4_HAVE_ZLIB 1
#else
    #define HALIGN4_HAVE_ZLIB 0
#endif

#include "kseq.h"

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
    struct KseqReader::Impl
    {

#if HALIGN4_HAVE_ZLIB
        gzFile fp{nullptr};
#else
        std::FILE* fp{nullptr};
#endif
        kseq_t* seq{nullptr};
        FilePath file_path; // store the input path for diagnostics and error messages
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

    KseqReader::KseqReader(const FilePath& file_path)
        : impl_(std::make_unique<Impl>())
    {
        impl_->file_path = file_path;

#if HALIGN4_HAVE_ZLIB
        // gzopen/gzread 是 kseq 常用组合。:contentReference[oaicite:2]{index=2}
        // 额外说明：gzopen 在读非 gzip 文件时也会直接读取（不解压），所以一套 reader 可覆盖 .gz / 非 .gz。:contentReference[oaicite:3]{index=3}
        impl_->fp = gzopen(file_path.string().c_str(), "rb");
        if (!impl_->fp) {
            throw makeIoError("failed to open input", file_path);
        }
        impl_->seq = kseq_init(impl_->fp);
#else
        impl_->fp = std::fopen(file_path.string().c_str(), "rb");
        if (!impl_->fp) {
            throw makeIoError("failed to open input", file_path);
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

    void FastaWriter::write(const SeqRecord& rec)
    {
        if (!out_) throw std::runtime_error("FastaWriter output stream is not ready");

        out_.put('>');
        out_.write(rec.id.data(), static_cast<std::streamsize>(rec.id.size()));
        if (!rec.desc.empty()) {
            out_.put(' ');
            out_.write(rec.desc.data(), static_cast<std::streamsize>(rec.desc.size()));
        }
        out_.put('\n');

        writeWrapped(out_, rec.seq, line_width_);
    }

    void FastaWriter::flush()
    {
        out_.flush();
    }

} // namespace seq_io
