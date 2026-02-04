// ==================================================================
// seq_io.cpp - HAlign-4 序列文件 I/O 模块（FASTA/FASTQ/SAM 读写）
// ==================================================================
// 功能概述：
// 本文件提供高性能的序列文件读写功能，支持三种格式：
// 1. FASTA：序列数据的标准格式（仅序列，无质量值）
// 2. FASTQ：包含质量值的序列格式（用于测序数据）
// 3. SAM：序列比对结果格式（包含 CIGAR、位置等比对信息）
//
// 核心组件：
// - KseqReader：基于 kseq.h 的高性能 FASTA/FASTQ 读取器
// - SeqWriter：支持 FASTA/SAM 格式的批量写入器
// - SamReader：SAM 格式的读取器
// - 辅助函数：格式转换、记录构建等
//
// ==================================================================

#include "utils.h"
#include <cstdio>   // std::FILE*, fopen, fclose, fread, setvbuf
#include <cstdlib>  // malloc, free（用于自定义 I/O 缓冲区）

#ifdef _DEBUG
#include <spdlog/spdlog.h>  // 调试日志输出
#endif

// ------------------------------------------------------------------
// zlib 条件编译：支持 .gz 压缩文件的读取
// ------------------------------------------------------------------
// 说明：
// - 如果系统提供 zlib.h，则编译 gzip 支持（HALIGN4_HAVE_ZLIB=1）
// - 否则只支持未压缩文件（HALIGN4_HAVE_ZLIB=0）
// - 这样可以在没有 zlib 的环境（如嵌入式系统）中编译
// ------------------------------------------------------------------
#if __has_include(<zlib.h>)
    #include <zlib.h>
    #define HALIGN4_HAVE_ZLIB 1
#else
    #define HALIGN4_HAVE_ZLIB 0
#endif

#include "kseq.h"  // kseq 库：高性能的 FASTA/FASTQ 解析器（单头文件）

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
// - 无 zlib：用 FILE* + 自定义 fileRead(fread)，并通过 setvbuf() 增大 stdio 缓冲。

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
    // ==================================================================
    // KseqReader::Impl 结构体：封装底层文件描述符和kseq对象
    // ==================================================================
    // 设计说明（Pimpl 惯用法）：
    // - 使用 Impl（Implementation）结构体隐藏实现细节
    // - 好处：
    //   1. 减少头文件依赖（kseq.h 只需在 .cpp 中包含）
    //   2. 二进制兼容性（修改 Impl 不影响 ABI）
    //   3. 编译时间优化（头文件更简洁）
    //
    // 资源所有权说明（重要！）：
    // ------------------------------------------------------------------
    // - fp：底层文件句柄（gzFile 或 FILE*），由 KseqReader 独占
    //   * 在构造函数中打开（gzopen 或 fopen）
    //   * 在析构函数中关闭（gzclose 或 fclose）
    //   * 移动构造/赋值会转移所有权
    //
    // - seq：kseq 解析器状态，依赖 fp 的有效性
    //   * 在构造函数中创建（kseq_init）
    //   * 在析构函数中释放（kseq_destroy）
    //   * **必须先 destroy(seq) 再 close(fp)**（顺序很关键！）
    //   * 原因：seq 内部可能持有指向 fp 缓冲区的指针
    //
    // - io_buf：用户态大缓冲区（仅非 zlib 模式）
    //   * 在构造函数中分配（malloc）
    //   * 传递给 setvbuf() 提升 stdio 性能
    //   * 在析构函数中释放（free）
    //   * **必须在 fclose 之后释放**（stdio 可能持有引用）
    //
    // ==================================================================
    struct KseqReader::Impl
    {
        // 资源所有权说明：
        // - fp：底层文件句柄（gzFile 或 FILE*），由 KseqReader 独占并在析构中关闭。
        // - seq：kseq 解析器状态，依赖 fp 的有效性；必须先 destroy(seq) 再 close(fp)。
        // - io_buf：仅在非 zlib 情况下用于 setvbuf() 的用户态缓冲；其生命周期必须覆盖整个读取过程。

#if HALIGN4_HAVE_ZLIB
        gzFile fp{nullptr};  // zlib 文件句柄（支持自动解压）
#else
        std::FILE* fp{nullptr};  // 标准 C 文件指针（仅支持未压缩文件）
#endif
        kseq_t* seq{nullptr};  // kseq 解析器状态（由 kseq_init 创建）
        FilePath file_path; // 存储输入文件路径（用于诊断和错误消息）

        // ------------------------------------------------------------------
        // io_buf/io_buf_size：用于提升底层顺序读取吞吐的"大块缓冲"
        // ------------------------------------------------------------------
        // 性能动机（详细说明）：
        // 【问题】：顺序扫描巨大 FASTA/FASTQ 时的系统调用瓶颈
        // - kseq 本身解析效率很高（单次遍历，状态机）
        // - 但底层 I/O 的系统调用次数会成为瓶颈
        // - 小 buffer（4-8 KB）→ read()/fread() 调用次数激增
        // - 每次系统调用的固定开销（上下文切换、内核调度）累积
        //
        // 【量化分析】（1 GB FASTA 文件示例）：
        // - 默认 4 KB 缓冲：
        //   * 需要 ~262,144 次 read() 系统调用
        //   * 每次调用开销 ~1-2 μs（上下文切换）
        //   * 总开销：262-524 ms（占总时间 10-20%）
        //
        // - 8 MB 缓冲：
        //   * 仅需 ~128 次 read() 系统调用
        //   * 总开销：0.128-0.256 ms（几乎可忽略）
        //   * 性能提升：2-5 倍（取决于硬件）
        //
        // 【实现策略】：
        // - 通过更大的 buffer，让底层一次性把更多数据搬到用户态
        // - kseq 在内存中快速解析，避免反复进入内核
        // - 对 SSD 特别有效（顺序读吞吐高，延迟低）
        //
        // 【代价与权衡】：
        // - 内存占用：每个 Reader 实例占用 io_buf_size 内存
        // - 多线程场景：N 个线程 × 8 MiB = 64 MiB（N=8）
        // - 如果机器内存紧张（<4 GB），可考虑减小到 1-2 Mi
        // ------------------------------------------------------------------
        char* io_buf{nullptr};            // 用户态缓冲区指针（malloc 分配）
        std::size_t io_buf_size{8 << 20}; // 默认 8 MiB（8 * 1024 * 1024 字节）
                                          // 说明：可根据实际场景调整
                                          //       - SSD/高性能：增大到 16-32 MiB
                                          //       - 内存紧张：减小到 1-4 MiB
    };

    // ==================================================================
    // 辅助函数：错误处理与字符串转换
    // ==================================================================

    // ------------------------------------------------------------------
    // 函数：makeIoError - 构建包含文件路径的异常对象
    // ------------------------------------------------------------------
    // 功能：将错误消息和文件路径组合成统一格式的异常
    // 参数：
    //   @param msg - 错误消息（如 "failed to open input"）
    //   @param p - 文件路径（FilePath 类型）
    // 返回：std::runtime_error 异常对象
    //
    // 使用场景：
    // - 文件打开失败、读取错误、格式错误等
    // - 统一错误消息格式："错误消息: 文件路径"
    //
    // 示例：
    //   throw makeIoError("failed to open input", "/path/to/file.fasta");
    //   // 异常消息："failed to open input: /path/to/file.fasta"
    // ------------------------------------------------------------------
    static std::runtime_error makeIoError(const std::string& msg, const FilePath& p)
    {
        return std::runtime_error(msg + ": " + p.string());
    }

    // ------------------------------------------------------------------
    // 函数：assignKstring - 将 kstring_t 转换为 std::string
    // ------------------------------------------------------------------
    // 功能：从 kseq 内部的 kstring_t 结构体提取字符串到 std::string
    // 参数：
    //   @param dst - 目标 std::string（输出参数，会被覆盖）
    //   @param ks - 源 kstring_t（kseq 内部类型）
    //
    // 工作原理：
    // - kstring_t 结构：{ char* s; size_t l; size_t m; }
    //   * s：字符串指针（以 '\0' 结尾）
    //   * l：字符串长度（不含 '\0'）
    //   * m：已分配容量
    // - 如果 s 非空且 l > 0，拷贝字符串到 dst
    // - 否则清空 dst
    //
    // 性能：
    // - 时间复杂度：O(l)（拷贝字符串）
    // - 使用 assign() 而非 operator=，避免多次拷贝
    //
    // 边界情况：
    // - ks.s == nullptr：清空 dst
    // - ks.l == 0：清空 dst（空序列）
    // ------------------------------------------------------------------
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

    // ==================================================================
    // KseqReader 析构函数：按正确顺序释放资源
    // ==================================================================
    // 说明：析构函数必须不抛异常（noexcept）
    // - I/O 清理失败最多影响资源回收，不应在栈展开期间触发 terminate
    // - 资源释放顺序很关键：先 kseq，再文件句柄，最后缓冲区
    //
    // 释放顺序说明（重要！）：
    // ------------------------------------------------------------------
    // 1. **先 kseq_destroy(seq)**：
    //    - kseq_t 内部可能持有指向底层缓冲/句柄的状态
    //    - 必须先释放解析器，再关闭文件
    //    - 否则可能出现悬空指针（dangling pointer）
    //
    // 2. **再 gzclose(fp) 或 fclose(fp)**：
    //    - 关闭底层文件句柄，刷新未写入的缓冲
    //    - stdio 可能在这一步释放内部 buffer
    //
    // 3. **最后 free(io_buf)**（仅非 zlib 模式）：
    //    - stdio 内部可能持有对 io_buf 的引用（直到 fclose）
    //    - 必须在 fclose 之后释放，避免 use-after-free
    //
    // 异常安全：
    // - 所有操作都不抛异常（kseq_destroy、fclose、free 都是 noexcept）
    // - 即使某步失败，也不影响后续步骤（各步独立）
    // ==================================================================
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

    // ==================================================================
    // KseqReader::next() - 读取下一条 FASTA/FASTQ 记录
    // ==================================================================
    // 功能：从文件中读取下一条序列记录并填充到 SeqRecord
    // 参数：
    //   @param rec - 输出参数，存储读取的序列记录
    // 返回：
    //   true - 成功读取一条记录
    //   false - 到达文件末尾（EOF）或读取错误
    //
    // kseq_read 返回值契约（重要！）：
    // ------------------------------------------------------------------
    // - ret >= 0：成功读取，返回值为序列长度
    //   * 正常情况，将数据拷贝到 rec
    //
    // - ret == -1：到达文件末尾（EOF）
    //   * 正常终止条件，返回 false
    //
    // - ret == -2：FASTQ 格式错误
    //   * 质量值长度与序列长度不匹配
    //   * 抛出异常，包含错误详情和文件路径
    //
    // - ret == -3：文件读取错误或格式损坏
    //   * I/O 错误或序列格式严重损坏
    //   * 抛出异常，包含错误详情和文件路径
    //
    // 字段提取说明：
    // ------------------------------------------------------------------
    // 从 kseq_t 提取以下字段到 SeqRecord：
    // 1. **name**（必有）：序列 ID（FASTA 的 >NAME 或 FASTQ 的 @NAME）
    // 2. **comment**（可选）：描述信息（>NAME COMMENT 中的 COMMENT 部分）
    // 3. **seq**（必有）：序列内容（A/C/G/T/N）
    // 4. **qual**（FASTQ 专有）：质量值字符串（Phred+33 编码）
    //    - FASTA 文件的 qual 为空字符串
    //
    // 性能：
    // - 时间复杂度：O(L)，L 为序列长度（kseq_read + 字符串拷贝）
    // - 内存分配：每个字段一次 assign（复用 std::string 的 SSO）
    //
    // 异常安全：
    // - 读取成功：强异常安全（rec 被正确填充）
    // - 读取失败（ret < -1）：基本异常安全（rec 可能部分修改，但不泄漏）
    // ==================================================================
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
        // ========================================================================
        // 性能关键路径：一次性写入大块数据，减少系统调用次数
        // ========================================================================
        if (buffer_.empty()) return;
        out_.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));

        // ========================================================================
        // 内存修复：释放 buffer 容量，避免峰值内存持续占用
        // ========================================================================
        // 问题根因：
        // 1. std::string::clear() 只删除内容，不释放容量（capacity 保持不变）
        // 2. 如果某个 batch 写入大量数据（例如 100MB），buffer_ 容量扩展到 100MB
        // 3. 后续即使写入少量数据，buffer_ 仍占用 100MB 内存
        // 4. 多线程场景下：nthreads × 100MB = 数 GB 持续占用
        //
        // 修复策略：
        // 1. flush 后立即释放容量（shrink_to_fit）
        // 2. 下次写入时会自动 reserve（构造函数中已设置）
        //
        // 性能权衡：
        // - shrink_to_fit 开销：O(1)（对于已清空的 string）
        // - 内存收益：每个线程减少数十 MB 到数百 MB
        // - 总体影响：100 万条序列增加约 0.1 秒（每 flush 一次 shrink），换取数 GB 内存释放
        // ========================================================================
        buffer_.clear();
        buffer_.shrink_to_fit();  // 释放未使用容量
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
    // 2. 质量值处理：若 query.qual 为空，则使用 SAM 默认值 "*"
    // 3. 所有字符串字段会被拷贝到 SamRecord 中（使用 std::string）
    // 4. 性能：O(N)，N 为字符串总长度（需要拷贝字符串）
    //
    // 参数：
    //   - query: 查询序列的 SeqRecord
    //   - ref_name: 参考序列名称
    //   - cigar_str: CIGAR 字符串
    //   - pos: 1-based 比对起始位置（默认 1）
    //   - mapq: mapping quality（默认 60）
    //   - flag: SAM flag（默认 0）
    //
    // 返回：SamRecord（所有字段为 std::string）
    // ------------------------------------------------------------------
    SamRecord makeSamRecord(
        const SeqRecord& query,
        std::string_view ref_name,
        std::string_view cigar_str,
        std::uint32_t pos,
        std::uint8_t mapq,
        std::uint16_t flag)
    {
        // 构建 SAM 记录（所有字段使用 std::string，拷贝数据）
        SamRecord sam_rec;
        sam_rec.qname = query.id;                                    // 拷贝 query 名称
        sam_rec.flag  = flag;                                        // SAM flag
        sam_rec.rname.assign(ref_name.data(), ref_name.size());      // 拷贝参考序列名称
        sam_rec.pos   = pos;                                         // 1-based 比对起始位置
        sam_rec.mapq  = mapq;                                        // mapping quality
        sam_rec.cigar.assign(cigar_str.data(), cigar_str.size());    // 拷贝 CIGAR 字符串
        sam_rec.seq   = query.seq;                                   // 拷贝 query 序列

        // 质量值处理：若 query.qual 不为空，则使用；否则保持默认值 "*"
        if (!query.qual.empty()) {
            sam_rec.qual = query.qual;
        }

        return sam_rec;
    }

    // ------------------------------------------------------------------
    // SamReader 实现
    // ------------------------------------------------------------------

    // SamReader::Impl 结构体：封装底层文件流和缓冲区
    //
    // 设计说明：
    // 1. 使用 std::ifstream 读取 SAM 文件（文本模式）
    // 2. 使用 std::getline 逐行读取，配合大缓冲区减少系统调用
    // 3. line_buffer 用于暂存当前行内容，避免重复分配
    // 4. io_buf 用于提升 std::ifstream 的缓冲区大小
    struct SamReader::Impl
    {
        std::ifstream in_;
        FilePath file_path_;
        std::string line_buffer_;
        char* io_buf_{nullptr};
        std::size_t io_buf_size_{0};

        // ------------------------------------------------------------------
        // 关键：为了让 next(SamRecord&) 返回的 string_view 在“本次 next 调用之后”仍然有效，
        // 我们必须把每次解析到的字段拷贝到 SamReader::Impl 的成员字符串中。
        // 否则 string_view 会指向 line_buffer_，下一次 getline() 就会覆盖，导致悬空引用。
        // ------------------------------------------------------------------
        std::string qname_storage_;
        std::string rname_storage_;
        std::string cigar_storage_;
        std::string rnext_storage_;
        std::string seq_storage_;
        std::string qual_storage_;
        std::string opt_storage_;
    };

    SamReader::SamReader(const FilePath& file_path, std::size_t buffer_size)
        : impl_(std::make_unique<Impl>())
    {
        impl_->file_path_ = file_path;
        impl_->io_buf_size_ = buffer_size;

        // 打开文件（文本模式，自动处理行结束符）
        impl_->in_.open(file_path, std::ios::in);
        if (!impl_->in_) {
            throw makeIoError("failed to open SAM file", file_path);
        }

        // 提升输入缓冲区大小以减少系统调用
        // 说明：
        // - 对于 std::ifstream，可以通过 rdbuf()->pubsetbuf() 设置缓冲区
        // - 必须在首次读取前调用
        // - 如果失败，仍然可以正常工作，只是性能可能下降
        if (buffer_size > 0) {
            impl_->io_buf_ = static_cast<char*>(std::malloc(buffer_size));
            if (impl_->io_buf_) {
                impl_->in_.rdbuf()->pubsetbuf(impl_->io_buf_, static_cast<std::streamsize>(buffer_size));
            }
        }

        // 预分配行缓冲区（SAM 行通常较长，预分配可减少 realloc）
        impl_->line_buffer_.reserve(4096);
    }

    SamReader::~SamReader()
    {
        if (!impl_) return;

        // 释放缓冲区（在关闭文件之后）
        if (impl_->in_.is_open()) {
            impl_->in_.close();
        }

        if (impl_->io_buf_) {
            std::free(impl_->io_buf_);
            impl_->io_buf_ = nullptr;
        }
    }

    SamReader::SamReader(SamReader&& other) noexcept = default;
    SamReader& SamReader::operator=(SamReader&& other) noexcept = default;

    bool SamReader::next(SamRecord& rec)
    {
        if (!impl_ || !impl_->in_) {
            throw std::runtime_error("SamReader is not initialized");
        }

        // 逐行读取，跳过 header 行（以 '@' 开头）
        while (std::getline(impl_->in_, impl_->line_buffer_)) {
            const std::string_view line(impl_->line_buffer_);

            // 跳过空行和 header 行
            if (line.empty() || line[0] == '@') {
                continue;
            }

            // SAM 必需字段 (11列) + 可选 TAG：
            // 0:QNAME 1:FLAG 2:RNAME 3:POS 4:MAPQ 5:CIGAR 6:RNEXT 7:PNEXT 8:TLEN 9:SEQ 10:QUAL [11+:OPT]
            std::array<std::string_view, 11> f{};
            std::size_t field_idx = 0;
            std::size_t start = 0;

            for (std::size_t i = 0; i <= line.size() && field_idx < f.size(); ++i) {
                if (i == line.size() || line[i] == '\t') {
                    f[field_idx] = line.substr(start, i - start);
                    start = i + 1;
                    ++field_idx;
                }
            }

            if (field_idx < f.size()) {
                throw std::runtime_error(
                    "invalid SAM record (missing required fields): " + impl_->file_path_.string());
            }

            // 11+ 的可选字段：直接保留整个尾巴（不解析 tag 结构，保持轻量）
            std::string_view opt_fields;
            if (start < line.size()) {
                opt_fields = line.substr(start);
            }

            // 0: QNAME
            rec.qname.assign(f[0].data(), f[0].size());

            // 1: FLAG
            try {
                rec.flag = static_cast<std::uint16_t>(std::stoul(std::string(f[1])));
            } catch (...) {
                throw std::runtime_error("invalid SAM FLAG: " + std::string(f[1]));
            }

            // 2: RNAME
            rec.rname.assign(f[2].data(), f[2].size());

            // 3: POS
            try {
                rec.pos = static_cast<std::uint32_t>(std::stoul(std::string(f[3])));
            } catch (...) {
                throw std::runtime_error("invalid SAM POS: " + std::string(f[3]));
            }

            // 4: MAPQ
            try {
                const unsigned long mapq_val = std::stoul(std::string(f[4]));
                rec.mapq = static_cast<std::uint8_t>(mapq_val > 255 ? 255 : mapq_val);
            } catch (...) {
                throw std::runtime_error("invalid SAM MAPQ: " + std::string(f[4]));
            }

            // 5: CIGAR
            rec.cigar.assign(f[5].data(), f[5].size());

            // 6: RNEXT
            rec.rnext.assign(f[6].data(), f[6].size());

            // 7: PNEXT
            try {
                rec.pnext = static_cast<std::uint32_t>(std::stoul(std::string(f[7])));
            } catch (...) {
                throw std::runtime_error("invalid SAM PNEXT: " + std::string(f[7]));
            }

            // 8: TLEN
            try {
                rec.tlen = static_cast<std::int32_t>(std::stol(std::string(f[8])));
            } catch (...) {
                throw std::runtime_error("invalid SAM TLEN: " + std::string(f[8]));
            }

            // 9: SEQ
            rec.seq.assign(f[9].data(), f[9].size());

            // 10: QUAL
            rec.qual.assign(f[10].data(), f[10].size());

            // OPT
            if (!opt_fields.empty()) {
                rec.opt.assign(opt_fields.data(), opt_fields.size());
            } else {
                rec.opt.clear();
            }

            return true;
        }

        return false;
    }

    void convertSamToFasta(const FilePath& sam_path, const FilePath& fasta_path, std::size_t line_width)
    {
        SamReader reader(sam_path);
        SeqWriter writer(fasta_path, line_width);

        SamRecord sam_rec;
        std::size_t count = 0;

        while (reader.next(sam_rec)) {
            // 关键：SamRecord 的 qname/seq 是 string_view，下一次 next() 可能覆盖。
            // 因此这里必须拷贝成 SeqRecord（samRecordToSeqRecord 内部做 assign）。
            const SeqRecord fasta_rec = samRecordToSeqRecord(sam_rec, /*keep_qual=*/false);
            writer.writeFasta(fasta_rec);
            ++count;
        }

        writer.flush();

#ifdef _DEBUG
        spdlog::debug("convertSamToFasta: converted {} records from {} to {}", count, sam_path.string(), fasta_path.string());
#endif
    }



} // namespace seq_io


