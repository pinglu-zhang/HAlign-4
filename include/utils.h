#ifndef HALIGN4_UTILS_H
#define HALIGN4_UTILS_H

// ================================================================
// utils.h - 文件与序列 I/O 公共声明（详细中文注释）
//
// 本头文件包含两个主要功能块：
//  1) file_io: 一组与文件/路径/网络获取相关的工具函数（copy/download/fetch/prep目录等），
//     目的在于把所有与文件系统交互的逻辑集中管理，便于统一错误处理与测试。
//  2) seq_io: 一组基于 kseq 的 FASTA/FASTQ 读写抽象（KseqReader / FastaWriter / SeqRecord），
//     并提供序列清洗函数（cleanSequence）和便捷的 reader/writer 构造器。
// 此处只声明接口；具体实现位于 src/utils/file_io.cpp 和 src/utils/seq_io.cpp 中。
// ================================================================

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <regex>
#include <array>
#include <cstdint>

#if __has_include(<curl/curl.h>)
#include <curl/curl.h>
#define HALIGN4_HAVE_LIBCURL 1
#endif

namespace fs = std::filesystem;
using FilePath = std::filesystem::path;


namespace file_io
{
    // 下面是 file_io 中提供的函数声明，并在每个函数前添加中文说明，帮助使用者理解语义与边界。

    // 格式化文件系统错误信息为可读字符串（带路径与错误消息）
    // 用法：在抛异常或记录日志时调用，统一错误输出格式
    std::string formatFsError(std::string_view msg,
                                     const FilePath& p,
                                     const std::error_code& ec);

    // 要求路径存在（任何类型）；若不存在或检查失败则抛出 runtime_error
    void requireExists(const FilePath& p, std::string_view what);

    // 要求路径是常规文件；若不是则抛出 runtime_error
    void requireRegularFile(const FilePath& p, std::string_view what);

    // 要求路径是目录；若不是则抛出 runtime_error
    void requireDirectory(const FilePath& p, std::string_view what);

    // 确保目录存在；若不存在则创建（包括父目录）。出错抛异常。
    void ensureDirectoryExists(const FilePath& p, std::string_view what = "directory");

    // 判断路径是否为空（对文件适用）
    bool isEmpty(const FilePath& p);

    // 核心：准备工作目录（不存在则创建；存在则要求为目录；可选要求为空）
    // - workdir: 需准备的目录路径
    // - must_be_empty: 若 true，则目录存在时必须为空，否则抛错
    void prepareEmptydir(const FilePath& workdir, bool must_be_empty);

    // 可选：确保输出文件的父目录存在（写文件前调用）
    void ensureParentDirExists(const FilePath& out_file);

    // 判断给定的 FilePath 是否表示一个远程 URL（例如 http://、https://、ftp://、s3:// 等）
    // 返回 true 表示应当通过网络下载（远程链接）；返回 false 表示为本地路径。
    // 注意：该判断为启发式；对 file:// 或其他边缘 scheme，调用方应根据需要额外处理。
    bool isUrl(const FilePath& p);

    // copyFile：复制文件并覆盖目标；优先使用 std::filesystem::copy_file。
    // 若遇到跨设备错误（EXDEV），回退到基于流的拷贝以确保可跨分区复制。
    // 出错抛异常。
    void copyFile(const FilePath& src, const FilePath& dst);

    // downloadFile：把远程 URL 下载到本地路径 dst
    // - 如果项目编译时检测到 libcurl，则使用其 API（更可靠、可控）
    // - 否则回退到命令行工具（curl/wget），此方式依赖运行时环境
    // 出错抛异常并尝试删除局部残留文件。
    void downloadFile(const std::string& url, const FilePath& dst);

    // fetchFile：如果 srcOrUrl 是 URL，则调用 downloadFile，否则调用 copyFile
    // 这是上层调用统一接口，屏蔽本地/远程差异
    void fetchFile(const FilePath& srcOrUrl, const FilePath& dst);

    // 递归删除路径（文件或目录），失败时抛出异常。慎用！
    void removeAll(const FilePath& p);

    // 读取整个小文件到 std::string（仅适用于文件体积适中时）
    std::string readFileToString(const FilePath& p);

}


namespace seq_io
{
#if defined(FilePath)
    using FilePath = ::FilePath;
#else
    using FilePath = std::filesystem::path;
#endif

    // SeqRecord：用于在内存中表示一条序列记录（FASTA/FASTQ）
    // 字段语义：
    //  - id: header 名称，通常是 '>' 后第一个空格前的 token
    //  - desc: header 剩余描述信息（可选）
    //  - seq: 序列字符串（A/C/G/T/U/-/N 等）
    //  - qual: FASTQ 质量字符串（若为 FASTQ，可填充；FASTA 情况可留空）
    struct SeqRecord
    {
        std::string id;     // header name（到空格前）
        std::string desc;   // header 其余部分
        std::string seq;    // 序列
        std::string qual;   // FASTQ 质量字符串（可选），实现中需要访问该字段
    };

    // makeCleanTable / clean_table:
    // - 用于把任意字符映射到一个受限的碱基字符集合（'A','C','G','T','U','N','-'），
    // - 实现为 256 大小的查表（针对 unsigned char），避免运行时的分支/条件判断，极大提升批量清洗性能；
    // - cleanSequence 提供了 in-place 清洗接口，对大序列重复使用时能减少临时分配。
    constexpr std::array<std::uint8_t, 256> makeCleanTable()
    {
        std::array<std::uint8_t, 256> table{};

        for (std::size_t i = 0; i < table.size(); ++i) {
            table[i] = static_cast<std::uint8_t>('N');
        }

        table[static_cast<unsigned char>('A')] = 'A';
        table[static_cast<unsigned char>('a')] = 'A';
        table[static_cast<unsigned char>('C')] = 'C';
        table[static_cast<unsigned char>('c')] = 'C';
        table[static_cast<unsigned char>('G')] = 'G';
        table[static_cast<unsigned char>('g')] = 'G';
        table[static_cast<unsigned char>('T')] = 'T';
        table[static_cast<unsigned char>('t')] = 'T';
        table[static_cast<unsigned char>('U')] = 'U';
        table[static_cast<unsigned char>('u')] = 'U';
        table[static_cast<unsigned char>('-')] = '-';

        return table;
    }

    // 头文件内的“单一定义”常量表：inline constexpr 是推荐写法
    inline constexpr auto clean_table = makeCleanTable();
    inline void cleanSequence(std::string& seq)
    {
        for (char& ch : seq) {
            ch = static_cast<char>(clean_table[static_cast<unsigned char>(ch)]);
        }
    }
    inline void cleanSequence(SeqRecord& seq)
    {
        for (char& ch : seq.seq) {
            ch = static_cast<char>(clean_table[static_cast<unsigned char>(ch)]);
        }
    }

    // ISequenceReader：抽象读取器接口，支持 next(SeqRecord&) 返回 false 表示 EOF
    class ISequenceReader
    {
    public:
        virtual ~ISequenceReader() = default;
        virtual bool next(SeqRecord& rec) = 0; // false => EOF
    };

    // KseqReader：基于 kseq 实现的高效序列读取器
    // 说明：实现位于 src/utils/seq_io.cpp，构造函数会根据是否启用 zlib 选择 gzopen/gzread 或 fopen/fread
    class KseqReader final : public ISequenceReader
    {
    public:
        explicit KseqReader(const FilePath& file_path);
        ~KseqReader() override;

        KseqReader(const KseqReader&) = delete;
        KseqReader& operator=(const KseqReader&) = delete;

        KseqReader(KseqReader&&) noexcept;
        KseqReader& operator=(KseqReader&&) noexcept;

        bool next(SeqRecord& rec) override;

    private:
        struct Impl;
        std::unique_ptr<Impl> impl_;
    };

    // FastaWriter：用于把 SeqRecord 以 FASTA 格式写出到磁盘
    // 性能提示：write 会把 header 与序列按行宽组织为一个临时缓冲并一次性写入，
    // 避免逐字符或逐行的小量写调用，从而提升磁盘写吞吐。
    class FastaWriter
    {
    public:
        explicit FastaWriter(const FilePath& file_path, std::size_t line_width = 80);

        void write(const SeqRecord& rec);
        void flush();

    private:
        std::ofstream out_;
        std::size_t line_width_{80};

        static void writeWrapped(std::ofstream& out, std::string_view s, std::size_t width);
    };

    inline std::unique_ptr<ISequenceReader> openKseqReader(const FilePath& file_path)
    {
        return std::make_unique<KseqReader>(file_path);
    }

} // namespace seq_io

namespace cmd
{
    struct BuildOptions
    {
        bool quiet = true;                  // 是否追加静默重定向
        bool close_stdin = true;            // 是否追加 "< /dev/null"
        bool detect_stdout_redirect = true; // 是否根据 '>' 决定静默策略
    };

    // cmd_template 必须包含 {input} 和 {output}；可选包含 {thread}
    // Linux-only 规则：
    // - 若模板中检测到 '>'（stdout 已被重定向）：追加 "2>/dev/null"
    // - 否则追加 "> /dev/null 2>&1"
    // - 可选追加 "< /dev/null" 关闭 stdin
    std::string buildCommand(std::string cmd_template,
                             const std::string& input_path,
                             const std::string& output_path,
                             int thread,
                             const BuildOptions& opt = BuildOptions{});

    // 执行命令并返回真实退出码：
    // - 正常退出：返回 exit code
    // - 信号终止：返回 128 + signal
    // - system() 自身失败：返回 -1
    int runCommand(const std::string& command);

    // 用 tiny.fasta 做自检，运行模板命令并检查输出文件是否生成
    bool testCommandTemplate(const std::string& cmd_template, const FilePath& workdir, int thread = 1);

} // namespace cmd

#endif //HALIGN4_UTILS_H

