#ifndef HALIGN4_UTILS_H
#define HALIGN4_UTILS_H

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
    // Declarations for file IO utilities. Implementations moved to src/utils/file_io.cpp.

    std::string formatFsError(std::string_view msg,
                                     const FilePath& p,
                                     const std::error_code& ec);

    void requireExists(const FilePath& p, std::string_view what);

    void requireRegularFile(const FilePath& p, std::string_view what);

    void requireDirectory(const FilePath& p, std::string_view what);

    void ensureDirectoryExists(const FilePath& p, std::string_view what = "directory");

    bool isEmpty(const FilePath& p);

    // 核心：准备工作目录（不存在则创建；存在则要求为目录；可选要求为空）
    void prepareEmptydir(const FilePath& workdir, bool must_be_empty);

    // 可选：确保输出文件的父目录存在
    void ensureParentDirExists(const FilePath& out_file);

    // 判断给定的 FilePath 是否表示一个远程 URL（例如 http://、https://、ftp://、s3:// 等）
    // 返回 true 表示应当通过网络下载（远程链接）；返回 false 表示为本地路径。
    bool isUrl(const FilePath& p);

    // copy file with overwrite; falls back to stream copy on cross-device errors
    void copyFile(const FilePath& src, const FilePath& dst);

    // 下载远程文件到本地文件（优先使用 libcurl；无 libcurl 时回退到命令行工具）
    void downloadFile(const std::string& url, const FilePath& dst);

    // 下载或复制：如果 srcOrUrl 表示远程 URL（由 isUrl 判断），则调用 downloadFile；否则调用 copyFile。
    void fetchFile(const FilePath& srcOrUrl, const FilePath& dst);

    // 可选：危险操作（如你未来需要清空 workdir）
    void removeAll(const FilePath& p);

    // 获取文件内容
    std::string readFileToString(const FilePath& p);

}


namespace seq_io
{
#if defined(FilePath)
    using FilePath = ::FilePath;
#else
    using FilePath = std::filesystem::path;
#endif

    struct SeqRecord
    {
        std::string id;     // header name（到空格前）
        std::string desc;   // header 其余部分
        std::string seq;    // 序列
        std::string qual;   // FASTQ 质量字符串（可选），实现中需要访问该字段
    };

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

    class ISequenceReader
    {
    public:
        virtual ~ISequenceReader() = default;
        virtual bool next(SeqRecord& rec) = 0; // false => EOF
    };

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

