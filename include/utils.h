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
    inline std::string formatFsError(std::string_view msg,
                                     const FilePath& p,
                                     const std::error_code& ec) {
        std::string s;
        s.reserve(msg.size() + p.string().size() + 64);
        s.append(msg);
        s.append(": ");
        s.append(p.string());
        if (ec) {
            s.append(" (");
            s.append(ec.message());
            s.append(")");
        }
        return s;
    }

    inline void requireExists(const FilePath& p, std::string_view what) {
        std::error_code ec;
        const bool ok = fs::exists(p, ec); // non-throwing overload clears/sets ec accordingly :contentReference[oaicite:1]{index=1}
        if (ec || !ok) {
            throw std::runtime_error(formatFsError(std::string(what) + " does not exist", p, ec));
        }
    }

    inline void requireRegularFile(const FilePath& p, std::string_view what) {
        std::error_code ec;
        const bool ok = fs::is_regular_file(p, ec); // non-throwing overload :contentReference[oaicite:2]{index=2}
        if (ec || !ok) {
            throw std::runtime_error(formatFsError(std::string(what) + " is not a regular file", p, ec));
        }
    }

    inline void requireDirectory(const FilePath& p, std::string_view what) {
        std::error_code ec;
        const bool ok = fs::is_directory(p, ec); // non-throwing overload :contentReference[oaicite:3]{index=3}
        if (ec || !ok) {
            throw std::runtime_error(formatFsError(std::string(what) + " is not a directory", p, ec));
        }
    }

    inline void ensureDirectoryExists(const FilePath& p, std::string_view what = "directory") {
        std::error_code ec;
        if (fs::exists(p, ec)) {
            if (ec) {
                throw std::runtime_error(formatFsError(std::string(what) + " status check failed", p, ec));
            }
            requireDirectory(p, what);
            return;
        }
        if (ec) {
            throw std::runtime_error(formatFsError(std::string(what) + " status check failed", p, ec));
        }

        // create_directories has error_code overload (non-throwing) :contentReference[oaicite:4]{index=4}
        fs::create_directories(p, ec);
        if (ec) {
            throw std::runtime_error(formatFsError(std::string("failed to create ") + std::string(what), p, ec));
        }
    }

    inline bool isEmpty(const FilePath& p) {
        std::error_code ec;
        const bool empty = fs::is_empty(p, ec); // checks empty file or directory :contentReference[oaicite:5]{index=5}
        if (ec) {
            throw std::runtime_error(formatFsError("failed to check emptiness", p, ec));
        }
        return empty;
    }

    // 核心：准备工作目录（不存在则创建；存在则要求为目录；可选要求为空）
    inline void prepareEmptydir(const FilePath& workdir, bool must_be_empty) {
        if (workdir.empty()) {
            throw std::runtime_error("workdir is empty");
        }

        ensureDirectoryExists(workdir, "workdir");

        if (must_be_empty && !isEmpty(workdir)) {
            throw std::runtime_error("workdir must be empty: " + workdir.string());
        }
    }

    // 可选：确保输出文件的父目录存在
    inline void ensureParentDirExists(const FilePath& out_file) {
        if (out_file.empty()) return;
        auto parent = out_file.parent_path();
        if (parent.empty()) return;
        ensureDirectoryExists(parent, "output parent dir");
    }

    // 判断给定的 FilePath 是否表示一个远程 URL（例如 http://、https://、ftp://、s3:// 等）
    // 返回 true 表示应当通过网络下载（远程链接）；返回 false 表示为本地路径。
    // 实现说明：使用正则匹配 RFC 风格的 scheme:// 前缀，另外也识别 protocol-relative 地址（以 // 开头）。
    inline bool isUrl(const FilePath& p) {
        const std::string s = p.string();
        if (s.size() >= 2 && s[0] == '/' && s[1] == '/') return true; // 协议相对 URL：//host/path
        static const std::regex scheme_re(R"(^[a-zA-Z][a-zA-Z0-9+.\-]*://)");
        return std::regex_search(s, scheme_re);
    }

    // copy file with overwrite; falls back to stream copy on cross-device errors
    // 中文注释：
    // 将 src 复制到 dst，行为说明：
    // - 要求 src 必须是常规文件（非目录或特殊文件），否则抛出异常。
    // - 确保 dst 的父目录存在（不存在则创建）。
    // - 优先使用 std::filesystem::copy_file 的快速实现，并使用 overwrite_existing 覆盖已存在的目标文件。
    // - 如果 copy_file 因为跨设备（EXDEV）而失败，则回退到基于流的拷贝（逐字节读取/写入）。
    // - 其他任何文件系统错误都会抛出 std::runtime_error，错误消息通过 formatFsError 生成，包含路径和错误信息。
    // - 函数为 inline，以便在头文件中使用。
    inline void copyFile(const FilePath& src, const FilePath& dst) {
        // 验证源是普通文件
        requireRegularFile(src, "source file");
        // 确保目标父目录存在
        ensureParentDirExists(dst);

        std::error_code ec;
        // 尝试快速复制（覆盖已存在的目标）
        if (fs::copy_file(src, dst, fs::copy_options::overwrite_existing, ec)) {
            // 成功返回
            return;
        }

        // copy_file 可能返回 false 且不设置 ec；如果此时 dst 已存在，则视为成功
        if (!ec) {
            std::error_code exists_ec;
            if (fs::exists(dst, exists_ec) && !exists_ec) return;
            throw std::runtime_error(formatFsError("failed to copy file (unknown reason)", dst, exists_ec));
        }

        // 跨设备（EXDEV）错误时的回退：使用流拷贝，逐字节复制文件内容
        if (ec == std::make_error_code(std::errc::cross_device_link)) {
            std::ifstream in(src, std::ios::binary);
            if (!in) {
                throw std::runtime_error(formatFsError("failed to open source for reading", src, std::make_error_code(std::errc::io_error)));
            }
            std::ofstream out(dst, std::ios::binary | std::ios::trunc);
            if (!out) {
                throw std::runtime_error(formatFsError("failed to open destination for writing", dst, std::make_error_code(std::errc::io_error)));
            }
            out << in.rdbuf();
            if (!out) {
                throw std::runtime_error(formatFsError("failed to write file", dst, std::make_error_code(std::errc::io_error)));
            }
            return;
        }

        // 其他错误统一抛出异常，包含错误信息
        throw std::runtime_error(formatFsError("failed to copy file", dst, ec));
    }

    // 下载远程文件到本地文件（优先使用 libcurl；无 libcurl 时回退到命令行工具）
    // 中文注释：
    // - 参数：url 为远程资源链接，dst 为要写入的本地路径（FilePath）。
    // - 行为：确保 dst 的父目录存在，下载数据并写入 dst（覆盖已存在的文件）。
    // - 实现：如果可用，使用 libcurl 的 easy 接口并直接写入 std::ofstream；否则尝试运行外部命令 curl 或 wget。
    // - 错误处理：在任何失败情况下抛出 std::runtime_error，消息包含有用的信息。
    inline void downloadFile(const std::string& url, const FilePath& dst) {
        if (url.empty()) {
            throw std::runtime_error("download url is empty");
        }
        ensureParentDirExists(dst);

#if defined(HALIGN4_HAVE_LIBCURL)
        // libcurl 写回调：将收到的数据写入传入的 std::ofstream
        auto write_cb = [](char* ptr, size_t size, size_t nmemb, void* userdata) -> size_t {
            auto out = static_cast<std::ofstream*>(userdata);
            out->write(ptr, static_cast<std::streamsize>(size * nmemb));
            return out->good() ? size * nmemb : 0;
        };

        std::ofstream out(dst, std::ios::binary | std::ios::trunc);
        if (!out) {
            throw std::runtime_error("failed to open destination for writing: " + dst.string());
        }

        CURL* curl = curl_easy_init();
        if (!curl) {
            throw std::runtime_error("failed to initialize libcurl");
        }
        CURLcode res;
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
        curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L); // treat HTTP >=400 as error
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, +write_cb);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &out);
        // perform
        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);
        if (res != CURLE_OK) {
            out.close();
            // 删除可能不完整的目标文件
            std::error_code rem_ec;
            fs::remove(dst, rem_ec);
            throw std::runtime_error(std::string("libcurl error: ") + curl_easy_strerror(res));
        }
        out.close();
#else
        // 回退：尝试使用系统上的 curl 或 wget
        // 优先使用 curl：-L 跟随重定向，-s 静默，-S 在出错时显示，-f 在 HTTP 错误时返回非零
        std::string cmd = "curl -L -sSf -o \"" + dst.string() + "\" \"" + url + "\"";
        int rc = std::system(cmd.c_str());
        if (rc == 0) return;
        // 如果 curl 不存在或失败，尝试 wget（较保守的静默模式）
        cmd = "wget -q -O \"" + dst.string() + "\" \"" + url + "\"";
        rc = std::system(cmd.c_str());
        if (rc == 0) return;

        // 两种方式都失败，尝试删除可能的临时/部分文件
        std::error_code rem_ec;
        fs::remove(dst, rem_ec);
        throw std::runtime_error("failed to download file using curl/wget (rc=" + std::to_string(rc) + ")");
#endif
    }

    // 下载或复制：如果 srcOrUrl 表示远程 URL（由 isUrl 判断），则调用 downloadFile；否则调用 copyFile。
    // 中文注释：
    // - 参数：srcOrUrl 可为本地文件路径或 URL（例如 "https://..."）；dst 为目标本地路径。
    // - 语义：对用户透明地将远程资源下载到本地，或将本地文件复制到目标位置。
    // - 错误处理：底层的 downloadFile/copyFile 会在失败时抛出异常，本函数不会吞掉异常。
    inline void fetchFile(const FilePath& srcOrUrl, const FilePath& dst) {
        if (isUrl(srcOrUrl)) {
            // 对远程 URL，使用 downloadFile（downloadFile 接受 string URL）
            downloadFile(srcOrUrl.string(), dst);
        } else {
            // 本地文件，直接复制（copyFile 会验证源为常规文件并处理跨设备等情况）
            copyFile(srcOrUrl, dst);
        }
    }

    // 可选：危险操作（如你未来需要清空 workdir）
    // remove_all 的语义与返回值见标准库文档 :contentReference[oaicite:6]{index=6}
    inline void removeAll(const FilePath& p) {
        std::error_code ec;
        fs::remove_all(p, ec);
        if (ec) {
            throw std::runtime_error(formatFsError("remove_all failed", p, ec));
        }
    }

    // 获取文件内容
    inline std::string readFileToString(const FilePath& p) {
        requireRegularFile(p, "input file");

        std::ifstream in(p, std::ios::binary);
        if (!in) {
            throw std::runtime_error("failed to open file for reading: " + p.string());
        }

        std::string content;
        in.seekg(0, std::ios::end);
        content.resize(static_cast<std::size_t>(in.tellg()));
        in.seekg(0, std::ios::beg);
        in.read(&content[0], static_cast<std::streamsize>(content.size()));
        if (!in) {
            throw std::runtime_error("failed to read file: " + p.string());
        }
        return content;
    }


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

