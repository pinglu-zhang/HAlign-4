#include "../include/utils.h"

#include <system_error>
#include <iostream>

namespace file_io {

    // 格式化文件系统错误信息为可读字符串（带路径与错误消息）
    std::string formatFsError(std::string_view msg,
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

    // 确保路径存在，否则抛出异常（用于输入检验）
    void requireExists(const FilePath& p, std::string_view what) {
        std::error_code ec;
        const bool ok = fs::exists(p, ec);
        if (ec || !ok) {
            throw std::runtime_error(formatFsError(std::string(what) + " does not exist", p, ec));
        }
    }

    // 确保路径是常规文件，否则抛出异常
    void requireRegularFile(const FilePath& p, std::string_view what) {
        std::error_code ec;
        const bool ok = fs::is_regular_file(p, ec);
        if (ec || !ok) {
            throw std::runtime_error(formatFsError(std::string(what) + " is not a regular file", p, ec));
        }
    }

    // 确保路径是目录，否则抛出异常
    void requireDirectory(const FilePath& p, std::string_view what) {
        std::error_code ec;
        const bool ok = fs::is_directory(p, ec);
        if (ec || !ok) {
            throw std::runtime_error(formatFsError(std::string(what) + " is not a directory", p, ec));
        }
    }

    // 确保目录存在；若不存在则创建
    void ensureDirectoryExists(const FilePath& p, std::string_view what) {
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

        fs::create_directories(p, ec);
        if (ec) {
            throw std::runtime_error(formatFsError(std::string("failed to create ") + std::string(what), p, ec));
        }
    }

    // 检查路径是否为空（文件为空）
    bool isEmpty(const FilePath& p) {
        std::error_code ec;
        const bool empty = fs::is_empty(p, ec);
        if (ec) {
            throw std::runtime_error(formatFsError("failed to check emptiness", p, ec));
        }
        return empty;
    }

    // 准备一个空目录：若不存在则创建；若要求为空但非空则抛错
    void prepareEmptydir(const FilePath& workdir, bool must_be_empty) {
        if (workdir.empty()) {
            throw std::runtime_error("workdir is empty");
        }

        ensureDirectoryExists(workdir, "workdir");

        if (must_be_empty && !isEmpty(workdir)) {
            throw std::runtime_error("workdir must be empty: " + workdir.string());
        }
    }

    // 确保输出文件的父目录存在（写文件前调用）
    void ensureParentDirExists(const FilePath& out_file) {
        if (out_file.empty()) return;
        auto parent = out_file.parent_path();
        if (parent.empty()) return;
        ensureDirectoryExists(parent, "output parent dir");
    }

    // 判断给定路径是否是 URL（简单判断 scheme:// 或 // 开头）
    bool isUrl(const FilePath& p) {
        const std::string s = p.string();
        if (s.size() >= 2 && s[0] == '/' && s[1] == '/') return true;
        static const std::regex scheme_re(R"(^[a-zA-Z][a-zA-Z0-9+.\-]*://)");
        return std::regex_search(s, scheme_re);
    }

    // 复制文件；优先使用 std::filesystem::copy_file，遇到跨设备错误时回退到流拷贝
    void copyFile(const FilePath& src, const FilePath& dst) {
        requireRegularFile(src, "source file");
        ensureParentDirExists(dst);

        std::error_code ec;
        if (fs::copy_file(src, dst, fs::copy_options::overwrite_existing, ec)) {
            return;
        }

        if (!ec) {
            std::error_code exists_ec;
            if (fs::exists(dst, exists_ec) && !exists_ec) return;
            throw std::runtime_error(formatFsError("failed to copy file (unknown reason)", dst, exists_ec));
        }

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

        throw std::runtime_error(formatFsError("failed to copy file", dst, ec));
    }

    // 下载远程文件到本地：优先使用 libcurl（如果可用），否则回退到调用 curl/wget 命令行
    void downloadFile(const std::string& url, const FilePath& dst) {
        if (url.empty()) {
            throw std::runtime_error("download url is empty");
        }
        ensureParentDirExists(dst);

    #if defined(HALIGN4_HAVE_LIBCURL)
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
        curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, +write_cb);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &out);
        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);
        if (res != CURLE_OK) {
            out.close();
            std::error_code rem_ec;
            fs::remove(dst, rem_ec);
            throw std::runtime_error(std::string("libcurl error: ") + curl_easy_strerror(res));
        }
        out.close();
    #else
        std::string cmd = "curl -L -sSf -o \"" + dst.string() + "\" \"" + url + "\"";
        int rc = std::system(cmd.c_str());
        if (rc == 0) return;
        cmd = "wget -q -O \"" + dst.string() + "\" \"" + url + "\"";
        rc = std::system(cmd.c_str());
        if (rc == 0) return;

        std::error_code rem_ec;
        fs::remove(dst, rem_ec);
        throw std::runtime_error("failed to download file using curl/wget (rc=" + std::to_string(rc) + ")");
    #endif
    }

    // 如果 srcOrUrl 是 URL 则下载，否则复制本地文件到目标
    void fetchFile(const FilePath& srcOrUrl, const FilePath& dst) {
        if (isUrl(srcOrUrl)) {
            downloadFile(srcOrUrl.string(), dst);
        } else {
            copyFile(srcOrUrl, dst);
        }
    }

    // 递归删除路径（文件或目录），失败时抛出异常
    void removeAll(const FilePath& p) {
        std::error_code ec;
        fs::remove_all(p, ec);
        if (ec) {
            throw std::runtime_error(formatFsError("remove_all failed", p, ec));
        }
    }

    // 读取整个文件到一个 std::string 中（用于小文件）
    std::string readFileToString(const FilePath& p) {
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

} // namespace file_io

