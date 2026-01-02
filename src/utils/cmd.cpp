#include "utils.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <sys/wait.h> // WIFEXITED/WEXITSTATUS/WIFSIGNALED/WTERMSIG

#include "config.hpp"

// 不使用 namespace detail：仅在本 cpp 内部可见的 helper（internal linkage）
static void replaceAll(std::string& s, const std::string& from, const std::string& to)
{
    if (from.empty()) return;
    std::size_t pos = 0;
    while ((pos = s.find(from, pos)) != std::string::npos) {
        s.replace(pos, from.size(), to);
        pos += to.size();
    }
}

static bool containsToken(const std::string& s, const std::string& tok)
{
    return s.find(tok) != std::string::npos;
}

namespace cmd
{
    std::string buildCommand(std::string cmd_template,
                             const std::string& input_path,
                             const std::string& output_path,
                             int thread,
                             const BuildOptions& opt)
    {
        // 必须包含 {input} 和 {output}
        if (!containsToken(cmd_template, "{input}")) {
            throw std::runtime_error("cmd template missing {input}");
        }
        if (!containsToken(cmd_template, "{output}")) {
            throw std::runtime_error("cmd template missing {output}");
        }

        // 占位符替换
        replaceAll(cmd_template, "{input}", input_path);
        replaceAll(cmd_template, "{output}", output_path);

        if (containsToken(cmd_template, "{thread}") && thread >= 0) {
            replaceAll(cmd_template, "{thread}", std::to_string(thread));
        }

        if (!opt.quiet && !opt.close_stdin) {
            return cmd_template;
        }

        const bool has_stdout_redirect =
            opt.detect_stdout_redirect ? (cmd_template.find('>') != std::string::npos) : false;

        // 静默策略（Linux）：
        // - 若用户模板已做 stdout 重定向：只丢弃 stderr（2>/dev/null）
        // - 否则：stdout 丢弃 + stderr -> stdout（> /dev/null 2>&1）
        if (opt.quiet) {
            if (has_stdout_redirect) {
                cmd_template.append(" 2>/dev/null");
            } else {
                cmd_template.append(" > /dev/null 2>&1");
            }
        }

        // 关闭 stdin，避免外部命令等待交互输入
        if (opt.close_stdin) {
            cmd_template.append(" < /dev/null");
        }

        return cmd_template;
    }

    // 执行命令并返回命令是否成功退出（exit code == 0）。
    // 之所以返回 bool 是因为调用方只关心命令是否成功。
    // - 若 system() 失败（fork/exec）：返回 false
    // - 若命令正常退出且 exit code == 0：返回 true
    // - 其它情况（非零退出码或被信号终止）：返回 false
    int runCommand(const std::string& command)
    {
        const int status = std::system(command.c_str());

        if (status == -1) {
            return status;
        }
        if (WIFEXITED(status)) {
            return WEXITSTATUS(status);
        }
        // 如果被信号终止或其他不可识别的返回，视为失败
        return status;
    }

    // 使用 tiny.fasta 做自检，运行模板命令并检查输出文件是否生成
    bool testCommandTemplate(const std::string& cmd_template, const FilePath& workdir, int thread)
    {
        const FilePath temp_dir = workdir / WORKDIR_TMP;
        const FilePath in_path  = temp_dir / "tiny.fasta";
        const FilePath out_path = temp_dir / "aligned.fasta";

        bool ok = true;

        do {
            // 使用 file_io 确保临时目录存在（会创建目录或在错误时抛出）
            try {
                file_io::ensureDirectoryExists(temp_dir, WORKDIR_TMP);
            } catch (const std::exception& e) {
                spdlog::error("Failed to create {}: {}", temp_dir.string(), e.what());
                ok = false; break;
            }

            // 写一个很小的 FASTA
            {
                // 确保写入文件的父目录存在
                file_io::ensureParentDirExists(in_path);
                std::ofstream ofs(in_path);
                if (!ofs) {
                    spdlog::error("Cannot open {} for writing.", in_path.string());
                    ok = false; break;
                }
                ofs <<
                    R"(>seq1
ACGTACGTGA
>seq2
ACGTTGCA
>seq3
ACGTACGA
))";
             }

             std::string cmd_line;
            try {
                // 确保输出文件父目录存在（以防用户模板写到子目录）
                file_io::ensureParentDirExists(out_path);
                cmd_line = buildCommand(cmd_template, in_path.string(), out_path.string(), thread);
            } catch (const std::exception& e) {
                spdlog::error("buildCommand failed: {}", e.what());
                ok = false; break;
            }

            spdlog::info("Running: {}", cmd_line);
            const int rc = runCommand(cmd_line);

            if (rc != 0) {
                spdlog::error("cmd failed");
                ok = false;
                break;
            }

            // 使用 file_io::requireExists 来判断输出是否产生（会在失败时抛出）
            try {
                file_io::requireExists(out_path, "command output");
            } catch (const std::exception& e) {
                spdlog::error("Output file not found: {} ({})", out_path.string(), e.what());
                ok = false;
                break;
            }

            spdlog::info("cmd finished successfully, output: {}", out_path.string());

        } while (false);

        // 删除临时目录，使用 file_io::removeAll
        try {
            file_io::removeAll(temp_dir);
        } catch (const std::exception& e) {
            spdlog::warn("Failed to remove {}: {}", temp_dir.string(), e.what());
        }

        return ok;
    }

} // namespace cmd
