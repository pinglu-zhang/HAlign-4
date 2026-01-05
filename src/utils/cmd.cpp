#include "utils.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <sys/wait.h> // WIFEXITED/WEXITSTATUS/WIFSIGNALED/WTERMSIG

#include "config.hpp"

// 不使用 namespace detail：仅在本 cpp 内部可见的 helper（internal linkage）
// replaceAll: 将字符串中所有出现的子串替换为另一个字符串。
// 注意：若 from 为空则不做任何事情；该实现会修改原字符串并支持被替换的目标长度变化。
static void replaceAll(std::string& s, const std::string& from, const std::string& to)
{
    if (from.empty()) return;
    std::size_t pos = 0;
    while ((pos = s.find(from, pos)) != std::string::npos) {
        s.replace(pos, from.size(), to);
        pos += to.size();
    }
}

// containsToken: 检查一个字符串是否包含指定的子串，主要用于快速判断模板里是否含有占位符
// 轻量工具，未做完整的模板解析，仅作存在性检测。
static bool containsToken(const std::string& s, const std::string& tok)
{
    return s.find(tok) != std::string::npos;
}

namespace cmd
{
    // buildCommand:
    // - cmd_template: 用户提供的命令模板，必须至少包含 {input} 和 {output} 占位符。
    // - input_path/output_path: 会替换模板中的 {input}/{output}。
    // - thread: 当模板含有 {thread} 且 thread >= 0 时会替换为线程数；传入 -1 表示不替换 {thread}。
    // - opt: BuildOptions 控制静默（quiet）与是否关闭 stdin（close_stdin）等策略。
    // 返回值：最终用于 system() 的命令行字符串（注意：是 shell 命令行，可能包含重定向）。
    // 安全提示：本函数简单做字符串替换，不对路径做 shell 转义；如果命令模板来自不可信来源，
    // 需用适当的转义或避免通过 shell 执行（如使用 execv 系列直接传参）。
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

        // 如果用户要求静默或关闭 stdin，则需要在命令后拼接重定向符号。
        // note: detect_stdout_redirect 只是检测是否有 '>' 字符，属于简单启发式判断，不完全可靠。
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

    // runCommand:
    // - 直接使用 std::system 执行 shell 命令并返回命令的 exit code（非布尔）
    // - 返回语义：
    //     * -1: system() 本身调用失败（如 fork 失败）；
    //     * 若 WIFEXITED(status) 为真，则返回子进程的退出码（WEXITSTATUS(status))；
    //     * 否则返回原始 status（例如被信号终止的编码）。
    // 注意：std::system 会调用 /bin/sh -c "command"，存在 shell 注入风险，命令参数应由调用者妥善转义。
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

    // testCommandTemplate:
    // - 在给定的 workdir 下创建临时目录 WORKDIR_TMP，并在其中写入一个小型的 FASTA（tiny.fasta）作为输入，
    //   随后构造并运行用户提供的命令模板，检查命令是否成功（退出码 0）并且输出文件是否生成。
    // - 参数：
    //     * cmd_template: 用户命令模板（会被 buildCommand 替换 {input}/{output}）
    //     * workdir: 基础工作目录（函数会在 workdir/WORKDIR_TMP 下创建临时文件）
    //     * thread: 传给 buildCommand 的线程替换值
    // - 返回：bool，true 表示命令成功执行并产生输出文件；false 表示任一步骤失败。
    // - 清理：函数最后会尝试删除临时目录（即使命令失败也会尝试清理），删除失败只会记录警告。
    // 安全与健壮性注意事项：
    // - 该函数用于“自检”模板是否能正常运行，命令执行有副作用（会在本地文件系统写入），请在可信环境运行。
    // - 对于长时间运行或交互式命令，建议在模板中加入超时或使用更安全的运行方式。
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
