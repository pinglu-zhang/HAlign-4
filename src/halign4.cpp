#include <config.hpp>
#include <utils.h>
#include <thread>
#include "preprocess.h"
#include "consensus.h"
#include <chrono>

#include "align.h"

// 该文件是程序入口（main），负责：
// 1) 解析命令行参数（使用 CLI11），并将值绑定到 Options 结构；
// 2) 校验参数（例如文件存在性、工作目录准备、外部 MSA 命令模板测试）；
// 3) 调用预处理（preprocessInputFasta）将原始输入复制/下载到工作目录并清洗；
// 4) 使用 RefAligner 将查询序列比对到参考序列（中心序列/共识序列）；
// 5) 合并比对结果并输出最终的对齐文件。
//
// 注：本文件尽量把流程拆成职责清晰的小步骤，便于定位错误并在 CI 中逐步测试。

// checkOption：参数再校验与工作目录准备
// 该函数做了多项关键工作：
// - 使用 file_io 模块统一校验文件存在性与路径类型（统一错误信息与行为）；
// - 检查数值参数的合理性（threads/kmer_size/kmer_window/cons_n）；
// - 准备工作目录：Debug 模式允许非空 workdir 以便快速迭代；Release 模式要求 workdir 为空以避免覆盖已有数据；
// - 对 MSA 命令模板做一次自检（cmd::testCommandTemplate）：在临时目录运行模板并验证是否产生输出，用于尽早发现模板占位符/命令行问题。
//   注意：testCommandTemplate 会实际运行外部命令，可能依赖外部工具或运行环境；若不希望在参数校验阶段执行外部命令，可将测试延后并在错误时给出更明确提示。
static void checkOption(Options& opt) {
    // 文件相关：统一调用 file_io
    file_io::requireRegularFile(opt.input, "input");

    if (!opt.center_path.empty()) {
        file_io::requireRegularFile(opt.center_path, "center_path");
    }

    // 数值参数相关：仍在这里检查（避免后续出现未定义行为或异常分支更难定位）
    if (opt.threads <= 0) throw std::runtime_error("threads must be > 0");
    if (opt.kmer_size <= 0) throw std::runtime_error("kmer_size must be > 0");
    if (opt.kmer_window <= 0) throw std::runtime_error("kmer_window must be > 0");
    if (opt.cons_n <= 0) throw std::runtime_error("cons_n must be > 0");

    // k 与 w 的常见约束（贴近 minimap2 等常用实现的可实现范围）
    if (opt.kmer_size > 31) throw std::runtime_error("kmer_size too large (must be <= 31)");
    if (opt.kmer_window >= 256) {
        spdlog::warn("kmer_window >= 256 may be slow/unsupported by the fast minimizer path; current value: {}", opt.kmer_window);
    }

    // workdir：Debug 下允许非空，Release 下必须空（保持现有逻辑不变）
#ifdef _DEBUG
    constexpr bool must_be_empty = false;
#else
    constexpr bool must_be_empty = true;
#endif
    file_io::prepareEmptydir(opt.workdir, must_be_empty);

    std::string msa_cmd_str = DEFALT_MSA_CMD;
    if (!opt.msa_cmd.empty()) {
        file_io::requireRegularFile(opt.msa_cmd, "msa_cmd");
        msa_cmd_str = opt.msa_cmd;
    }

    if (cmd::testCommandTemplate(msa_cmd_str, opt.workdir, opt.threads)) {
        spdlog::info("msa_cmd template test passed.");
    } else {
        throw std::runtime_error("msa_cmd template test failed.");
    }
    opt.msa_cmd = msa_cmd_str;

    // 可选：确保输出父目录存在（如果你希望自动创建）
    // file_io::ensureParentDirExists(opt.output);
}

// cleanupWorkdir：根据 save_workdir 开关决定是否删除工作目录
// 说明：
// - 当 save_workdir=false（默认）时，程序在成功完成后会尝试删除 workdir，以避免堆积中间文件；
// - 当 save_workdir=true 时，保留 workdir 便于 debug/复现/检查中间产物；
// - 删除失败不会影响主流程返回值：这里用 warn 记录原因，避免“成功结果被清理失败掩盖”。
static void cleanupWorkdir(const Options& opt) {
    if (!opt.save_workdir) {
        try {
            spdlog::info("正在删除工作目录: {}", opt.workdir);
            file_io::removeAll(FilePath(opt.workdir));
            spdlog::info("工作目录删除成功");
        } catch (const std::exception& e) {
            spdlog::warn("工作目录删除失败: {}", e.what());
        }
    } else {
        spdlog::info("保留工作目录: {}", opt.workdir);
    }
}

int main(int argc, char** argv) {
    // main 的整体流程（高层概览）：
    // 1) 初始化日志与线程池
    // 2) 定义 CLI 并解析参数
    // 3) 打印并校验参数
    // 4) 将日志输出到工作目录
    // 5) 执行预处理（fetch/copy + 清洗 + Top-K 选择）
    // 6) 若用户指定 center_path，则拷贝到 workdir；否则使用预处理生成的共识序列
    // 7) 若输入序列总数 <= cons_n，则调用外部 MSA 对齐并直接输出退出（快速路径）
    // 8) 否则使用 RefAligner 比对（alignQueryToRef）并合并（mergeAlignedResults）
    // 9) 输出最终对齐结果
    // 10) 根据 --save-workdir 决定是否清理 workdir

    try
    {
        spdlog::init_thread_pool(8192, 1);
        setupLogger();

        Options opt;
        CLI::App app{"halign4"};


        setupCli(app, opt);

        // CLI11 的解析必须在 main 中进行，以便直接处理 argc/argv
        CLI11_PARSE(app, argc, argv);

        // 解析成功后，opt 已被填充
        logParsedOptions(opt);
        // 程序启动时打印版本号
        spdlog::info("Starting halign4 version {}...", VERSION);

        // 校验参数（文件存在性、工作目录、msa 模板自检等）
        checkOption(opt);

        // 将日志也写到 workdir 中（便于收集运行时日志）；此调用假定 workdir 已存在并已经被准备好
        setupLoggerWithFile(opt.workdir);

        // ---------------- 预处理阶段（IO 密集） ----------------
        // preprocessInputFasta 的职责：
        // - 如果 input 是 URL 则下载到 workdir/data/raw，否则拷贝到该目录；
        // - 逐条读取、清洗并写出到 data/clean，同时选出最长的 cons_n 条序列写入共识输入文件；
        // - 返回处理的总序列数（用于决定下一步是否需要进一步的合并）。
        uint_t preproc_count = preprocessInputFasta(opt.input, opt.workdir, opt.cons_n);
        spdlog::info("Preprocessing produced {} records", preproc_count);

        // ---------------- 共识对齐阶段（调用外部 MSA） ----------------
        // alignConsensusSequence：
        // - 执行外部 MSA 命令（由用户在 opt.msa_cmd 指定），将未对齐的共识序列文件对齐为 aligned 文件；
        // - 该函数会打印运行时间并对输出文件做基本检查（存在性与文件大小），但不会抛出异常以中断流程；
        FilePath consensus_unaligned_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_UNALIGNED;
        FilePath consensus_aligned_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_ALIGNED;
        FilePath consensus_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_FASTA;
        FilePath consensus_json_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_JSON;

        if (!opt.center_path.empty())
        {
            // 如果用户直接指定了中心序列路径，则将其拷贝到工作目录中。
            // 说明：
            // 1. 保持用户原始文件不变（不删除 opt.center_path）；
            // 2. 删除工作目录中可能存在的旧文件（consensus_unaligned_file）；
            // 3. 将用户指定的文件拷贝到工作目录，以便后续流程统一处理
            spdlog::info("Using user-specified center sequence: {}", opt.center_path);

            // 删除工作目录中旧文件（如果存在）
            if (std::filesystem::exists(consensus_unaligned_file)) {
                file_io::removeAll(consensus_unaligned_file);
            }

            // 从用户指定路径拷贝到工作目录
            file_io::copyFile(FilePath(opt.center_path), consensus_unaligned_file);
            spdlog::info("Center sequence copied to: {}", consensus_unaligned_file.string());
        }



        // 如果输入序列总数小于等于用于生成共识的数量（cons_n），说明我们已在预处理阶段处理完毕。
        // 此时直接把对齐输出拷贝到用户指定的最终输出文件并退出；这是常见的快速路径，避免无谓的合并工作。
        if (opt.center_path.empty() && preproc_count <= opt.cons_n)
        {
            // 调用外部 MSA （可能是 mafft/clustalo/其他高质量工具），会阻塞直到命令完成
            alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, opt.msa_cmd,  opt.threads);
            file_io::copyFile(consensus_aligned_file,FilePath(opt.output));
            spdlog::info("All sequences processed; final output written to {}", opt.output);

            // 清理工作目录后退出
            cleanupWorkdir(opt);

            spdlog::info("halign4 End!");
            return 0;
        }



        // ---------------- 比对阶段（使用 RefAligner） ----------------
        // 说明：
        // 1. 若用户指定了 center_path，则从该文件读取作为参考序列；
        // 2. 否则使用预处理生成的 consensus_file 作为参考序列；
        // 3. RefAligner 内部会持有参考序列并完成后续比对与合并流程。
        FilePath ref_path = opt.center_path.empty() ? consensus_file : FilePath(opt.center_path);

        align::RefAligner ref_aligner(opt, ref_path);
        ref_aligner.alignQueryToRef(opt.input);
        ref_aligner.mergeAlignedResults(consensus_aligned_file, opt.msa_cmd);
        // 调用 RefAligner 完成后续的比对与合并工作

        FilePath final_output_path = FilePath(opt.workdir) / RESULTS_DIR / FINAL_ALIGNED_FASTA;
        file_io::copyFile(final_output_path, FilePath(opt.output));
        spdlog::info("Final aligned output written to {}", opt.output);

        // 成功完成后，根据 --save-workdir 决定是否清理 workdir
        cleanupWorkdir(opt);

        spdlog::info("halign4 End!");

        return 0;
    } catch (const std::exception &e) {
        spdlog::error("Fatal error: {}", e.what());
        spdlog::error("halign4 End!");
        return 1;
    } catch (...) {
        spdlog::error("Fatal error: unknown exception");
        spdlog::error("halign4 End!");
        return 1;
    }
}
