#include <config.hpp>
#include <utils.h>
#include <thread>
#include "preprocess.h"
#include "consensus.h"
#include <chrono>

#include "align.h"

// 该文件是程序的入口（main），负责：
// 1) 解析命令行参数（使用 CLI11），并把值绑定到 Options 结构；
// 2) 校验参数（例如文件存在性、工作目录准备、外部 MSA 命令模板测试）；
// 3) 调用预处理（preprocessInputFasta）把原始输入复制/下载到工作目录并清洗；
// 4) 执行针对共识序列的多序列比对（外部工具，通过命令模板），并产生对齐结果；
// 5) 以多线程并行/分批的方式生成最终共识序列并写出结果文件。
//
// 注：本文件中对各步骤尽量保持小而明确的职责划分，使得错误更易定位并便于在 CI 中逐步测试。
// checkOption：参数再校验与工作目录准备
// 该函数做了多项关键工作：
// - 使用 `file_io` 模块统一校验文件存在性与路径类型（以中央化错误信息与行为）；
// - 检查数值参数的合理性（threads/kmer_size/cons_n）；
// - 准备工作目录：在 Debug 模式允许非空工作目录以便快速迭代，在 Release 下要求 workdir 为空以避免覆盖已有数据；
// - 对 MSA 命令模板做一次自检（调用 cmd::testCommandTemplate），该函数会在临时目录运行模板并验证是否产生输出，
//   以便提前发现模板占位符或命令行问题。注意：testCommandTemplate 会实际运行外部命令，可能需要联网或外部工具存在；
//   如果不希望在参数校验阶段执行外部命令，可把此测试移到运行时并在出错时给用户明确说明。
static void checkOption(Options& opt) {
    // 文件相关：统一调用 file_io
    file_io::requireRegularFile(opt.input, "input");

    if (!opt.center_path.empty()) {
        file_io::requireRegularFile(opt.center_path, "center_path");
    }

    // 测试

    // 数值参数相关：仍在这里检查
    if (opt.threads <= 0) throw std::runtime_error("threads must be > 0");
    if (opt.kmer_size <= 0) throw std::runtime_error("kmer_size must be > 0");
    if (opt.kmer_window <= 0) throw std::runtime_error("kmer_window must be > 0");
    if (opt.cons_n <= 0) throw std::runtime_error("cons_n must be > 0");

    // k 与 w 的常见约束（贴近 minimap2 的可实现范围）
    if (opt.kmer_size > 31) throw std::runtime_error("kmer_size too large (must be <= 31)");
    if (opt.kmer_window >= 256) {
        spdlog::warn("kmer_window >= 256 may be slow/unsupported by the fast minimizer path; current value: {}", opt.kmer_window);
    }

    // workdir：Debug 下允许非空，Release 下必须空（与你现有逻辑一致）
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

int main(int argc, char** argv) {
    // main 的整体流程（高层概览）：
    // 1) 初始化日志与线程池
    // 2) 定义 CLI 并解析参数
    // 3) 打印并校验参数
    // 4) 准备 logger 输出到工作目录
    // 5) 运行预处理（fetch/copy + 清洗 + Top-K 选择）
    // 6) 对选出的共识未对齐序列调用外部 MSA（alignConsensusSequence）
    // 7) 若输入序列总数 <= cons_n，则直接复制对齐结果为最终输出并退出
    // 8) 否则调用 consensus::generateConsensusSequence（多线程/分批/可能双缓冲实现）产出最终共识
    // 9) 后续步骤：分组、双序列比对、合并等（TODO）

    try
    {
        spdlog::init_thread_pool(8192, 1);
        setupLogger();
        spdlog::info("Starting halign4...");

        Options opt;
        CLI::App app{"halign4"};

        setupCli(app, opt);

        // CLI11 的解析必须在 main 中进行，以便直接处理 argc/argv
        CLI11_PARSE(app, argc, argv);


        // 解析成功后，opt 已被填充
        logParsedOptions(opt);

        // 校验参数（文件存在性、工作目录、msa 模板自检等）
        checkOption(opt);

        // 将日志也写到 workdir 中（便于收集运行时日志）；此调用假定 workdir 已存在并已经被准备好
        setupLoggerWithFile(opt.workdir);

        // ---------------- 预处理阶段（IO 密集） ----------------
        // preprocessInputFasta 的职责：
        // - 如果 input 是 URL 则 download 到 workdir/data/raw，否则 copy 到该目录；
        // - 逐条读取、清洗并写出到 data/clean，同时选出最长的 cons_n 条序列写入共识输入文件；
        // - 返回处理的总序列数（用于决定下一步是否需要进一步的合并）
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
            // 如果用户直接指定了中心序列路径，则将其拷贝到工作目录中
            // 说明：
            // 1. 保持用户原始文件不变（不删除 opt.center_path）
            // 2. 删除工作目录中可能存在的旧文件（consensus_unaligned_file）
            // 3. 将用户指定的文件拷贝到工作目录，以便后续流程统一处理
            spdlog::info("Using user-specified center sequence: {}", opt.center_path);

            // 删除工作目录中的旧文件（如果存在）
            if (std::filesystem::exists(consensus_unaligned_file)) {
                file_io::removeAll(consensus_unaligned_file);
            }

            // 从用户指定的路径拷贝到工作目录
            file_io::copyFile(FilePath(opt.center_path), consensus_unaligned_file);
            spdlog::info("Center sequence copied to: {}", consensus_unaligned_file.string());
        }



        // 如果输入序列总数小于等于用于生成共识的数量（cons_n），说明我们已经在预处理阶段就处理完毕，
        // 此时直接把 align 输出拷贝到用户指定的最终输出文件并退出；这是一个常见的快速路径，避免无谓的合并工作。
        if (opt.center_path.empty() && preproc_count <= opt.cons_n)
        {
            // 调用外部 MSA（可能是 mafft/clustalo/其他高质量工具），会阻塞直到命令完成
            alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, opt.msa_cmd,  opt.threads);
            file_io::copyFile(consensus_aligned_file,FilePath(opt.output));
            spdlog::info("All sequences processed; final output written to {}", opt.output);
            spdlog::info("halign4 End!");
            return 0;
        }



        // ---------------- 比对阶段（使用 RefAligner） ----------------
        // 说明：
        // 1. 如果用户指定了 center_path，则从文件读取（consensus_string 为空）
        // 2. 否则使用内存中的 consensus_string（避免重复读取文件）
        // 3. consensus_string 会被保存到 RefAligner 的成员变量中
        FilePath ref_path = opt.center_path.empty() ? consensus_file : FilePath(opt.center_path);

        align::RefAligner ref_aligner(opt, ref_path);
        ref_aligner.alignQueryToRef(opt.input);
        ref_aligner.mergeAlignedResults(consensus_aligned_file, opt.msa_cmd);
        // 调用RefAligner进行后续的对齐和合并工作

        FilePath final_output_path = FilePath(opt.workdir) / RESULTS_DIR / FINAL_ALIGNED_FASTA;
        file_io::copyFile(final_output_path, FilePath(opt.output));
        spdlog::info("Final aligned output written to {}", opt.output);
        // 后续流程（概要）：
        // - 使用 minimizer / k-mer 方法估算序列间相似度并分组；
        // - 对于每个分组使用多序列比对并生成局部共识；

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
