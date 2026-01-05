#include <config.hpp>
#include <utils.h>
#include <thread>
#include "preprocess.h"
#include "consensus.h"
#include <chrono>

// 该文件是程序的入口（main），负责：
// 1) 解析命令行参数（使用 CLI11），并把值绑定到 Options 结构；
// 2) 校验参数（例如文件存在性、工作目录准备、外部 MSA 命令模板测试）；
// 3) 调用预处理（preprocessInputFasta）把原始输入复制/下载到工作目录并清洗；
// 4) 执行针对共识序列的多序列比对（外部工具，通过命令模板），并产生对齐结果；
// 5) 以多线程并行/分批的方式生成最终共识序列并写出结果文件。
//
// 注：本文件中对各步骤尽量保持小而明确的职责划分，使得错误更易定位并便于在 CI 中逐步测试。

struct Options {
    // 输入/输出与工作目录
    std::string input;          // -i：输入序列文件（路径或压缩文件）
    std::string output;         // -o：最终输出文件（写入位置）
    std::string workdir;        // -w：工作目录，所有中间文件（data/raw, data/clean 等）放在该目录下

    // 可选参数：中心序列、MSA 命令模板
    std::string center_path;    // -c：可选，指定中心序列文件路径，若指定则绕过自动选择
    std::string msa_cmd;        // -p：用于对共识序列做 MSA 的命令模板（可以包含 {input} {output} {thread} 占位符）

    // 并行与算法参数
    int threads = [](){ unsigned int hc = std::thread::hardware_concurrency(); return static_cast<int>(hc ? hc : 1u); }(); // -t：线程数，默认为 CPU 核心数
    int kmer_size = 15;         // --kmer-size：用于归类/聚类的 k-mer 大小（后续步骤使用）
    int cons_n = 1000;          // --cons-n：挑选用于共识计算的序列数量（Top-K by length）

    bool keep_length = false;   // --keep-length：是否在后续处理保持原始序列长度（影响共识/校正策略）
};

// setupCli：定义 CLI 参数并绑定到 Options
// 注释说明：
// - 使用 CLI11 库实现参数解析，支持短参数/长参数与基本校验（例如 ExistingFile）
// - 对于像 -p/--msa-cmd 这类可能是可执行名（而非完整路径）的参数，ExistingFile 会拒绝仅命令名的情况；
//   如果希望允许命令名（在 PATH 中解析），可以去掉 check(CLI::ExistingFile) 或改为用户层面的更宽容判断。
static void setupCli(CLI::App& app, Options& opt) {
    app.description("HAlign4 / MSA tool");

   // 必须参数（同时支持短参数和长参数）
    app.add_option("-i,--input", opt.input, "Input sequences (file path)")
        ->required()
        ->check(CLI::ExistingFile);

    app.add_option("-o,--output", opt.output, "Output result (file path)")
        ->required();

    app.add_option("-w,--workdir", opt.workdir, "Working directory")
        ->required();

    // 可选参数（增加长参数形式）
    app.add_option("-c,--center-path", opt.center_path, "Center sequence path")
        ->check(CLI::ExistingFile);

    // 如果 -p 是“可执行文件路径”，ExistingFile 通常也能用；
    // 若你希望允许仅命令名（在 PATH 中），这里就不要 check
    app.add_option("-p,--msa-cmd", opt.msa_cmd, "High-quality method command path")
        ->check(CLI::ExistingFile);

    app.add_option("-t,--thread", opt.threads, "Number of threads")
        ->default_val([](){ unsigned int hc = std::thread::hardware_concurrency(); return static_cast<int>(hc ? hc : 1u); }())
        ->check(CLI::Range(1, 100000));

    app.add_option("--kmer-size", opt.kmer_size, "K-mer size")
        ->default_val(15)
        ->check(CLI::Range(1, 4096));

    app.add_option("--cons-n", opt.cons_n, "Number of sequences for consensus")
        ->default_val(1000)
        ->check(CLI::Range(1, 1000000));

    app.add_flag("--keep-length", opt.keep_length, "Keep sequence length unchanged");
}

// logParsedOptions：把解析后的参数以漂亮的表格形式输出到日志
// 说明：此处的输出用于帮助用户和调试（打印被截断的长字符串、boolean 友好显示等），
// 不影响程序行为。若程序在无头环境运行（服务/容器），日志也便于审计和复现运行参数。
static void logParsedOptions(const Options& opt) {
    // Helper to convert values and truncate long strings for tidy display
    auto toString = [](const std::string& s, size_t maxLen) -> std::string {
        if (s.empty()) return "(empty)";
        if (s.size() <= maxLen) return s;
        return s.substr(0, maxLen - 3) + "...";
    };

    auto boolToStr = [](bool b) { return b ? "true" : "false"; };

    const size_t keyW = 14;
    const size_t valW = 60;
    const size_t innerW = keyW + 3 + valW; // "key : value"

    std::vector<std::pair<std::string, std::string>> rows = {
        {"input", toString(opt.input, valW)},
        {"output", toString(opt.output, valW)},
        {"workdir", toString(opt.workdir, valW)},
        {"center_path", toString(opt.center_path, valW)},
        {"msa_cmd", toString(opt.msa_cmd, valW)},
        {"threads", std::to_string(opt.threads)},
        {"kmer_size", std::to_string(opt.kmer_size)},
        {"cons_n", std::to_string(opt.cons_n)},
        {"keep_length", boolToStr(opt.keep_length)}
    };

    std::ostringstream oss;

    // top border
    oss << "+" << std::string(innerW, '-') << "+\n";

    // title centered
    const std::string title = " Parsed options ";
    size_t paddingLeft = 0;
    if (innerW > title.size()) paddingLeft = (innerW - title.size()) / 2;
    oss << "|" << std::string(paddingLeft, ' ') << title
        << std::string(innerW - paddingLeft - title.size(), ' ') << "|\n";

    // separator
    oss << "+" << std::string(innerW, '-') << "+\n";

    // rows
    for (auto &kv : rows) {
        oss << "| " << std::left << std::setw(keyW) << kv.first << " : "
            << std::setw(valW) << kv.second << "|\n";
    }

    // bottom border
    oss << "+" << std::string(innerW, '-') << "+";

    spdlog::info("\n{}", oss.str());
}

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
    if (opt.cons_n <= 0) throw std::runtime_error("cons_n must be > 0");

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
            // 如果用户直接指定了中心序列路径，则用该路径替换自动生成的共识输入，不再从预处理选择
            consensus_unaligned_file = opt.center_path;
        }

        // 调用外部 MSA（可能是 mafft/clustalo/其他高质量工具），会阻塞直到命令完成
        alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, opt.msa_cmd, opt.workdir, opt.threads);

        // 如果输入序列总数小于等于用于生成共识的数量（cons_n），说明我们已经在预处理阶段就处理完毕，
        // 此时直接把 align 输出拷贝到用户指定的最终输出文件并退出；这是一个常见的快速路径，避免无谓的合并工作。
        if (preproc_count <= opt.cons_n)
        {
            file_io::copyFile(consensus_aligned_file,FilePath(opt.output));
            spdlog::info("All sequences processed; final output written to {}", opt.output);
            spdlog::info("halign4 End!");
            return 0;
        }

        // ---------------- 共识生成阶段（计算密集/内存密集） ----------------
        // 说明：generateConsensusSequence 可能有多种实现（单线程/多线程/双缓冲/批处理/SoA），
        // 这里调用的是现有命名空间下的实现并传入线程数。实现应当考虑如下优化点：
        // - 按批读取对齐序列，避免一次性把整个文件加载到内存；
        // - 对每个位置进行并行计数（按列并行或线程本地计数后合并）以避免频繁原子操作与 false sharing；
        // - 使用 SoA（Structure-of-Arrays）便于向量化与缓存友好；
        // - 在 Release 模式启用 -O3 与 -march=native 等编译器优化，并根据目标架构调整线程亲和性（affinity）；
        // - 记录运行时间以便基准分析。
        const std::size_t batch_size = 512; // 可调：每轮处理的序列数（若实现使用批次）；对内存与并发有直接影响
        spdlog::info("Starting consensus generation (double-buffered), batch_size={}...", batch_size);
        auto t_start = std::chrono::steady_clock::now();
        std::string consensus_string = consensus::generateConsensusSequence(
            consensus_aligned_file,
            consensus_file, // 覆盖写回
            consensus_json_file,
            0, // 不限制数量
            opt.threads
        );

        auto t_end = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
        spdlog::info("Consensus generation finished, elapsed: {:.3f} s", elapsed);


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
