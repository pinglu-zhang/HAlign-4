#include <config.hpp>
#include <utils.h>
#include <thread>
#include "preprocess.h"
#include "consensus.h"
#include <chrono>

struct Options {
    std::string input;          // -i
    std::string output;         // -o
    std::string workdir;        // -w

    std::string center_path;    // -c
    std::string msa_cmd;    // -p

    int threads = [](){ unsigned int hc = std::thread::hardware_concurrency(); return static_cast<int>(hc ? hc : 1u); }();            // -t
    int kmer_size = 15;         // --kmer-size
    int cons_n = 1000;         // --cons-n

    bool keep_length = false;   // --keep-length
};

// 这里不再返回 Options，也不调用 CLI11_PARSE
// 只负责：把参数定义到 app，并绑定到 opt
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

// 新增：将解析后参数的日志打印抽象为一个函数，便于后续增加参数时扩展
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
    // 检验所有的参数是否合法，如果不合法打印并且退出
    try
    {
        spdlog::init_thread_pool(8192, 1);
        setupLogger();
        spdlog::info("Starting halign4...");

        Options opt;
        CLI::App app{"halign4"};

        setupCli(app, opt);

        // 关键：CLI11_PARSE 必须在返回 int 的函数里用（典型就是 main）
        CLI11_PARSE(app, argc, argv);


        // 解析成功后，opt 已被填充
        logParsedOptions(opt);

        // 校验参数
        checkOption(opt);

        setupLoggerWithFile(opt.workdir);

        // TODO: 这里开始调用你的算法 pipeline
        // 预处理原始数据
        // 输入文件路径，工作目录
        // 输出清理好的数据和共识序列fasta文件
        uint_t preproc_count = preprocessInputFasta(opt.input, opt.workdir, opt.cons_n);
        spdlog::info("Preprocessing produced {} records", preproc_count);

        // 获取共识序列
        // 输入共识序列fasta文件路径，调用别的方法比对
        // 最后返回共识序列
        FilePath consensus_unaligned_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_UNALIGNED;
        FilePath consensus_aligned_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_ALIGNED;
        FilePath consensus_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_FASTA;
        FilePath consensus_json_file = FilePath(opt.workdir) / WORKDIR_DATA / DATA_CLEAN/ CLEAN_CONS_JSON;
        if (!opt.center_path.empty())
        {
            consensus_unaligned_file = opt.center_path;
        }
        alignConsensusSequence(consensus_unaligned_file, consensus_aligned_file, opt.msa_cmd, opt.workdir, opt.threads);
        if (preproc_count <= opt.cons_n)
        {
            file_io::copyFile(consensus_aligned_file,FilePath(opt.output));
            spdlog::info("All sequences processed; final output written to {}", opt.output);
            spdlog::info("halign4 End!");
            return 0;
        }

        // 计时：测量 generateConsensusSequenceDoubleBuffered 的执行时间
        const std::size_t batch_size = 512; // 可调：每轮处理的序列数
        spdlog::info("Starting consensus generation (double-buffered), batch_size={}...", batch_size);
        auto t_start = std::chrono::steady_clock::now();
        std::string consensus_string = consensus::generateConsensusSequenceDoubleBuffered(
            consensus_aligned_file,
            consensus_file, // 覆盖写回
            consensus_json_file,
            0, // 不限制数量
            opt.threads,
            batch_size
        );

        auto t_end = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
        spdlog::info("Consensus generation finished, elapsed: {:.3f} s", elapsed);


        // 提取minimzer估算相似度分组
        // 输入清理好的数据fasta文件路径，工作目录，kmer_size
        // 输出相似度和fasta文件路径

        // 分组多序列比对
        // 输入相似度和fasta文件路径
        // 写出所有双序列比对结果和共识序列到文件

        // 合并得到最终结果

        // run_msa(opt);

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
