//
// Created by 30451 on 2026/1/1.
//

#include <iostream>
#include <config.hpp>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <stdexcept>

struct Options {
    std::string input;          // -i
    std::string output;         // -o
    std::string workdir;        // -w

    std::string center_path;    // -c
    std::string msa_cmd_path;    // -p

    int threads = 1;            // -t
    int kmer_size = 15;         // --kmer-size

    bool keep_length = false;   // --keep-length
};

// 这里不再返回 Options，也不调用 CLI11_PARSE
// 只负责：把参数定义到 app，并绑定到 opt
static void setup_cli(CLI::App& app, Options& opt) {
    app.description("HAlign4 / MSA tool");

    // 必须参数
    app.add_option("-i", opt.input, "Input sequences (file path)")
        ->required()
        ->check(CLI::ExistingFile);

    app.add_option("-o", opt.output, "Output result (file path)")
        ->required();

    app.add_option("-w", opt.workdir, "Working directory")
        ->required();

    // 可选参数
    app.add_option("-c", opt.center_path, "Center sequence path")
        ->check(CLI::ExistingFile);

    // 如果 -p 是“可执行文件路径”，ExistingFile 通常也能用；
    // 若你希望允许仅命令名（在 PATH 中），这里就不要 check
    app.add_option("-p", opt.msa_cmd_path, "High-quality method command path")
        ->check(CLI::ExistingFile);

    app.add_option("-t", opt.threads, "Number of threads")
        ->default_val(1)
        ->check(CLI::Range(1, 100000));

    app.add_option("--kmer-size", opt.kmer_size, "K-mer size")
        ->default_val(15)
        ->check(CLI::Range(1, 4096));

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
        {"msa_cmd_path", toString(opt.msa_cmd_path, valW)},
        {"threads", std::to_string(opt.threads)},
        {"kmer_size", std::to_string(opt.kmer_size)},
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

// 新增：检查选项的合法性并处理 workdir（必须为空文件夹；如果不存在则创建）
static void checkOption(const Options& opt) {
    namespace fs = std::filesystem;

    // input must exist
    if (opt.input.empty()) {
        throw std::runtime_error("input is empty");
    }
    if (!fs::exists(opt.input)) {
        throw std::runtime_error("input file does not exist: " + opt.input);
    }

    // optional center_path must exist if provided
    if (!opt.center_path.empty() && !fs::exists(opt.center_path)) {
        throw std::runtime_error("center_path does not exist: " + opt.center_path);
    }

    // optional msa_cmd_path must exist if provided
    if (!opt.msa_cmd_path.empty() && !fs::exists(opt.msa_cmd_path)) {
        throw std::runtime_error("msa_cmd_path does not exist: " + opt.msa_cmd_path);
    }

    if (opt.threads <= 0) {
        throw std::runtime_error("threads must be > 0");
    }
    if (opt.kmer_size <= 0) {
        throw std::runtime_error("kmer_size must be > 0");
    }

    if (opt.workdir.empty()) {
        throw std::runtime_error("workdir is empty");
    }

    fs::path wd{opt.workdir};

    std::error_code ec;
    if (!fs::exists(wd, ec)) {
        // try to create
        if (!fs::create_directories(wd, ec) || ec) {
            throw std::runtime_error("failed to create workdir: " + opt.workdir + " (" + ec.message() + ")");
        }
        spdlog::info("Created workdir: {}", opt.workdir);
    } else {
        if (!fs::is_directory(wd, ec) || ec) {
            throw std::runtime_error("workdir exists but is not a directory: " + opt.workdir);
        }
        // ensure empty
        auto it = fs::directory_iterator(wd, ec);
        if (ec) {
            throw std::runtime_error("cannot read workdir: " + opt.workdir + " (" + ec.message() + ")");
        }
#ifndef _DEBUG
        if (it != fs::end(it)) {
            throw std::runtime_error("workdir must be empty: " + opt.workdir);
        }
#endif
    }
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

        setup_cli(app, opt);

        // 关键：CLI11_PARSE 必须在返回 int 的函数里用（典型就是 main）
        CLI11_PARSE(app, argc, argv);
        // 解析成功后，opt 已被填充
        logParsedOptions(opt);

        // 校验参数
        checkOption(opt);

        // TODO: 这里开始调用你的算法 pipeline
        // run_msa(opt);

        return 0;
    } catch (const std::exception &e) {
        spdlog::error("Fatal error: {}", e.what());
        spdlog::error("halign4 End!");
        return 1;
    } catch (...) {
        spdlog::error("Fatal error: unknown exception");
        spdlog::error("halign4 End!");
        return 2;
    }

    return 0;
}
