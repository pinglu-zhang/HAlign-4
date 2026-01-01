//
// Created by 30451 on 2026/1/1.
//

#include <iostream>
#include <CLI/CLI.hpp>
#include <string>
#include <cstdint>

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
        ->check(CLI::Range(1, 1024));

    app.add_option("--kmer-size", opt.kmer_size, "K-mer size")
        ->default_val(15)
        ->check(CLI::Range(1, 4096));

    app.add_flag("--keep-length", opt.keep_length, "Keep sequence length unchanged");
}

int main(int argc, char** argv) {
    Options opt;
    CLI::App app{"HAlign4"};

    setup_cli(app, opt);

    // 关键：CLI11_PARSE 必须在返回 int 的函数里用（典型就是 main）
    CLI11_PARSE(app, argc, argv);

    // 解析成功后，opt 已被填充
    std::cout << "Parsed options:\n"
              << "  -i input       : " << opt.input << "\n"
              << "  -o output      : " << opt.output << "\n"
              << "  -w workdir     : " << opt.workdir << "\n"
              << "  -c center_path : " << (opt.center_path.empty() ? "(empty)" : opt.center_path) << "\n"
              << "  -p hq_cmd_path : " << (opt.msa_cmd_path.empty() ? "(empty)" : opt.msa_cmd_path) << "\n"
              << "  -t threads     : " << opt.threads << "\n"
              << "  --kmer-size    : " << opt.kmer_size << "\n"
              << "  --keep-length  : " << (opt.keep_length ? "true" : "false") << "\n";

    // TODO: 这里开始调用你的算法 pipeline
    // run_msa(opt);

    return 0;
}
