#ifndef CONFIG_HPP
#define CONFIG_HPP

// ------------------------------------------------------------------
// config.hpp
// 说明（详细中文注释）
//
// 本文件集中定义了全局配置常量、日志初始化函数、以及与命令行解析（CLI11）相关的
// 辅助类型与工具（例如自定义格式器与 validator）。
//
// 目的：
// - 为项目提供单一入口的“配置库”，便于在代码各处引用一致的符号（例如工作目录结构、默认命令模板、日志文件名等）；
// - 提供便捷的日志初始化函数 `setupLogger` / `setupLoggerWithFile`，方便在 main 中统一配置日志输出到控制台与文件；
// - 提供 CLI 美化与输入修剪（trim_whitespace），提高命令行体验与容错性。
//
// 注意事项：
// - 这里的默认 MSA 命令模板（DEFALT_MSA_CMD）仅做占位和演示用途；生产环境请替换为实际可用的 MSA 工具与参数。
// - 所有路径常量为字符串字面量（相对路径），在使用时通常与 `workdir` 进行拼接以形成绝对或工作相对路径。
// ------------------------------------------------------------------

// ------------------------------------------------------------------
// 引入头文件：功能模块包括线程池、命令行解析、日志系统、序列化库等
// ------------------------------------------------------------------
#include <CLI/CLI.hpp>                       // CLI11 命令行解析库（本地版本）
#include "spdlog/spdlog.h"                       // spdlog 主头文件
#include "spdlog/sinks/stdout_color_sinks.h"     // 控制台彩色输出 sink
#include "spdlog/sinks/basic_file_sink.h"        // 文件输出 sink
#include "spdlog/async.h"                        // 异步日志支持

#include <filesystem>
#include <sstream>
#include <cinttypes>

// ------------------------------------------------------------------
// 通用配置常量
// ------------------------------------------------------------------
#define VERSION "1.2.0"                   // 版本号，程序启动时可打印以便追踪
#define LOGGER_NAME "logger"              // 默认日志器名称（用于 spdlog 注册）
#define LOGGER_FILE "halign4.log"         // 默认日志文件名（相对于工作目录）
#define CONFIG_FILE "config.json"         // 默认配置文件路径（如果将来支持外部配置）


const std::string DEFALT_MSA_CMD = "minipoa {input} -S -t {thread} -r1 > {output}"; // 默认多序列比对命令模板

// 工作目录体系
const std::string WORKDIR_DATA = "data";         // 原始数据目录
const std::string WORKDIR_TMP = "temp";         // 临时目录（短生命周期文件）
const std::string RESULTS_DIR = "result";      // 最终结果目录

const std::string DATA_RAW = "raw_data";        // 原始数据子目录（下载/拷贝来的原始输入）
const std::string DATA_CLEAN = "clean_data";    // 清理后数据子目录（预处理后写出的序列）

// 共识相关的文件名（相对于 DATA_CLEAN）
const std::string CLEAN_CONS_UNALIGNED = "consensus_unaligned.fasta"; // 共识序列文件名（未对齐）
const std::string CLEAN_CONS_ALIGNED = "consensus_aligned.fasta";     // 共识序列文件名（已对齐）

const std::string CLEAN_CONS_FASTA = "consensus.fasta"; // 最终共识序列 FASTA 文件名
const std::string CLEAN_CONS_JSON = "consensus.json";   // 共识统计/计数输出（JSON）

// ------------------------------------------------------------------
// 对齐输出相关的文件名（相对于 RESULTS_DIR）
// ------------------------------------------------------------------
// 说明：这些文件名用于多序列比对（MSA）和插入序列处理的中间文件与最终输出文件
// 目的：集中管理文件名常量，避免在代码中硬编码字符串字面量，便于统一修改与维护

// 最终对齐结果文件名（所有序列对齐后的 MSA 输出）
#define FINAL_ALIGNED_FASTA "final_aligned.fasta"

// 插入序列相关文件名
#define ALL_INSERTION_FASTA "all_insertion.fasta"           // 合并的所有插入序列（未对齐）
#define ALIGNED_INSERTION_FASTA "aligned_insertion.fasta"   // 对齐后的插入序列（MSA 结果）

// 线程级别输出文件名模板（用于并行写出）
// 说明：每个线程独立写出 SAM 文件，避免线程竞争；最后由主线程合并
#define THREAD_SAM_PREFIX "thread"                          // 线程 SAM 文件前缀（"thread" + tid + ".sam"）
#define THREAD_SAM_SUFFIX ".sam"                            // 线程 SAM 文件后缀
#define THREAD_INSERTION_SAM_SUFFIX "_insertion.sam"        // 线程插入序列 SAM 文件后缀（"thread" + tid + "_insertion.sam"）

// ------------------------------------------------------------------
// 调试与整数精度配置
// ------------------------------------------------------------------
#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef M64
#define M64 0
#endif

// 根据宏定义切换 32 位或 64 位整数
// 目的：在处理非常大量的数据计数时，使用 64-bit 能避免溢出；开发/轻量运行可用 32-bit 节省内存。
#if M64
typedef int64_t	int_t;
typedef uint64_t uint_t;
#define PRIdN	PRId64
#define U_MAX	UINT64_MAX
#define I_MAX	INT64_MAX
#define I_MIN	INT64_MIN
#else
typedef int32_t int_t;
typedef uint32_t uint_t;
#define PRIdN	PRId32
#define U_MAX	UINT32_MAX
#define I_MAX	INT32_MAX
#define I_MIN	INT32_MIN
#endif

// 获取硬件并发线程数（兜底 1）
static int get_default_threads() {
    unsigned int hc = std::thread::hardware_concurrency();
    return static_cast<int>(hc ? hc : 1u);
}

struct Options {
	// 输入/输出与工作目录
	std::string input;          // -i：输入序列文件（路径或压缩文件）
	std::string output;         // -o：最终输出文件（写入位置）
	std::string workdir;        // -w：工作目录，所有中间文件（data/raw, data/clean 等）放在该目录下

	// 可选参数：中心序列、MSA 命令模板
	std::string center_path;    // -c：可选，指定中心序列文件路径，若指定则绕过自动选择
	std::string msa_cmd;        // -p：用于对共识序列做 MSA 的命令模板（可以包含 {input} {output} {thread} 占位符）

	// 并行与算法参数
	int threads = get_default_threads(); // -t：线程数，默认为 CPU 核心数
	int kmer_size = 15;         // --kmer-size：用于归类/聚类的 k-mer 大小（后续步骤使用）
	int kmer_window = 10;       // --kmer-window：minimizer 窗口大小 w（以 k-mer 为单位）
	int cons_n = 1000;          // --cons-n：挑选用于共识计算的序列数量（Top-K by length）
	int sketch_size = 2000;     // --sketch-size：用于 sketch 的大小（默认 2000）

	// keep length 相关开关：
	// - keep_first_length：仅保持“第一条/中心序列”的长度不变（其余序列允许按对齐结果变化/填充），适用于只关心输出共识/中心序列长度的场景。
	// - keep_all_length  ：保持“所有中心序列”的长度不变，适用于后续流程严格依赖原始长度坐标系（代价通常更高）。
	bool keep_first_length = false; // --keep-first-length
	bool keep_all_length = false;   // --keep-all-length

	// workdir 清理开关：
	// - save_workdir：若为 true，则在完成比对/输出后保留工作目录；若为 false（默认），则在成功完成后删除工作目录。
	bool save_workdir = false;      // --save-workdir
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
        ->default_val(get_default_threads())
        ->check(CLI::Range(1, 100000));

    app.add_option("--kmer-size", opt.kmer_size, "K-mer size")
        ->default_val(15)
        ->check(CLI::Range(4, 31));

    app.add_option("--kmer-window", opt.kmer_window, "Minimizer window size (w, in number of k-mers)")
        ->default_val(10)
        ->check(CLI::Range(1, 1000000));

    app.add_option("--cons-n", opt.cons_n, "Number of sequences for consensus")
        ->default_val(1000)
        ->check(CLI::Range(1, 1000000));

    // 新增 sketch-size 参数
    app.add_option("--sketch-size", opt.sketch_size, "Sketch size")
        ->default_val(2000)
        ->check(CLI::Range(1, 10000000));

    // keep length：拆分为两个互斥开关（更精确地表达需求）
    app.add_flag("--keep-first-length", opt.keep_first_length,
        "Keep only the first/center sequence length unchanged");
    app.add_flag("--keep-all-length", opt.keep_all_length,
        "Keep all center sequences length unchanged");

    // workdir 管理：是否在完成后保留工作目录
    app.add_flag("--save-workdir", opt.save_workdir,
        "Keep the working directory after alignment completion (default: remove)");
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
        {"kmer_window", std::to_string(opt.kmer_window)},
        {"cons_n", std::to_string(opt.cons_n)},
        {"sketch_size", std::to_string(opt.sketch_size)},
        {"keep_first_length", boolToStr(opt.keep_first_length)},
        {"keep_all_length", boolToStr(opt.keep_all_length)},
        {"save_workdir", boolToStr(opt.save_workdir)}
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

// ------------------------------------------------------------------
// CLI11 自定义格式器（美化选项输出）
// 说明：
// - 自定义 `make_option_opts` 可以在帮助中显示参数类型与默认值，便于用户理解；
// - 自定义 `make_usage` 提供更友好的使用示例与说明，方便新手快速上手。
// ------------------------------------------------------------------
class CustomFormatter : public CLI::Formatter {
public:
	CustomFormatter() : Formatter() {}

	// 自定义参数展示样式（带默认值）
	std::string make_option_opts(const CLI::Option* opt) const override {
		if (opt->get_type_size() == 0) return "";
		std::ostringstream out;
		out << " " << opt->get_type_name();
		if (!opt->get_default_str().empty())
			out << " (default: " << opt->get_default_str() << ")";
		return out.str();
	}

	// 提供使用示例；此处的 Example 可根据项目实际可执行名称更新
	std::string make_usage(const CLI::App* app, std::string name) const override {
		std::ostringstream out;
		out << "Usage:\n"
			<< "  ./halign4 -i <ref.fa> -o <output.fa> -w </path/to/workdir> [options]\n\n"
			<< "Example:\n"
			<< "  ./halign4 -i ref.fa -o output.fa -w ./tmp -t 8\n\n";
		return out.str();
	}
};

// ------------------------------------------------------------------
// CLI11 自定义 validator：自动去除参数两侧空白
// 说明：有时用户在命令行中误加空格或复制粘贴带有换行，trim_whitespace 可以提高健壮性
// ------------------------------------------------------------------
inline CLI::Validator trim_whitespace = CLI::Validator(
	[](std::string& s) {
		auto start = s.find_first_not_of(" \t\n\r");
		auto end = s.find_last_not_of(" \t\n\r");
		s = (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
		return std::string();  // 空字符串表示验证通过
	}, ""
);



// ------------------------------------------------------------------
// 日志系统初始化
// 说明（中文注释）：
// - `setupLoggerWithFile(path)` 会创建一个异步的 spdlog 日志器，输出到控制台与指定目录下的日志文件；
// - 异步日志（async_logger）在高并发日志场景下能降低阻塞与 I/O 延迟，但需要在程序启动时配置线程池；
// - `setupLogger()` 仅输出到控制台，适合交互式调试或短期运行；
// - 两个函数都会设置默认日志级别（Debug/Info）并定期刷盘（flush_every），这有助于在崩溃时保留日志。
// ------------------------------------------------------------------
inline void setupLoggerWithFile(std::filesystem::path log_dir) {
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] [%l] %v%$");

	std::filesystem::path log_file = log_dir / LOGGER_FILE;
	auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file.string(), true);
	file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");

	spdlog::sinks_init_list sinks = { console_sink, file_sink };
	auto logger = std::make_shared<spdlog::async_logger>(
		LOGGER_NAME, sinks.begin(), sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);

	spdlog::set_default_logger(logger);
#ifdef _DEBUG
	spdlog::set_level(spdlog::level::debug);
#else
	spdlog::set_level(spdlog::level::info);
#endif
	spdlog::flush_every(std::chrono::seconds(3));
}

// 控制台日志（用于开发/调试）
inline void setupLogger() {
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] [%l] %v%$");

	spdlog::sinks_init_list sinks = { console_sink };
	auto logger = std::make_shared<spdlog::async_logger>(
		LOGGER_NAME, sinks.begin(), sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);

	spdlog::set_default_logger(logger);
#ifdef _DEBUG
	spdlog::set_level(spdlog::level::debug);
#else
	spdlog::set_level(spdlog::level::info);
#endif

	spdlog::flush_every(std::chrono::seconds(3));
}

// 获取完整命令行字符串（用于日志记录和重现）
inline std::string getCommandLine(int argc, char** argv) {
	std::ostringstream cmd;
	for (int i = 0; i < argc; ++i) {
		cmd << argv[i];
		if (i != argc - 1) cmd << " ";
	}
	return cmd.str();
}

#endif // CONFIG_HPP

