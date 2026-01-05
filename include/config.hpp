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

const std::string DATA_RAW = "raw_data";        // 原始数据子目录（下载/拷贝来的原始输入）
const std::string DATA_CLEAN = "clean_data";    // 清理后数据子目录（预处理后写出的序列）

// 共识相关的文件名（相对于 DATA_CLEAN）
const std::string CLEAN_CONS_UNALIGNED = "consensus_unaligned.fasta"; // 共识序列文件名（未对齐）
const std::string CLEAN_CONS_ALIGNED = "consensus_aligned.fasta";     // 共识序列文件名（已对齐）

const std::string CLEAN_CONS_FASTA = "consensus.fasta"; // 最终共识序列 FASTA 文件名
const std::string CLEAN_CONS_JSON = "consensus.json";   // 共识统计/计数输出（JSON）

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

