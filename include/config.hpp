#ifndef CONFIG_HPP
#define CONFIG_HPP

// ------------------------------------------------------------------
// 引入头文件：功能模块包括线程池、命令行解析、日志系统、序列化库等
// ------------------------------------------------------------------
#include <CLI/CLI.hpp>                       // CLI11 命令行解析库（本地版本）
#include "spdlog/spdlog.h"                       // spdlog 主头文件
#include "spdlog/sinks/stdout_color_sinks.h"     // 控制台彩色输出 sink
#include "spdlog/sinks/basic_file_sink.h"        // 文件输出 sink
#include "spdlog/async.h"                        // 异步日志支持

// ------------------------------------------------------------------
// 通用配置常量
// ------------------------------------------------------------------
#define VERSION "1.0.0"                   // 版本号
#define LOGGER_NAME "logger"              // 默认日志器名称
#define LOGGER_FILE "halign4.log"          // 默认日志文件名
#define CONFIG_FILE "config.json"         // 配置文件路径


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

	// 提供使用示例
	std::string make_usage(const CLI::App* app, std::string name) const override {
		std::ostringstream out;
		out << "Usage:\n"
			<< "  ./halign4 -i <ref.fa> -o <output.fa> -w </path/to/workdir> [options]\n\n"
			<< "Example:\n"
			<< "  ./RaMAx -r ref.fa -q query.fa -o output/ -w ./tmp -t 8\n\n";
		return out.str();
	}
};

// ------------------------------------------------------------------
// CLI11 自定义 validator：自动去除参数两侧空白
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
// ------------------------------------------------------------------

// 控制台 + 文件日志
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
	spdlog::set_level(spdlog::level::trace);
	spdlog::flush_every(std::chrono::seconds(3));
}

// 控制台日志
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