#include "preprocess.h"
#include <chrono>
// 新增头文件，用于调用外部命令、构造字符串和检查文件存在性
#include <cstdlib>
#include <sstream>
#include <filesystem>

// 本文件包含输入 FASTA 的预处理逻辑：
// - 将输入文件（可以是本地路径或 URL）获取到工作目录的 data/raw 下；
// - 对序列做简单清洗（大写化，非 AGCTU 字符替换为 N）；
// - 将清洗后的序列写入 data/clean 下；
// - 维护一个 Top-K 选择器，选择长度最长的前 K 条序列（稳定保留较早出现的序列）；
// - 最终把清洗后的数据与选中的 consensus 输出到指定文件（在此文件中只负责选择与写出接口，具体写出由调用者处理）。

// 解释：此处使用的 FilePath 为 file_io::FilePath（即 std::filesystem::path 的别名），
// seq_io 命名空间封装了读取/写入 FASTA 的细节（openKseqReader / FastaWriter / SeqRecord 等）。

uint_t preprocessInputFasta(const std::string input_path, const std::string workdir, const int cons_n) {
    // 参数说明：
    // - input_path: 输入 FASTA 的路径或 URL（字符串）。
    // - workdir: 工作目录路径（应当已经准备好或由调用方保证），在此目录下会创建 data/raw 和 data/clean 等子目录。
    // - cons_n: 要保留的最长序列数量（Top-K）。

    // 计时用于日志记录，帮助性能分析
    const auto t_start = std::chrono::steady_clock::now();

    spdlog::info("Preprocessing input FASTA file: {}", input_path);
    spdlog::info("Working directory: {}", workdir);

    // ---------- 目录准备 ----------
    // 1) 确保工作目录下存在 data 文件夹。
    //    该目录用于存放原始与清洗后的数据：data/raw 和 data/clean。
    FilePath data_dir = FilePath(workdir) / WORKDIR_DATA;
    file_io::ensureDirectoryExists(data_dir);
    spdlog::info("Ensured data directory exists: {}", data_dir.string());

    // 2) 在 data 下创建 raw_data 文件夹，用于保存原始（未清洗）输入。
    FilePath raw_data_dir = data_dir / DATA_RAW;
    file_io::ensureDirectoryExists(raw_data_dir);
    spdlog::info("Ensured raw data directory exists: {}", raw_data_dir.string());

    // 3) 在 data 下创建 clean_data 文件夹，用于保存清洗后的输出。
    FilePath clean_data_dir = data_dir / DATA_CLEAN;
    file_io::ensureDirectoryExists(clean_data_dir);
    spdlog::info("Ensured clean data directory exists: {}", clean_data_dir.string());

    // ---------- 获取输入文件（支持本地/远程） ----------
    // 4) 将输入文件复制或下载到 raw_data 下。
    //    这里调用 file_io::fetchFile，函数内部会判断是 URL 还是本地路径并做相应操作（download 或 copy）。
    FilePath input_file = FilePath(input_path);
    FilePath raw_dest_file = raw_data_dir / input_file.filename();

    spdlog::info("Fetching input to working raw path: {} -> {}", input_file.string(), raw_dest_file.string());
    // 说明：fetchFile 在遇到远程 URL 时会调用 downloadFile 将数据写入本地；在本地路径时会调用 copyFile（包含跨设备回退等逻辑）。
    // 这里可能抛异常（例如网络错误、权限问题），调用方应在上层捕获并处理。预处理阶段假定 fetchFile 成功。
    file_io::fetchFile(input_file,raw_dest_file);
    spdlog::info("Input available at: {}", raw_dest_file.string());

    // ---------- 准备输出文件名 ----------
    // 5) 打开 raw 文件并逐条读取；对每条序列进行清洗（cleanSequence），写入 clean_data
    //    同时维护 TopKLongestSelector，选择最长的 cons_n 条序列（用于之后的 consensus 生成）。
    // handle input filenames like `sample.fasta.gz` -> `sample.fasta`
    FilePath in_fname = input_file.filename();
    std::string in_name = in_fname.string();
    const std::string comp_ext = ".gz";
    if (in_name.size() > comp_ext.size() &&
        in_name.compare(in_name.size() - comp_ext.size(), comp_ext.size(), comp_ext) == 0) {
        in_name.resize(in_name.size() - comp_ext.size());
        spdlog::info("Detected compressed input; using output name: {}", in_name);
    }
    FilePath clean_dest_file = clean_data_dir / FilePath(in_name);
    FilePath consensus_file = clean_data_dir / CLEAN_CONS_UNALIGNED;
    spdlog::info("Clean output: {} ; Consensus output: {}", clean_dest_file.string(), consensus_file.string());

    // ---------- 打开 reader/writer 与 TopK 选择器 ----------
    // seq_io::openKseqReader 返回一个 reader 指针（抽象），用于逐条读取序列；
    // seq_io::FastaWriter 用于把清洗后的序列写入到目标文件。
    auto reader = seq_io::openKseqReader(raw_dest_file);
    seq_io::FastaWriter clean_writer(clean_dest_file);
    TopKLongestSelector selector(cons_n);

    // 处理循环：读取 -> 清洗 -> 写出 -> 交给 TopK 选择器
    seq_io::SeqRecord rec;
    std::size_t total_records = 0;
    const std::size_t log_interval = 10000;

    // 重要说明（性能与正确性）：
    // - 该循环为预处理的热路径，若输入很大（成千上万 / 几百万条序列），要关注 IO 与内存占用。
    // - 性能优化点：
    //    * 使用 seq_io 的 KseqReader（基于 fread/gzread）并配合大缓冲可以显著提升读取吞吐；
    //    * FastaWriter::write 会把一条记录的 header 与折行后的序列缓存在临时字符串中一次性写出，
    //      避免逐字符写入带来的系统调用开销；这对写大文件非常重要；
    //    * TopKLongestSelector 应该实现为维护一个大小为 K 的最小堆，插入/替换成本为 O(log K)，适合 K 远小于记录总数的场景；
    // - 内存权衡：TopK 的实现会保留 K 条完整记录（占用内存 O(K * avg_len)），若 K 很大需注意内存使用。

    while (reader->next(rec)) {
        ++total_records;
        // 对序列进行规范化清洗：例如把字母转为大写，非 AGCTU 替换为 N（具体实现由 seq_io::cleanSequence 提供）。
        // cleanSequence 就地修改 rec.seq，尽量避免重复复制以节省内存带宽。
        seq_io::cleanSequence(rec.seq);

        // 将清洗后的记录写入 clean_data 文件夹
        // 注：FastaWriter::write 已经在内部做了拼接与一次性写出的优化，性能友好。
        clean_writer.write(rec);

        // 将当前记录交给 TopK 选择器进行考虑（内部维护堆以保证 O(log K) 的替换成本）
        // 注意：selector.consider 应复制或接管必要的字段（例如 id/seq），以免后续 rec 被复用/覆盖导致数据错误。
        selector.consider(rec);

        if ((total_records % log_interval) == 0) {
            spdlog::info("Processed {} records so far...", total_records);
        }
    }

    // ---------- 将 TopK 结果写成共识输入文件 ----------
    // 说明：takeSortedDesc 返回按长度降序排序的记录列表（一般用于选取最长的 N 条序列作为共识计算输入）
    seq_io::FastaWriter cons_writer(consensus_file);
    auto cons_seqs = selector.takeSortedDesc();

    // 将选出的序列写到 consensus 文件；注意保持一致的换行宽度等格式规则，FastaWriter 负责这些细节。
    for (const auto& cons_rec : cons_seqs) {
        // 这里写出的 cons_rec 应该是一个深拷贝的 SeqRecord（由 selector 返回以保证安全），如果不是需要在 selector 中做拷贝。
        cons_writer.write(cons_rec);
    }

    // ---------- 统计与返回值 ----------
    const auto t_end = std::chrono::steady_clock::now();
    const double elapsed_s = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();

    spdlog::info("Preprocessing completed. Total records processed: {}. Selected top {} sequences: {}. Elapsed: {:.2f} s",
                 total_records, cons_n, cons_seqs.size(), elapsed_s);
    // 将 size_t 转为项目级别的 uint_t（在 config.hpp 中定义）；防止溢出则截断到 U_MAX
    if (total_records > static_cast<std::size_t>(U_MAX)) {
        spdlog::warn("Processed records ({}) exceed U_MAX ({}); truncating to U_MAX", total_records, U_MAX);
        return static_cast<uint_t>(U_MAX);
    }

    return static_cast<uint_t>(total_records);
}


void alignConsensusSequence(const FilePath& input_file, const FilePath& output_file,
                            const std::string& msa_cmd, const std::string& workdir, int threads)
{

    // 检查输入文件是否存在
    if (!std::filesystem::exists(input_file)) {
        spdlog::warn("Consensus unaligned file not found: {}", input_file.string());
        return;
    }

    // 记录开始时间
    const auto t_start = std::chrono::steady_clock::now();

    spdlog::info("Starting consensus alignment");
    spdlog::info("  input : {}", input_file.string());
    spdlog::info("  output: {}", output_file.string());
    spdlog::info("  tool  : {}", msa_cmd);
    spdlog::info("  thrs  : {}", threads);

    // 尝试记录输入文件大小（若可访问）
    try {
        if (std::filesystem::exists(input_file)) {
            auto in_size = std::filesystem::file_size(input_file);
            spdlog::info("Input file size: {} bytes", in_size);
        }
    } catch (const std::exception &e) {
        spdlog::warn("Failed to stat input file {}: {}", input_file.string(), e.what());
    }

    // 组装命令。默认把 -i / -o / -t 作为参数传入，便于未来统一替换为 cmd 模块的调用。
    cmd::BuildOptions build_opt;
    const std::string cmd_str = cmd::buildCommand(msa_cmd,input_file.string(),output_file.string(),threads, build_opt);
    spdlog::info("Built MSA command (length {}): {}", cmd_str.size(), cmd_str);
    spdlog::info("MSA command (escaped): {}", cmd_str);

    // 调用外部命令（当前使用 std::system；若以后替换为项目内 cmd 接口，只需修改此处）
    try {
        const auto cmd_start = std::chrono::steady_clock::now();
        int rc = cmd::runCommand(cmd_str);
        const auto cmd_end = std::chrono::steady_clock::now();
        const double cmd_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(cmd_end - cmd_start).count();

        if (rc != 0) {
            spdlog::error("MSA command failed (exit code {}): {}", rc, cmd_str);
        } else {
            spdlog::info("MSA command exited with code 0 (success). Elapsed: {:.3f} s", cmd_elapsed);
        }

        // 检查输出文件
        try {
            if (std::filesystem::exists(output_file)) {
                auto out_size = std::filesystem::file_size(output_file);
                spdlog::info("Aligned consensus output exists: {} ({} bytes)", output_file.string(), out_size);
                if (out_size == 0) {
                    spdlog::warn("Aligned consensus output is empty: {}", input_file.string());
                }
            } else {
                spdlog::warn("Aligned consensus output not found after running MSA command: {}", input_file.string());
            }
        } catch (const std::exception &e) {
            spdlog::warn("Failed to stat output file {}: {}", input_file.string(), e.what());
        }

    } catch (const std::exception &e) {
        spdlog::error("Exception while running MSA command for {}: {}", input_file.string(), e.what());
    }

    const auto t_end = std::chrono::steady_clock::now();
    const double elapsed_s = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    spdlog::info("Finished consensus alignment. Total elapsed: {:.3f} s", elapsed_s);
}
