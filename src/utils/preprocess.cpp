#include "preprocess.h"
#include <chrono>

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

    const auto t_start = std::chrono::steady_clock::now();

    spdlog::info("Preprocessing input FASTA file: {}", input_path);
    spdlog::info("Working directory: {}", workdir);

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

    // 4) 将输入文件复制或下载到 raw_data 下。
    //    这里调用 file_io::fetchFile，函数内部会判断是 URL 还是本地路径并做相应操作（download 或 copy）。
    FilePath input_file = FilePath(input_path);
    FilePath raw_dest_file = raw_data_dir / input_file.filename();

    spdlog::info("Fetching input to working raw path: {} -> {}", input_file.string(), raw_dest_file.string());
    // 说明：fetchFile 在遇到远程 URL 时会调用 downloadFile 将数据写入本地；在本地路径时会调用 copyFile（包含跨设备回退等逻辑）。
    file_io::fetchFile(input_file,raw_dest_file);
    spdlog::info("Input available at: {}", raw_dest_file.string());

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
    FilePath consensus_file = clean_data_dir / CLEAN_UNALIGNED;
    spdlog::info("Clean output: {} ; Consensus output: {}", clean_dest_file.string(), consensus_file.string());

    // seq_io::openKseqReader 返回一个 reader 指针（抽象），用于逐条读取序列；
    // seq_io::FastaWriter 用于把清洗后的序列写入到目标文件。
    auto reader = seq_io::openKseqReader(raw_dest_file);
    seq_io::FastaWriter clean_writer(clean_dest_file);
    TopKLongestSelector selector(cons_n);

    seq_io::SeqRecord rec;
    std::size_t total_records = 0;
    const std::size_t log_interval = 10000;
    while (reader->next(rec)) {
        ++total_records;
        // 对序列进行规范化清洗：例如把字母转为大写，非 AGCTU 替换为 N（具体实现由 seq_io::cleanSequence 提供）。
        seq_io::cleanSequence(rec.seq);
        // 将清洗后的记录写入 clean_data 文件夹
        clean_writer.write(rec);
        // 将当前记录交给 TopK 选择器进行考虑（内部维护堆以保证 O(log K) 的替换成本）
        selector.consider(rec);

        if ((total_records % log_interval) == 0) {
            spdlog::info("Processed {} records so far...", total_records);
        }
    }

    // 把选出的序列写入到文件中
    seq_io::FastaWriter cons_writer(consensus_file);
    auto cons_seqs = selector.takeSortedDesc();
    for (const auto& cons_rec : cons_seqs) {
        cons_writer.write(cons_rec);
    }

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

// ---------------- TopKLongestSelector 实现说明 ----------------
// TopKLongestSelector 用于在单次顺序扫描中选出长度最长的 K 条序列。
// 设计目标：
// - 空间复杂度 O(K)，只保留 K 条记录（使用堆存储）；
// - 时间复杂度：每个元素的插入/维护为 O(log K)；遍历 N 条记录总体为 O(N log K)；
// - 稳定性：对相同长度的序列，优先保留较早出现的序列（即出现更早的序列被视为更“好”）。
//
// 实现细节：
// - 使用一个最小堆（min-heap），heap_[0] 始终为当前已收集 K 条中的最差（即最短或最晚出现的）元素；
// - 当堆未满（size < K）时直接加入并上浮（siftUp）；
// - 当堆已满时，比较候选项与 heap_[0]：如果候选项更好（长度更长或相同长度但出现更早），则替换 heap_[0] 并下沉（siftDown）；
// - 比较规则（betterThan / worseThan）使用 (length, order) 的字典序：长度为主，出现顺序为次（更早更优）。

TopKLongestSelector::TopKLongestSelector(std::size_t k)
        : k_(k)
    {
        heap_.reserve(k_);
    }

    void TopKLongestSelector::reset(std::size_t k)
    {
        k_ = k;
        order_counter_ = 0; // order_counter_ 用于记录元素出现的先后顺序，以便在长度相同时保持稳定性
        heap_.clear();
        heap_.reserve(k_);
    }

    std::size_t TopKLongestSelector::size() const
    {
        return heap_.size();
    }

    std::size_t TopKLongestSelector::capacity() const
    {
        return k_;
    }

    bool TopKLongestSelector::empty() const
    {
        return heap_.empty();
    }

    // 比较函数：判断 a 是否比 b 更 "差"
    // 语义：在最小堆中，“更差”的元素会被置于堆顶（即应当被替换或弹出）。
    bool TopKLongestSelector::worseThan(const Item& a, const Item& b)
    {
        if (a.len != b.len) return a.len < b.len;      // 长度更短的被认为更差
        return a.order > b.order;                      // 同长度时，后来出现的（order 值更大）被认为更差
    }

    // 判断候选项 cand 是否比当前最坏项 worst 更好（用于决定是否替换堆顶）
    bool TopKLongestSelector::betterThan(const Item& cand, const Item& worst)
    {
        if (cand.len != worst.len) return cand.len > worst.len; // 更长则更好
        // 同长度时，我们偏向保留更早出现的（即 order 更小的项优先）
        return cand.order < worst.order;
    }

    // 上浮操作：当新元素插入到堆尾时调用，若其比父节点更差则与父节点交换直到堆序恢复
    // 说明：这里实现的是一个 min-heap 的上浮（维护堆顶为最差元素），使得堆顶始终为当前最差
    void TopKLongestSelector::siftUp(std::size_t idx)
    {
        while (idx > 0) {
            const std::size_t parent = (idx - 1) / 2;
            // min-heap：如果当前更“差”（更小），就上浮（与父节点交换）
            if (worseThan(heap_[idx], heap_[parent])) {
                std::swap(heap_[idx], heap_[parent]);
                idx = parent;
            } else {
                break;
            }
        }
    }

    // 下沉操作：当堆顶被替换后调用，把新的堆顶下沉到合适位置
    // 通过比较左右孩子找到更"差"的孩子并与当前交换，直到堆序恢复
    void TopKLongestSelector::siftDown(std::size_t idx)
    {
        const std::size_t n = heap_.size();
        while (true) {
            const std::size_t left = idx * 2 + 1;
            if (left >= n) break; // 无孩子，停止

            const std::size_t right = left + 1;

            // 选择更“差”的孩子（min-heap 中更差的孩子应当成为交换对象）
            std::size_t worst_child = left;
            if (right < n && worseThan(heap_[right], heap_[left])) {
                worst_child = right;
            }

            // 如果孩子比当前更“差”，则下沉（交换）
            if (worseThan(heap_[worst_child], heap_[idx])) {
                std::swap(heap_[idx], heap_[worst_child]);
                idx = worst_child;
            } else {
                break;
            }
        }
    }

    // 考虑一条新记录：将其包装为 Item 并尝试加入堆中
    // 逻辑：
    // - 如果 k_ == 0，直接返回（不保存任何记录）
    // - 如果堆未满，则直接 push 并上浮
    // - 否则比较候选与堆顶（最差）；若候选更好，则替换堆顶并下沉
    void TopKLongestSelector::consider(seq_io::SeqRecord rec)
    {
        if (k_ == 0) return;

        Item cand;
        cand.len = rec.seq.size();
        cand.order = order_counter_++;
        cand.rec = std::move(rec);

        if (heap_.size() < k_) {
            heap_.push_back(std::move(cand));
            siftUp(heap_.size() - 1);
            return;
        }

        // heap_[0] 是当前“最差”的那条（堆顶）
        if (betterThan(cand, heap_[0])) {
            heap_[0] = std::move(cand);
            siftDown(0);
        }
    }

    // 将保留的堆内容按长度降序输出（同长度按出现顺序升序，保证稳定性）
    // 返回值：vector<SeqRecord>，其中包含已选择的记录（按长度从长到短排序）
    std::vector<seq_io::SeqRecord> TopKLongestSelector::takeSortedDesc()
    {
        // 为避免复制大量序列数据，move 出 heap_ 的内容到一个临时 vector 上进行排序
        std::vector<Item> items = std::move(heap_);
        heap_.clear();
        heap_.shrink_to_fit();
        heap_.reserve(k_);

        // 自定义排序：长度降序；若长度相同则按原始出现顺序升序（更早的排在前面）
        std::sort(items.begin(), items.end(),
                  [](const Item& a, const Item& b) {
                      if (a.len != b.len) return a.len > b.len;
                      return a.order < b.order;
                  });

        // 将 SeqRecord 从 items 移出到返回向量中
        std::vector<seq_io::SeqRecord> out;
        out.reserve(items.size());
        for (auto& it : items) {
            out.push_back(std::move(it.rec));
        }
        return out;
    }

