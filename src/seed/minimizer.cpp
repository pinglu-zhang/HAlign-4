#include "seed.h"

#include <cstdint>
#include <string>
#include <vector>

namespace minimizer
{
    // splitmix64：非常常用的 64-bit mixer，速度快且分布不错。
    // 在 minimap2 里也会对 k-mer 编码做 hash64 混洗，理念相同。
    static inline constexpr std::uint64_t splitmix64(std::uint64_t x) noexcept
    {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    // Cand：窗口中的候选 (hash, pos)
    // pos 是 k-mer 的起始位置（0-based）
    struct Cand
    {
        std::uint64_t h;
        std::uint32_t pos;
    };

    // =========================================================
    // RingMinQueue：固定容量的“环形单调队列”（维护窗口最小值）
    // =========================================================
    class RingMinQueue
    {
    public:
        explicit RingMinQueue(std::uint32_t capacity)
            : buf_(capacity), cap_(capacity)
        {
        }

        void clear() noexcept
        {
            head_ = 0;
            size_ = 0;
        }

        bool empty() const noexcept { return size_ == 0; }

        void push(std::uint64_t h, std::uint32_t pos) noexcept
        {
            const Cand c{h, pos};
            // 弹出所有 >= 当前 hash 的队尾元素，保持单调递增。
            while (size_ && back().h >= c.h) pop_back();
            buf_[idx(size_)] = c;
            ++size_;
        }

        void popExpired(std::uint32_t win_start) noexcept
        {
            while (size_ && front().pos < win_start) pop_front();
        }

        std::uint64_t minHash() const noexcept { return front().h; }

    private:
        std::vector<Cand> buf_;
        std::uint32_t cap_{0};
        std::uint32_t head_{0};
        std::uint32_t size_{0};

        std::uint32_t idx(std::uint32_t off) const noexcept
        {
            std::uint32_t i = head_ + off;
            if (i >= cap_) i -= cap_;
            return i;
        }

        Cand& front() noexcept { return buf_[head_]; }
        const Cand& front() const noexcept { return buf_[head_]; }

        Cand& back() noexcept { return buf_[idx(size_ - 1)]; }
        const Cand& back() const noexcept { return buf_[idx(size_ - 1)]; }

        void pop_front() noexcept
        {
            head_ = idx(1);
            --size_;
        }

        void pop_back() noexcept
        {
            --size_;
        }
    };


    // =============================================================
    // 目标：从单条序列中提取 minimizer（hash-only）。
    //
    // 设计思路（尽量贴近 minimap2 的 sketch 过程）：
    // 1) rolling 2-bit 编码生成每个位置的 k-mer（O(n)）。
    // 2) 对 k-mer 编码做 64-bit mixer（splitmix64）。
    // 3) winnowing：每个窗口(w 个 k-mer)取 hash 最小者作为 minimizer。
    // 4) 相邻窗口去重。
    // =============================================================
    MinimizerHashes extractMinimizerHash(const std::string& seq,
                                         std::size_t k,
                                         std::size_t w,
                                         bool is_forward)
    {
        MinimizerHashes out;

        const std::uint32_t n = static_cast<std::uint32_t>(seq.size());
        if (k == 0 || w == 0 || n < k) return out;
        if (k > 31) return out;
        if (w >= 256) return out;

        const auto& nt4_table_ref = minimizer::nt4_table;

        const std::uint32_t total_kmer = n - static_cast<std::uint32_t>(k) + 1;
        const std::uint32_t win = std::min<std::uint32_t>(static_cast<std::uint32_t>(w), total_kmer);
        if (win == 0) return out;

        out.reserve(std::max<std::uint32_t>(1, n / win));

        const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
        const std::uint64_t shift = 2ULL * (k - 1);

        std::uint64_t fwd = 0;
        std::uint64_t rev = 0;
        std::uint32_t valid = 0;

        RingMinQueue q(win);

        std::uint64_t last_out = 0;
        bool has_last = false;

        for (std::uint32_t i = 0; i < n; ++i) {
            // 1) 查表将字符转为 2-bit（非法字符 -> 4）
            const std::uint8_t c = nt4_table_ref[static_cast<std::uint8_t>(seq[i])];
            if (c >= 4) {
                // 非法字符会打断 k-mer：清空滚动状态与窗口结构
                fwd = rev = 0;
                valid = 0;
                q.clear();
                continue;
            }

            // 2) rolling 更新 fwd/rev
            fwd = ((fwd << 2) | c) & mask;
            rev = (rev >> 2) | (std::uint64_t(3U ^ c) << shift);

            // 3) 只有连续合法字符达到 k 后，才能形成 k-mer
            if (valid < k) ++valid;
            if (valid < k) continue;

            const std::uint32_t pos = i + 1 - static_cast<std::uint32_t>(k);

            // 4) 根据 is_forward 选择正向或反向互补编码
            const std::uint64_t code = is_forward ? fwd : rev;
            const std::uint64_t h = splitmix64(code);

            // 5) 放入窗口最小值结构
            q.push(h, pos);

            // 6) 窗口满 win 个 k-mer 后才输出 minimizer
            if (pos + 1 < win) continue;
            const std::uint32_t win_start = pos + 1 - win;

            // 7) 弹出过期元素
            q.popExpired(win_start);

            // 8) 输出窗口最小值，并做相邻去重
            if (!q.empty()) {
                const std::uint64_t mh = q.minHash();
                if (!has_last || mh != last_out) {
                    out.push_back(MinimizerHash{mh});
                    last_out = mh;
                    has_last = true;
                }
            }
        }

        return out;
    }

} // namespace minimizer

