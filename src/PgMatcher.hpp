#ifndef CURC_PGMATCHER_HPP
#define CURC_PGMATCHER_HPP

#include <bits/stdc++.h>
#include <parallel/algorithm>

template <size_t K>
inline std::uint32_t maRushPrime1HashSimplified(const char *str) {
    std::uint64_t hash = K;
    for (std::uint32_t j = 0; j < K/4; ) {
        std::uint32_t k;
        memcpy(&k, str, 4);
        k += j++;
        hash ^= k;
        hash *= 171717;
        str += 4;
    }
    return (std::uint32_t)(hash);
}

struct PgMatchResult {
    uint64_t posSrc;
    uint64_t length;
    uint64_t posDest;

    PgMatchResult(uint64_t posSrc, uint64_t length, uint64_t posDest) : posSrc(posSrc), length(length), posDest(posDest) {}

    bool operator==(const PgMatchResult &rhs) const {
        return posSrc == rhs.posSrc &&
               length == rhs.length &&
               posDest == rhs.posDest;
    }

    bool operator!=(const PgMatchResult &rhs) const {
        return !(rhs == *this);
    }

    bool operator<(const PgMatchResult &rhs) const {
        if (posDest < rhs.posDest)
            return true;
        if (rhs.posDest < posDest)
            return false;
        if (posSrc < rhs.posSrc)
            return true;
        if (rhs.posSrc < posSrc)
            return false;
        return length < rhs.length;
    }
};

/// 注意使用PgMatcher的时候, 先匹配 LQ ref 和 N ref, 最后再匹配 HQ ref, 因为要留到最后一步再匹配删除 HQ ref 的部分.
class PgMatcher {
    static constexpr size_t src_step = 5;
    static constexpr size_t dest_step = 3;
    static constexpr size_t kmer_size = 36;
    static constexpr size_t target_match_len = 50;
    static constexpr size_t bucket_size_limit = 12;
    static char complement_base(char b) {
        switch(b){
            case 'a': return 't';
            case 'A': return 'T';
            case 'c': return 'g';
            case 'C': return 'G';
            case 'g': return 'c';
            case 'G': return 'C';
            case 't': return 'a';
            case 'T': return 'A';
            default: return 'N';
        }
    }

    static void reverse_complement(std::string& seq) {
        size_t L = seq.size();
#pragma omp parallel for
        for (size_t i = 0; i < L / 2; ++i) {
            char a = complement_base(seq[i]);
            char b = complement_base(seq[L - i - 1]);
            seq[i] = b;
            seq[L - i - 1] = a;
        }
        if (L % 2) seq[L / 2] = complement_base(seq[L / 2]);
    }

    static uint32_t get_hash_size(size_t L, size_t step)  {
        constexpr uint8_t HASH_SIZE_MIN_ORDER = 24, HASH_SIZE_MAX_ORDER = 31;
        uint32_t hash_size;
        uint8_t i = HASH_SIZE_MIN_ORDER;
        do {
            hash_size = ((uint32_t) 1) << (i++);
        } while (i <= HASH_SIZE_MAX_ORDER && hash_size < L / step);
        return hash_size;
    }

    static string getTotalMatchStat(uint64_t totalMatchLength, uint64_t destPgLength) {
        return PgSAHelpers::toString(totalMatchLength) + " (" + PgSAHelpers::toString((totalMatchLength * 100.0) / destPgLength, 1) + "%)";
    }

    const std::string& src; // 持有的是引用
    uint32_t hash_size_minus_one;
    std::vector<uint64_t> key_ranges;
    std::vector<uint64_t> pos_array;
    std::vector<PgMatchResult> result;
public:
    explicit PgMatcher(const std::string& src) : src(src) {
        auto hash_size = get_hash_size(src.size(), src_step);
        hash_size_minus_one = hash_size - 1;
        std::vector<uint8_t> counts(hash_size + 2, 0);
        size_t index_count = (src.size() - kmer_size + 1 + src_step - 1) / src_step;
#pragma omp parallel for
        for (size_t i = 0; i < index_count; ++i) {
            uint64_t pos = i * src_step;
            uint32_t hash_pos = maRushPrime1HashSimplified<kmer_size>(src.data() + pos) & hash_size_minus_one;
// #pragma omp critical
            // {
            if (counts[hash_pos] < bucket_size_limit) ++counts[hash_pos];
            // }
        }

        key_ranges.resize(hash_size + 2, 0);
        for (size_t i = 1; i < hash_size + 2; ++i) {
            key_ranges[i] = key_ranges[i - 1] + counts[i - 1];
        }

        uint64_t hash_count = key_ranges[hash_size + 1];
        printf("hash count : %zu \n", hash_count);
        pos_array.resize(hash_count + 2);

#pragma omp parallel for
        for (size_t i = 0; i < index_count; ++i) {
            uint64_t pos = i * src_step;
            uint32_t hash_pos = maRushPrime1HashSimplified<kmer_size>(src.data() + pos) & hash_size_minus_one;
// #pragma omp critical
            // {
            if (counts[hash_pos]) pos_array[key_ranges[hash_pos] + (--counts[hash_pos])] = pos;
            // }
        }
    }

    void match(const std::string& dest, bool destIsSrcPg) {
        const char* src_begin = src.data();
        const char* dest_begin = dest.data();
        const char* src_end = src.data() + src.size();
        const char* dest_end = dest.data() + dest.size();
        size_t match_count = (dest.size() - kmer_size + 1 + dest_step - 1) / dest_step;
#pragma omp parallel for
        for (size_t i = 0; i < match_count; ++i) {
            uint64_t dest_pos = i * dest_step;
            uint32_t hash_pos = maRushPrime1HashSimplified<kmer_size>(dest.data() + dest_pos) & hash_size_minus_one;
            uint64_t begin = key_ranges[hash_pos], end = key_ranges[hash_pos + 1];
            if (begin != end) {
                for (size_t j = begin; j < end; ++j) {
                    uint64_t src_pos = pos_array[j];
                    if (destIsSrcPg && (dest.size() - src_pos < dest_pos)) continue;

                    const char * dest_curr = dest_begin + dest_pos;
                    const char * src_curr = src_begin + src_pos;
                    const char* p1 = dest_curr + kmer_size - 1;
                    const char* p2 = src_curr + kmer_size - 1;
                    while (++p1 != dest_end && ++p2 != src_end && *p1 == *p2);
                    const char* right = p1;
                    p1 = dest_curr;
                    p2 = src_curr;
                    while (p1 != dest_begin && p2 != src_begin && *--p1 == *--p2);
                    if (right - p1 > target_match_len && memcmp(dest_curr, src_curr, kmer_size) == 0) {
#pragma omp critical
                        result.emplace_back(p2 + 1 - src_begin, right - p1 - 1, p1 + 1 - dest_begin);
                        break;
                    }
                }
            }
        }
    }

    void markAndRemoveExactMatches(bool destIsSrcPg, std::string& destPg, std::string& resPgMapOff, std::string& resPgMapLen) {
        uint64_t destPgLength = destPg.size();

        if (destIsSrcPg) {
            std::string dest = destPg;
            reverse_complement(dest);
            match(dest, destIsSrcPg);
        } else {
            reverse_complement(destPg);
            match(destPg, destIsSrcPg);
            reverse_complement(destPg);
        }

        printf("match result size: %zu\n", result.size());

        // correct dest pos due to rev complement matching
#pragma omp parallel for
        for (size_t i = 0; i < result.size(); ++i) {
            result[i].posDest = destPg.size() - (result[i].posDest + result[i].length);
        }

        // resolve mapping collision in the same text
        if (destIsSrcPg) {
#pragma omp parallel for
            for (size_t i = 0; i < result.size(); ++i) {
                if (result[i].posSrc > result[i].posDest) {
                    uint64_t tmp = result[i].posSrc;
                    result[i].posSrc = result[i].posDest;
                    result[i].posDest = tmp;
                }
                uint64_t endPosSrc = (result[i].posSrc + result[i].length);
                if (endPosSrc > result[i].posDest) {
                    uint64_t margin = (endPosSrc - result[i].posDest + 1) / 2;
                    result[i].length -= margin;
                    result[i].posDest += margin;
                }
            }
        }

        std::ostringstream pgMapOffDest;
        std::ostringstream pgMapLenDest;
        PgSAHelpers::writeUIntByteFrugal(pgMapLenDest, target_match_len);
        __gnu_parallel::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
        printf("Unique exact matches: %zu\n", result.size());

        char * dest_ptr = &destPg[0]; // destPg.data(); // gcc >= 7.5 才支持返回非const的data pointer
        uint64_t pos = 0, npos = 0, totalDestOverlap = 0, totalMatched = 0;
        bool isPgLengthStd = src.size() <= UINT32_MAX; /// 注意后续文件格式也要写入 isPgLengthStd.
        for (auto & match : result) {
            if (match.posDest < pos) {
                uint64_t overflow = pos - match.posDest;
                if (overflow >= match.length) {
                    totalDestOverlap += match.length;
                    match.length = 0;
                    continue;
                }
                totalDestOverlap += overflow;
                match.length -= overflow;
                match.posDest += overflow;
            }
            if (match.length < target_match_len) {
                totalDestOverlap += match.length;
                continue;
            }
            totalMatched += match.length;
            uint64_t length = match.posDest - pos;
            std::memmove(dest_ptr + npos, dest_ptr + pos, length);
            npos += length;
            destPg[npos++] = '%'; // mark
            if (isPgLengthStd) {
                PgSAHelpers::writeValue<uint32_t>(pgMapOffDest, match.posSrc);
            } else {
                PgSAHelpers::writeValue<uint64_t>(pgMapOffDest, match.posSrc);
            }
            PgSAHelpers::writeUIntByteFrugal(pgMapLenDest, match.length - target_match_len);
            pos = match.posDest + match.length;
        }
        uint64_t length = destPg.size() - pos;
        std::memmove(dest_ptr + npos, dest_ptr + pos, length);
        npos += length;
        destPg.resize(npos);

        result.clear();
        resPgMapOff = pgMapOffDest.str();
        pgMapOffDest.clear();
        resPgMapLen = pgMapLenDest.str();
        pgMapLenDest.clear();

        printf("Final size of Pg: %zu (remove: %s ; %zu chars in overlapped dest symbol)\n",
               npos, getTotalMatchStat(totalMatched, destPgLength).c_str(), totalDestOverlap);

    }
};

#endif //CURC_PGMATCHER_HPP
