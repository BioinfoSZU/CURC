#ifndef CURC_UTIL_CUH
#define CURC_UTIL_CUH

#include <bits/stdc++.h>
#include <experimental/filesystem>
#include <parallel/algorithm>
#ifndef ABS
#define ABS std::abs
#endif
#include <fqzcomp/clr.cdr>
#include <fqzcomp/simple_model.h>
#include <lzma/VarLenDNACoder/LzmaLib.h>
#include "fast_lzma2_helper.hpp"

namespace fs = std::experimental::filesystem;

static inline void gpu_assert(cudaError_t code, const char *file, int line) {
    if (code != cudaSuccess) {
        fprintf(stderr,"GPU assert: %s %s %d\n", cudaGetErrorString(code), file, line);
        exit(code);
    }
}
#define gpuErrorCheck(ans) { gpu_assert((ans), __FILE__, __LINE__); }

template<size_t unit_size>
static inline size_t ceil(size_t value) {
    return (value + unit_size - 1) / unit_size;
}

static inline uint32_t get_hash_size(size_t item_count) {
    constexpr uint8_t HASH_SIZE_MIN_ORDER = 24, HASH_SIZE_MAX_ORDER = 31;
    uint32_t hash_size;
    uint8_t i = HASH_SIZE_MIN_ORDER;
    do {
        hash_size = ((uint32_t) 1) << (i++);
    } while (i <= HASH_SIZE_MAX_ORDER && hash_size < item_count);
    return hash_size;
}

__device__ __forceinline__ uint64_t murmur_hash64(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccd;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53;
    x ^= x >> 33;
    return x;
}

/** ----------------- CPU util function --------------------- */
static constexpr char base_decode_table_cpu[] = "ACGTN-acgtn*";
static constexpr uint8_t base_encode_table_cpu[256] = { // A : 0; C : 1; G : 2; T : 3
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
static inline uint64_t extract_word_cpu(const uint64_t * bits, uint32_t start_index, uint64_t off, uint8_t base_number) {
    uint64_t offset1 = off >> 5ULL;
    uint32_t offset2 = (off & 0x1FULL) << 1ULL;
    uint64_t left_part = bits[start_index + offset1] << offset2;
    if ((64ULL - offset2) >= (base_number << 1ULL)) return left_part >> ((32ULL - base_number) << 1ULL);
    uint64_t right_part = (bits[start_index + offset1 + 1] >> (62ULL - offset2)) >> 2ULL;
    uint64_t kmer_right_shift = (32ULL - base_number) << 1ULL;
    return (left_part | right_part) >> kmer_right_shift;
}
static inline uint64_t extract_word_unaligned_cpu(const uint64_t * bits, uint64_t base_offset, uint8_t base_number){
    uint64_t offset1 = base_offset >> 5ULL;
    uint32_t offset2 = (base_offset & 0x1FULL) << 1ULL;
    uint64_t left_part = bits[offset1] << offset2;
    uint64_t right_part = (bits[offset1 + 1] >> (62ULL - offset2)) >> 2ULL;
    uint64_t kmer_right_shift = (32ULL - base_number) << 1ULL;
    return (left_part | right_part) >> kmer_right_shift;
}
static inline char extract_base_cpu(const uint64_t * reads_db, uint32_t read_start_index, uint64_t off) {
    return base_decode_table_cpu[(reads_db[read_start_index + (off >> 5ULL)] >> (((~off) & 0x1FULL) << 1ULL)) & 0x03ULL];
}
static inline char reverse_dna_base_cpu(char b) {
    switch(b) {
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


/** ----------------- GPU util function ------------------- */
__device__ __constant__ char base_decode_table_mgpu[] = "ACGTN-acgtn*";
__device__ __constant__ uint8_t base_encode_table_mgpu[256] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
__device__ __forceinline__ char reverse_dna_base(char b){
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

__device__ __forceinline__ uint64_t extract_word_gpu(const uint64_t * bits, uint32_t start_index, uint64_t off, uint8_t base_number) {
    uint64_t offset1 = off >> 5ULL;
    uint32_t offset2 = (off & 0x1FULL) << 1ULL;
    uint64_t left_part = bits[start_index + offset1] << offset2;
    if ((64ULL - offset2) >= (base_number << 1ULL)) return left_part >> ((32ULL - base_number) << 1ULL);
    uint64_t right_part = (bits[start_index + offset1 + 1] >> (62ULL - offset2)) >> 2ULL;
    uint64_t kmer_right_shift = (32ULL - base_number) << 1ULL;
    return (left_part | right_part) >> kmer_right_shift;
}
__device__ __forceinline__ uint64_t extract_word_unaligned_gpu(const uint64_t * bits, uint64_t base_offset, uint8_t base_number){
    uint64_t offset1 = base_offset >> 5ULL;
    uint32_t offset2 = (base_offset & 0x1FULL) << 1ULL;
    uint64_t left_part = bits[offset1] << offset2;
    uint64_t right_part = (bits[offset1 + 1] >> (62ULL - offset2)) >> 2ULL;
    uint64_t kmer_right_shift = (32ULL - base_number) << 1ULL;
    return (left_part | right_part) >> kmer_right_shift;
}
__device__ __forceinline__ char extract_base_gpu(const uint64_t * reads_db, uint32_t read_start_index, uint64_t off) {
    return base_decode_table_mgpu[(reads_db[read_start_index + (off >> 5ULL)] >> (((~off) & 0x1FULL) << 1ULL)) & 0x03ULL];
}
__device__ __forceinline__ char extract_base_unaligned_gpu(const uint64_t * bits, uint64_t base_offset) {
    uint64_t offset1 = base_offset >> 5ULL;
    uint32_t offset2 = (((~(base_offset)) & 0x1FULL) << 1ULL);
    return base_decode_table_mgpu[(bits[offset1] >> offset2) & 0x03ULL];
}

__device__ __forceinline__ uint64_t byte_swap_64_cuda(uint64_t value) {
    uint32_t left = value >> 32ULL;
    uint32_t right = value;
    uint32_t left_byte_swap  = __byte_perm (left, 0, 0x0123);  // 32bit byte swap
    uint32_t right_byte_swap = __byte_perm (right, 0, 0x0123); // 32bit byte swap
    uint64_t left_64 = left_byte_swap, right_64 = right_byte_swap;
    return left_64 | (right_64 << 32ULL);
}
__global__ void test_byte_swap_64_cuda () {
    uint64_t value = 0x1233456678992311;
    printf ("input:  %08lx\n", value);
    value = byte_swap_64_cuda(value);
    printf ("output : %08lx\n", value);
}
__device__ __forceinline__ uint64_t get_reverse_kmer(uint64_t kmer, uint8_t kmer_size) {
    kmer = ~kmer;
    kmer = ((kmer & 0x3333333333333333ULL) << 2ULL) | ((kmer & 0xCCCCCCCCCCCCCCCCULL) >> 2ULL);
    kmer = ((kmer & 0x0F0F0F0F0F0F0F0FULL) << 4ULL) | ((kmer & 0xF0F0F0F0F0F0F0F0ULL) >> 4ULL);
    kmer = byte_swap_64_cuda(kmer);
    return (kmer >> ((32ULL - kmer_size) << 1ULL));
}

__device__ __forceinline__ int mismatch_count_in_two_kmers(uint64_t kmer1, uint64_t kmer2) {
    uint64_t kmer_xor = kmer1 ^ kmer2;
    uint64_t ret = ((kmer_xor & 0xAAAAAAAAAAAAAAAAULL) >> 1ULL) | (kmer_xor & 0x5555555555555555ULL);
    return __popcll(ret);
}
__device__ __forceinline__ int forward_mismatch_count(const uint64_t * data1, uint64_t offset1, // ref
                                                      const uint64_t * data2, uint64_t offset2, // reads_db
                                                      uint32_t length, uint32_t max_mismatch_count) {
    int mismatch_count = 0;
    uint64_t kmer_32_a, kmer_32_b;
    size_t i = 0;
    for (i = 0; (i + 32 < length) && (mismatch_count <= max_mismatch_count); i += 32) {
        kmer_32_a = extract_word_unaligned_gpu(data1, offset1 + i, 32);
        kmer_32_b = extract_word_gpu(data2, offset2, i, 32);
        mismatch_count += mismatch_count_in_two_kmers(kmer_32_a, kmer_32_b);
    }
    if ((i < length) && (mismatch_count <= max_mismatch_count)) {
        kmer_32_a = extract_word_unaligned_gpu(data1, offset1 + i, length - i);
        kmer_32_b = extract_word_gpu(data2, offset2, i, length - i);
        mismatch_count += mismatch_count_in_two_kmers(kmer_32_a, kmer_32_b);
    }
    return mismatch_count;
}
__device__ __forceinline__ int reverse_complement_mismatch_count(const uint64_t * data1, uint64_t offset1, // ref
                                                                 const uint64_t * data2, uint64_t offset2, // reads_db
                                                                 uint32_t length, uint32_t max_mismatch_count) {
    int mismatch_count = 0;
    uint64_t kmer_32_a, kmer_32_b;
    size_t i = 0;
    for (i = 0; (i + 32 < length) && (mismatch_count <= max_mismatch_count); i += 32) {
        kmer_32_a = extract_word_unaligned_gpu(data1, offset1 + i, 32);
        kmer_32_b = extract_word_gpu(data2, offset2, length - (i + 32), 32);
        kmer_32_b = get_reverse_kmer(kmer_32_b, 32);
        mismatch_count += mismatch_count_in_two_kmers(kmer_32_a, kmer_32_b);
    }
    if ((i < length) && (mismatch_count <= max_mismatch_count)) {
        kmer_32_a = extract_word_unaligned_gpu(data1, offset1 + i, length - i);
        kmer_32_b = extract_word_gpu(data2, offset2, 0, length - i);
        kmer_32_b = get_reverse_kmer(kmer_32_b, length - i);
        mismatch_count += mismatch_count_in_two_kmers(kmer_32_a, kmer_32_b);
    }
    return mismatch_count;
}
__device__ __forceinline__ uint64_t get_kmer_integer (const char* data, uint64_t base_number) {
    uint64_t bits = 0, c;
    uint8_t idx2 = 0;
    size_t kmer_bits = base_number * 2;
    for (size_t i = 0; i < base_number; ++i) {
        c = base_encode_table_mgpu[data[i]] & 0x03U;
        bits |= (c << (kmer_bits - 2 - idx2));
        idx2 = (idx2 + 2U) & 0x3FU;
    }
    return bits;
}

template<typename T>
static inline void write_data_to_binary_file(const T *data, size_t data_size, const char *filename) {
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::binary);
    if (data_size > 0) {
        file.write(reinterpret_cast<const char *>(data), sizeof(T) * data_size);
    }
}
template<typename T>
static inline void write_vector_to_binary_file(const std::vector<T> &data, const std::string &filename) {
    write_data_to_binary_file(data.data(), data.size(), filename.c_str());
}
template<typename T>
static inline void read_vector_from_binary_file(std::vector<T> &data, const fs::path &filename) {
    if (!fs::exists(filename)) throw std::runtime_error(filename.string() + " isn't exist");
    std::ifstream file;
    file.open(filename, std::ios::in | std::ios::binary);
    file.seekg(0, std::ios::end);
    auto size = file.tellg();
    file.seekg(0);
    data.resize(size / sizeof(T));
    file.read(reinterpret_cast<char *>(data.data()), size);
}

template<size_t max_read_len>
void mismatch_offset_compress(const std::vector<uint16_t> & mismatch_offset_array,
                              size_t mismatch_count,
                              const std::string& name,
                              const std::string& suffix = ".bin") {
    std::vector<char> compressed_output;
    compressed_output.resize(mismatch_offset_array.size() * sizeof(uint16_t) + 1024);
    char * output = compressed_output.data();

    RangeCoder rc{};
    rc.output(output);
    rc.StartEncode();

    SIMPLE_MODEL<2> flag;
    SIMPLE_MODEL<max_read_len / 2 + 1> first_value;
    SIMPLE_MODEL<max_read_len> value;

    size_t i = 0;
    while (i < mismatch_offset_array.size()) {
        flag.encodeSymbol(&rc, mismatch_offset_array[i++]);
        first_value.encodeSymbol(&rc, mismatch_offset_array[i++]);
        if (mismatch_count > 1) {
            if (mismatch_count == 2) {
                value.encodeSymbol(&rc, mismatch_offset_array[i++]);
            } else {
                value.encodeSymbol(&rc, mismatch_offset_array[i++]);
                for (size_t m = 1; m < mismatch_count; ++m) {
                    value.encodeSymbol(&rc, mismatch_offset_array[i++]);
                }
            }
        }
    }

    rc.FinishEncode();
    compressed_output.resize(rc.size_out());
    write_vector_to_binary_file(compressed_output, name + suffix);
}

template<typename RefLenType>
void paired_end_id_to_pos_compression(const std::vector<uint64_t>& id_to_pos, std::ofstream& out, const fs::path& working_path, int level, int t_num) {
    auto read_file_and_output = [&](const fs::path &filename) {
        std::vector<uint8_t> data;
        read_vector_from_binary_file(data, filename);
        uint64_t size = data.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof (uint64_t));
        PgSAHelpers::writeArray(out, data.data(), size);
    };

    vector<RefLenType> basePairPos;
    vector<uint8_t> offsetInUint16Flag;
    vector<uint8_t> offsetIsBaseFirstFlag;
    vector<uint16_t> offsetInUint16Value;
    vector<uint8_t> deltaInInt16Flag;
    vector<uint8_t> deltaIsBaseFirstFlag;
    vector<int16_t> deltaInInt16Value;
    vector<RefLenType> notBasePairPos;
    const uint32_t pairsCount = id_to_pos.size() / 2;
    basePairPos.reserve(pairsCount);
    offsetInUint16Flag.reserve(pairsCount);
    offsetIsBaseFirstFlag.reserve(pairsCount);
    offsetInUint16Value.reserve(pairsCount);
    deltaInInt16Flag.reserve(pairsCount / 2);
    deltaIsBaseFirstFlag.reserve(pairsCount / 8);
    deltaInInt16Value.reserve(pairsCount / 8);
    notBasePairPos.reserve(pairsCount / 4);

    vector<uint32_t> bppRank;
    bppRank.reserve(pairsCount);
    for (uint32_t i = 0; i < id_to_pos.size(); i += 2) {
        basePairPos.push_back(id_to_pos[i]);
        bppRank.push_back(i >> 1);
    }
    __gnu_parallel::stable_sort(bppRank.begin(), bppRank.end(),
                                [&](const uint32_t &idx1, const uint32_t &idx2) -> bool
                                { return basePairPos[idx1] < basePairPos[idx2]; });
    int64_t refPrev = 0;
    int64_t prev = 0;
    bool match = false;
    for (uint32_t p = 0; p < pairsCount; p++) {
        uint32_t i = bppRank[p] * 2;

        bool isBaseBefore = id_to_pos[i] < id_to_pos[i + 1];
        RefLenType relativeAbsOffset = isBaseBefore?(id_to_pos[i + 1] - id_to_pos[i]):
                                        id_to_pos[i] - id_to_pos[i + 1];
        const bool isOffsetInUint16 = relativeAbsOffset <= UINT16_MAX;
        offsetInUint16Flag.push_back(isOffsetInUint16 ? 1 : 0);
        if (isOffsetInUint16) {
            offsetIsBaseFirstFlag.push_back(isBaseBefore?1:0);
            offsetInUint16Value.push_back((uint16_t) relativeAbsOffset);
            continue;
        }
        {
            const int64_t delta = relativeAbsOffset - refPrev;
            const bool isDeltaInInt16 = delta <= INT16_MAX && delta >= INT16_MIN;
            deltaInInt16Flag.push_back((uint8_t) isDeltaInInt16);
            if (isDeltaInInt16) {
                match = true;
                deltaIsBaseFirstFlag.push_back(isBaseBefore ? 1 : 0);
                deltaInInt16Value.push_back((int16_t) delta);
                refPrev = relativeAbsOffset;
            } else {
                if (!match || refPrev != prev)
                    refPrev = relativeAbsOffset;
                notBasePairPos.push_back((RefLenType) id_to_pos[i + 1]);
                match = false;
            }
            prev = relativeAbsOffset;
        }
    }

//    uint32_t total_reads_count = id_to_pos.size();
//    out.write(reinterpret_cast<const char*>(&total_reads_count), 4);

    std::ofstream base_pair_pos_file(working_path / "base_pair_pos.bin");
    PgSAHelpers::writeArray(base_pair_pos_file, basePairPos.data(), basePairPos.size() * sizeof(RefLenType));
    lzma2::lzma2_compress((working_path / "base_pair_pos.bin").c_str(), (working_path / "base_pair_pos.lzma2").c_str(), level, t_num);
    read_file_and_output(working_path / "base_pair_pos.lzma2");

    writeCompressed(out, (char *) offsetInUint16Flag.data(), offsetInUint16Flag.size() * sizeof(uint8_t),
                    PPMD7_CODER, 3, 3, COMPRESSION_ESTIMATION_UINT8_BITMAP);
    writeCompressed(out, (char *) offsetIsBaseFirstFlag.data(), offsetIsBaseFirstFlag.size() * sizeof(uint8_t),
                    PPMD7_CODER, 3, 3, COMPRESSION_ESTIMATION_UINT8_BITMAP);
    writeCompressed(out, (char*) offsetInUint16Value.data(), offsetInUint16Value.size() * sizeof(uint16_t),
                    PPMD7_CODER, 3, 3, 1);
    writeCompressed(out, (char *) deltaInInt16Flag.data(), deltaInInt16Flag.size() * sizeof(uint8_t),
                    PPMD7_CODER, 3, 3, COMPRESSION_ESTIMATION_UINT8_BITMAP);
    writeCompressed(out, (char *) deltaIsBaseFirstFlag.data(), deltaIsBaseFirstFlag.size() * sizeof(uint8_t),
                    PPMD7_CODER, 3, 3, COMPRESSION_ESTIMATION_UINT8_BITMAP);
    writeCompressed(out, (char *) deltaInInt16Value.data(), deltaInInt16Value.size() * sizeof(int16_t),
                    PPMD7_CODER, 3, 3, 1);

    std::ofstream not_base_pair_pos_file(working_path / "not_base_pair_pos.bin");
    PgSAHelpers::writeArray(not_base_pair_pos_file, notBasePairPos.data(), notBasePairPos.size() * sizeof(RefLenType));
    lzma2::lzma2_compress((working_path / "not_base_pair_pos.bin").c_str(), (working_path / "not_base_pair_pos.lzma2").c_str(), level, t_num);
    read_file_and_output(working_path / "not_base_pair_pos.lzma2");
}

void paired_end_reads_order_compress(const std::vector<uint32_t>& pe_reads_order, std::ofstream& out, const fs::path& working_path, int level, int t_num) {
    auto read_file_and_output = [&](const fs::path &filename) {
        std::vector<uint8_t> data;
        read_vector_from_binary_file(data, filename);
        uint64_t size = data.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof (uint64_t));
        PgSAHelpers::writeArray(out, data.data(), size);
    };

    std::vector<uint8_t> offsetPairFlag;
    std::vector<uint8_t> nonOffsetPairFlag;
    std::vector<uint8_t> offsetInUint8Flag;
    std::vector<uint8_t> offsetInUint8Value;
    std::vector<uint8_t> deltaInInt8Flag;
    std::vector<int8_t>  deltaInInt8Value;
    std::vector<uint32_t> offset;
    offsetPairFlag.reserve(pe_reads_order.size() / 2);
    nonOffsetPairFlag.reserve(pe_reads_order.size() / 4);
    offsetInUint8Flag.reserve(pe_reads_order.size() / 2);
    offsetInUint8Value.reserve(pe_reads_order.size() / 2);
    deltaInInt8Flag.reserve(pe_reads_order.size() / 4);
    deltaInInt8Value.reserve(pe_reads_order.size() / 8);

    std::vector<bool> visit(pe_reads_order.size(), false);
    std::vector<uint32_t> pe_reads_order_index(pe_reads_order.size());
    for (uint32_t i = 0; i < pe_reads_order.size(); ++i) {
        pe_reads_order_index[pe_reads_order[i]] = i;
    }
    int64_t ref_prev_offset = 0;
    int64_t prev_offset = 0;
    bool match = false;
    for (uint32_t order_index = 0; order_index < pe_reads_order.size(); ++order_index) {
        if (visit[order_index]) continue;
        uint32_t order = pe_reads_order[order_index];
        uint32_t pair_order = order % 2 ? (order - 1) : (order + 1);
        uint32_t pair_order_index = pe_reads_order_index[pair_order];
        visit[pair_order_index] = true;
        int64_t pair_index_offset = pair_order_index - order_index;
        offsetInUint8Flag.push_back((uint8_t)(pair_index_offset <= UINT8_MAX));
        if (pair_index_offset <= UINT8_MAX) {
            offsetInUint8Value.push_back((uint8_t) pair_index_offset);
            offsetPairFlag.push_back(order % 2);
            continue;
        }
        nonOffsetPairFlag.push_back(order % 2);
        const int64_t delta = pair_index_offset - ref_prev_offset;
        const bool isDeltaInInt8 = INT8_MIN <= delta && delta <= INT8_MAX;
        deltaInInt8Flag.push_back((uint8_t) isDeltaInInt8);
        if (isDeltaInInt8) {
            match = true;
            deltaInInt8Value.push_back((int8_t) delta);
            ref_prev_offset = pair_index_offset;
        } else {
            if (!match || ref_prev_offset != prev_offset) {
                ref_prev_offset = pair_index_offset;
            }
            offset.push_back(pair_index_offset);
            match = false;
        }
        prev_offset = pair_index_offset;
    }
    writeCompressed(out, (char*) offsetInUint8Flag.data(), offsetInUint8Flag.size() * sizeof(uint8_t),
                    PPMD7_CODER, 3, 3, COMPRESSION_ESTIMATION_UINT8_BITMAP);
    writeCompressed(out, (char*) offsetInUint8Value.data(), offsetInUint8Value.size() * sizeof(uint8_t),
                    PPMD7_CODER, 3, 2, 1);
    writeCompressed(out, (char*) deltaInInt8Flag.data(), deltaInInt8Flag.size() * sizeof(uint8_t),
                    PPMD7_CODER, 3, 3, COMPRESSION_ESTIMATION_UINT8_BITMAP);
    writeCompressed(out, (char*) offsetPairFlag.data(), offsetPairFlag.size() * sizeof(uint8_t),
                    PPMD7_CODER, 3, 2, COMPRESSION_ESTIMATION_UINT8_BITMAP);
    writeCompressed(out, (char*) nonOffsetPairFlag.data(), nonOffsetPairFlag.size() * sizeof(uint8_t),
                    PPMD7_CODER, 3, 2, COMPRESSION_ESTIMATION_UINT8_BITMAP);

    std::ofstream delta_file(working_path / "delta.bin");
    PgSAHelpers::writeArray(delta_file, deltaInInt8Value.data(), deltaInInt8Value.size() * sizeof(int8_t));
    lzma2::lzma2_compress((working_path / "delta.bin").c_str(), (working_path / "delta.lzma2").c_str(), level, t_num);
    read_file_and_output(working_path / "delta.lzma2");

    std::ofstream offset_file(working_path / "offset.bin");
    PgSAHelpers::writeArray(offset_file, offset.data(), offset.size() * sizeof(uint32_t));
    lzma2::lzma2_compress((working_path / "offset.bin").c_str(), (working_path / "offset.lzma2").c_str(), level, t_num);
    read_file_and_output(working_path / "offset.lzma2");
}

static inline void bsc_compress(const char *input, const char *output) {
    fs::path input_path(input);
    fs::path output_path(output);
    std::string bsc_command = std::string("./bsc") + " e "
                              + fs::absolute(input_path).c_str() + " "
                              + fs::absolute(output_path).c_str() + " -pm0 -e2 > /dev/null";
    int status = std::system(bsc_command.c_str());
    if (status != 0) throw std::runtime_error("Error occurred during bsc compress.");
}

template <class T>
struct PinnedMemoryAllocator {
    typedef T value_type;

    PinnedMemoryAllocator () = default;
    template <class U> constexpr explicit PinnedMemoryAllocator (const PinnedMemoryAllocator <U>&) noexcept {}

    [[nodiscard]] T* allocate(std::size_t n) {
#if defined(DEBUG) || defined(_DEBUG)
        if (n > std::numeric_limits<std::size_t>::max() / sizeof(T)) throw std::bad_array_new_length();
#endif
        T * p;
        cudaError_t ret = cudaMallocHost((void**)&p, n * sizeof(T));
#if defined(DEBUG) || defined(_DEBUG)
        if (ret != cudaSuccess) {
            fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(ret));
            throw std::bad_alloc();
        }
#endif
        return p;
    }

    void deallocate(T* p, std::size_t n) noexcept {
        cudaFreeHost(p);
    }

private:

};

template <class T, class U>
bool operator==(const PinnedMemoryAllocator <T>&, const PinnedMemoryAllocator <U>&) { return true; }
template <class T, class U>
bool operator!=(const PinnedMemoryAllocator <T>&, const PinnedMemoryAllocator <U>&) { return false; }


#endif //CURC_UTIL_CUH
