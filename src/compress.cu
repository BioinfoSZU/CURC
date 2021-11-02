#include <immintrin.h>
#include <thrust/device_vector.h>
#include <thrust/binary_search.h>
#include <thrust/system/cuda/execution_policy.h>
#include <thrust/system/omp/execution_policy.h>
#include <mio/mio.hpp>
#include "compress.cuh"
#include "constant.hpp"
#include "util.cuh"

std::mutex gpu_lock[100]; // device count require <= 100
#define LOCK_START                                               \
finish = false;                                                  \
while(!finish) {                                                 \
    for (device_id = 0; device_id < device_count; ++device_id) { \
        if (gpu_lock[device_id].try_lock()) {                    \

#define LOCK_END                                                 \
            gpu_lock[device_id].unlock();                        \
            finish = true;                                       \
            break;                                               \
        }                                                        \
    }                                                            \
}                                                                \

struct AlignmentRecord {
    uint32_t read_id;
    uint64_t is_contain_N : 1,
             strand_id    : 1,
             pos          : 62;

    AlignmentRecord() : read_id(0), is_contain_N(0), strand_id(0), pos(0) {}
    AlignmentRecord(uint32_t read_id, uint64_t is_contain_N, uint64_t strand_id, uint64_t pos)
       : read_id(read_id), is_contain_N(is_contain_N), strand_id(strand_id), pos(pos) {}
};

// template<size_t read_unit_size>
struct RefRecord {
    std::vector<AlignmentRecord> record;
    std::string ref_string;
    uint64_t * binary_ref_string = nullptr;
    size_t binary_ref_size = 0;

    void check_binary_ref(const uint64_t * reads_db, size_t reads_count, size_t read_len, size_t read_unit_size) {
        std::vector<bool> had_access(reads_count, false);
        for (size_t i = 0; i < record.size(); ++i) {
            uint32_t read_id = record[i].read_id;
            if (had_access[read_id]) {
                printf("ref check fail because had access\n");
                std::exit(-1);
            }
            had_access[read_id] = true;
            for (size_t w = 0; w < read_len; w += 32) {
                auto base_number = std::min(32ul, read_len - w);
                auto w_a = extract_word_cpu(reads_db, read_id * read_unit_size, w, base_number);
                auto w_b = extract_word_unaligned_cpu(binary_ref_string, record[i].pos + w, base_number);
                if (w_a != w_b) {
                    printf("ref check fail because error match\n");
                    std::exit(-1);
                }
            }
        }
        std::cout << "ref check success \n";
    }

    void assembly (
            const uint64_t * reads_db,
            const uint32_t * ref_reads_id,
            uint32_t ref_reads_count,
            const uint32_t * next,
            const uint32_t * prev,
            const uint16_t * offset,
            size_t reads_count,
            size_t read_len,
            size_t read_unit_size) {
        std::vector<bool> visit(reads_count, false);
        for (uint32_t r = 0; r < ref_reads_count; ++r) {
            auto read_id = ref_reads_id[r];
            if (visit[read_id]) continue;
            auto i = read_id;
            while (prev[i] != INVALID && prev[i] != read_id) {
                i = prev[i];
            }
            if (prev[i] == read_id) {
                uint16_t max_off = offset[i];
                uint32_t max_off_read_id = i;
                while (i != read_id) {
                    visit[i] = true;
                    i = next[i];
                    if (offset[i] > max_off) {
                        max_off = offset[i];
                        max_off_read_id = i;
                    }
                }
                visit[read_id] = true;
                uint32_t start = next[max_off_read_id];
                while (start != max_off_read_id) {
                    record.emplace_back(start, 0, 0, ref_string.size());
                    for (size_t j = 0; j < offset[start]; ++j) {
                        ref_string.push_back(extract_base_cpu(reads_db, start * read_unit_size, j));
                    }
                    start = next[start];
                }
                record.emplace_back(max_off_read_id, 0, 0, ref_string.size());
                for (size_t j = 0; j < read_len; ++j) {
                    ref_string.push_back(extract_base_cpu(reads_db, max_off_read_id * read_unit_size, j));
                }
            } else { // prev[i] == INVALID
                while (next[i] != INVALID) {
                    visit[i] = true;
                    record.emplace_back(i, 0, 0, ref_string.size());
                    for (size_t j = 0; j < offset[i]; ++j) {
                        ref_string.push_back(extract_base_cpu(reads_db, i * read_unit_size, j));
                    }
                    i = next[i];
                }
                visit[i] = true;
                record.emplace_back(i, 0, 0, ref_string.size());
                for (size_t j = 0; j < read_len; ++j) {
                    ref_string.push_back(extract_base_cpu(reads_db, i * read_unit_size, j));
                }
            }
        }
    }

    void assembly (
            const char* N_reads_db,
            const uint32_t * N_ref_reads_id,
            uint32_t N_ref_reads_count,
            const uint32_t * next,
            const uint32_t * prev,
            const uint16_t * offset,
            size_t N_reads_count,
            size_t read_len) {
        std::vector<bool> visit(N_reads_count, false);
        for (uint32_t r = 0; r < N_ref_reads_count; ++r) {
            auto read_id = N_ref_reads_id[r];
            if (visit[read_id]) continue;
            auto i = read_id;
            while (prev[i] != INVALID && prev[i] != read_id) {
                i = prev[i];
            }
            if (prev[i] == read_id) {
                uint16_t max_off = offset[i];
                uint32_t max_off_read_id = i;
                while (i != read_id) {
                    visit[i] = true;
                    i = next[i];
                    if (offset[i] > max_off) {
                        max_off = offset[i];
                        max_off_read_id = i;
                    }
                }
                visit[read_id] = true;
                uint32_t start = next[max_off_read_id];
                while (start != max_off_read_id) {
                    record.emplace_back(start, 1, 0, ref_string.size());
                    ref_string += std::string(N_reads_db + start * read_len, offset[start]);
                    start = next[start];
                }
                record.emplace_back(max_off_read_id, 1, 0, ref_string.size());
                ref_string += std::string(N_reads_db + max_off_read_id * read_len, read_len);
            } else { // prev[i] == INVALID
                while (next[i] != INVALID) {
                    visit[i] = true;
                    record.emplace_back(i, 1, 0, ref_string.size());
                    ref_string += std::string(N_reads_db + i * read_len, offset[i]);
                    i = next[i];
                }
                visit[i] = true;
                record.emplace_back(i, 1, 0, ref_string.size());
                ref_string += std::string(N_reads_db + i * read_len, read_len);
            }

        }
    }

    void allocate_binary_ref() {
        binary_ref_size = ceil<32>(ref_string.size()) + 1;
        cudaHostAlloc((void**) &binary_ref_string, binary_ref_size * sizeof(uint64_t), cudaHostAllocPortable);
        std::uninitialized_fill(binary_ref_string, binary_ref_string + binary_ref_size, 0);
    }

    void reset_binary_ref() {
        std::uninitialized_fill(binary_ref_string, binary_ref_string + binary_ref_size, 0);
    }

    void compute_binary_ref() {
        if (binary_ref_string == nullptr) {
            printf("ERROR: binary ref is null\n");
            std::exit(-1);
        }
#pragma omp parallel for
        for (size_t i = 0; i < ref_string.size(); i += 32) {
            uint64_t idx1 = i / 32, c;
            uint8_t idx2 = 0;

            size_t start = i, end = std::min(i + 32, ref_string.size());
            for (size_t j = start; j < end; j++) {
                c = base_encode_table_cpu[ref_string[j]] & 0x03U;
                binary_ref_string[idx1] |= (c << (62U - idx2));
                idx2 += 2;
            }
        }
    }

    void erase_binary_ref() {
        if (binary_ref_string != nullptr) {
            cudaFreeHost(binary_ref_string);
            binary_ref_size = 0;
            binary_ref_string = nullptr;
        }
    }

    void print(const std::string& ref_name, std::ostream& out = std::cout) const {
        out << ref_name << " {\n";
        out << "\treference reads count is " << record.size() << "\n";
        out << "\treference length is " << ref_string.size() << "\n";
        out << "}\n\n";
    }

    ~RefRecord() {
        erase_binary_ref();
    }
};

// template<size_t read_unit_size>
struct ReadMapper {
private:
    int device_id;
    size_t ref_len;
    size_t unmatch_reads_count;
    size_t N_reads_count;
    const Param& param;

    uint64_t * reads_db;
    char * N_reads_db;
    uint64_t * ref;

    AlignmentRecord * matches;
    uint8_t * last_mismatch_count;

    uint32_t * key_ranges;
    uint32_t * pos_array;
    uint32_t hash_size_minus_one;

public:
    ReadMapper(int device_id, size_t ref_len, size_t unmatch_reads_count, size_t N_reads_count, const Param& param)
    : device_id(device_id), ref_len(ref_len), unmatch_reads_count(unmatch_reads_count), N_reads_count(N_reads_count), param(param) {}

    void init(const uint64_t * reads_db_host,  const char* N_reads_db_host, const uint32_t * unmatch_reads_id,
              const uint64_t * binary_ref, size_t binary_ref_size) {
        cudaSetDevice(device_id);
        size_t read_unit_size = param.read_unit_size;
        uint64_t * reads_db_buffer;
        cudaMallocHost((void**) &reads_db_buffer, unmatch_reads_count * read_unit_size * sizeof(uint64_t));
#pragma omp parallel for
        for (size_t i = 0; i < unmatch_reads_count; ++i) {
#pragma omp simd
            for (size_t j = 0; j < read_unit_size; ++j) {
                reads_db_buffer[i * read_unit_size + j] = reads_db_host[unmatch_reads_id[i] * read_unit_size + j];
            }
        }
        gpuErrorCheck(cudaMalloc((void**) &reads_db, unmatch_reads_count * read_unit_size * sizeof(uint64_t)));
        cudaMemcpy(reads_db, reads_db_buffer, unmatch_reads_count * read_unit_size * sizeof(uint64_t), cudaMemcpyHostToDevice);
        cudaFreeHost(reads_db_buffer);

        gpuErrorCheck(cudaMalloc((void**) &N_reads_db, N_reads_count * param.read_len));
        cudaMemcpy(N_reads_db, N_reads_db_host, N_reads_count * param.read_len, cudaMemcpyHostToDevice);

        gpuErrorCheck(cudaMalloc((void**) &ref, binary_ref_size * sizeof(uint64_t)));
        cudaMemcpy(ref, binary_ref, binary_ref_size * sizeof(uint64_t), cudaMemcpyHostToDevice);

        gpuErrorCheck(cudaMalloc((void**) &matches, (unmatch_reads_count + N_reads_count) * sizeof(AlignmentRecord)));
        gpuErrorCheck(cudaMalloc((void**) &last_mismatch_count, (unmatch_reads_count + N_reads_count) * sizeof(uint8_t)));
        thrust::uninitialized_fill(thrust::device,
                                   last_mismatch_count,
                                   last_mismatch_count + unmatch_reads_count + N_reads_count,
                                   param.max_mismatch_count + 1);
        thrust::for_each(thrust::device,
                         thrust::counting_iterator<uint32_t>(0),
                         thrust::counting_iterator<uint32_t>(unmatch_reads_count + N_reads_count),
                         [unmatch_reads_count = unmatch_reads_count, matches = matches] __device__(uint32_t i) {
                             if (i < unmatch_reads_count) {
                                 matches[i].read_id = i;
                                 matches[i].is_contain_N = 0;
                             } else {
                                 matches[i].read_id = i - unmatch_reads_count;
                                 matches[i].is_contain_N = 1;
                             }
                         });
    }

    template<bool is_reverse_complement>
    void index_and_mapping(size_t k, size_t step, size_t bucket_limit, size_t read_index_step, size_t target_mismatch_count) {
        cudaSetDevice(device_id);
        size_t index_count = (ref_len - k + 1 + step - 1) / step;
        if (index_count > (1ULL << 30ULL)) { // two-step build index and mapping
            size_t first_index_count = index_count / 2, second_index_count = index_count - first_index_count;
            size_t first_index_start_pos = 0, second_index_start_pos = first_index_count * step;
            build_index<is_reverse_complement>(k, step, bucket_limit, first_index_start_pos, first_index_count);
            mapping<is_reverse_complement>(k, step, read_index_step, first_index_start_pos, target_mismatch_count);
            cudaFree(key_ranges);
            cudaFree(pos_array);

            build_index<is_reverse_complement>(k, step, bucket_limit, second_index_start_pos, second_index_count);
            mapping<is_reverse_complement>(k, step, read_index_step, second_index_start_pos, target_mismatch_count);
            cudaFree(key_ranges);
            cudaFree(pos_array);
        } else {
            build_index<is_reverse_complement>(k, step, bucket_limit, 0, index_count);
            mapping<is_reverse_complement>(k, step, read_index_step, 0, target_mismatch_count);
            cudaFree(key_ranges);
            cudaFree(pos_array);
        }
    }

    template<bool is_reverse_complement>
    void build_index(size_t k, size_t step, size_t bucket_limit, size_t index_start_pos, size_t index_count) {
        cudaSetDevice(device_id);
        uint32_t hash_size = get_hash_size(index_count);
        hash_size_minus_one = hash_size - 1;
        uint8_t * counts;
        gpuErrorCheck(cudaMalloc((void**) &counts, hash_size + 2));
        thrust::uninitialized_fill(thrust::device, counts, counts + hash_size + 2, 0);
        thrust::for_each(thrust::device,
                         thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(index_count),
                         [ref = ref, hash_size_minus_one = hash_size_minus_one, step, k, bucket_limit, index_start_pos, counts]
                         __device__ (size_t i) {
                             uint64_t pos = index_start_pos + i * step;
                             uint64_t kmer = extract_word_unaligned_gpu(ref, pos, k);
                             if /*constexpr*/ (is_reverse_complement) {
                                 kmer = get_reverse_kmer(kmer, k);
                             }
                             uint32_t hash_pos = murmur_hash64(kmer) & hash_size_minus_one;
                             if (counts[hash_pos] < bucket_limit) ++counts[hash_pos];
                         });

        gpuErrorCheck(cudaMalloc((void**) &key_ranges, (hash_size + 2) * sizeof(uint32_t)));
        thrust::exclusive_scan(thrust::device, counts, counts + hash_size + 2, key_ranges, 0);

        uint32_t hash_count;
        cudaMemcpy(&hash_count, key_ranges + hash_size + 1, sizeof(uint32_t), cudaMemcpyDeviceToHost);
        gpuErrorCheck(cudaMalloc((void**) &pos_array, (hash_count + 2) * sizeof(uint32_t)));

        thrust::for_each(thrust::device,
                         thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(index_count),
                         [ref=ref, key_ranges=key_ranges, pos_array=pos_array, hash_size_minus_one = hash_size_minus_one, step, k, index_start_pos, counts]
                         __device__ (size_t i) {
                             uint64_t pos = index_start_pos + i * step;
                             uint64_t kmer = extract_word_unaligned_gpu(ref, pos, k);
                             if /*constexpr*/ (is_reverse_complement) {
                                 kmer = get_reverse_kmer(kmer, k);
                             }
                             uint32_t hash_pos = murmur_hash64(kmer) & hash_size_minus_one;
                             if (counts[hash_pos]) {
                                 pos_array[key_ranges[hash_pos] + (--counts[hash_pos])] = i;
                             }
                         });

        cudaFree(counts);
    }

    template<bool is_reverse_complement>
    void mapping (size_t k, size_t step, size_t read_index_step, size_t index_start_pos, size_t target_mismatch_count) {
        cudaSetDevice(device_id);
        thrust::for_each(thrust::device,
                         thrust::counting_iterator<uint32_t>(0),
                         thrust::counting_iterator<uint32_t>(unmatch_reads_count + N_reads_count),
                                 [last_mismatch_count = last_mismatch_count, matches = matches,
                                 ref = ref, reads_db = reads_db, N_reads_db = N_reads_db,
                                 key_ranges = key_ranges, pos_array = pos_array, hash_size_minus_one = hash_size_minus_one,
                                 read_len = param.read_len, read_unit_size = param.read_unit_size, ref_len = ref_len,
                                 k, step, read_index_step, index_start_pos, target_mismatch_count] __device__(uint32_t i) {
                             if (last_mismatch_count[i] <= target_mismatch_count) return;
                             uint64_t is_contain_N = matches[i].is_contain_N;
                             if (!is_contain_N) {
                                 uint32_t read_id = matches[i].read_id;
                                 for (size_t j = 0; j < read_len - k + 1; j += read_index_step) {
                                     size_t read_position = j;
                                     uint64_t read_key = extract_word_gpu(reads_db, read_id * read_unit_size, read_position, k);
                                     auto hash_key_pos = murmur_hash64(read_key) & hash_size_minus_one;
                                     auto begin = key_ranges[hash_key_pos], end = key_ranges[hash_key_pos + 1];
                                     if (begin != end) {
                                         end = thrust::min(begin + 12, end);
                                         for (size_t p = begin; p < end; ++p) {
                                             uint64_t ref_position = index_start_pos + (uint64_t) pos_array[p] * step;
                                             uint64_t left_shift = read_position;
                                             if /*constexpr*/ (is_reverse_complement)
                                                 left_shift = read_len - k - read_position;
                                             if (ref_position < left_shift) continue;
                                             uint64_t match_pos = ref_position - left_shift;
                                             if (match_pos + read_len > ref_len) continue;

                                             if /*constexpr*/ (!is_reverse_complement) { // same direction
                                                 auto mismatch_count = forward_mismatch_count(ref, match_pos, reads_db, read_id * read_unit_size, read_len, last_mismatch_count[i] - 1);
                                                 if (mismatch_count < last_mismatch_count[i]) {
                                                     matches[i].strand_id = 0;
                                                     matches[i].pos = match_pos;
                                                     last_mismatch_count[i] = mismatch_count;
                                                     if (last_mismatch_count[i] <= target_mismatch_count)
                                                         return;
                                                 }
                                             } else {  // different direction
                                                 auto mismatch_count = reverse_complement_mismatch_count(ref, match_pos, reads_db, read_id * read_unit_size, read_len, last_mismatch_count[i] - 1);
                                                 if (mismatch_count < last_mismatch_count[i]) {
                                                     matches[i].strand_id = 1;
                                                     matches[i].pos = match_pos;
                                                     last_mismatch_count[i] = mismatch_count;
                                                     if (last_mismatch_count[i] <= target_mismatch_count)
                                                         return;
                                                 }
                                             }
                                         }
                                     }
                                 }
                             } else {
                                 uint32_t N_read_id = matches[i].read_id;
                                 for (size_t j = 0; j < read_len - k + 1; j += read_index_step) {
                                     size_t read_position = j;
                                     uint64_t read_key = get_kmer_integer(N_reads_db + N_read_id * read_len + read_position, k);
                                     auto hash_key_pos = murmur_hash64(read_key) & hash_size_minus_one;
                                     auto begin = key_ranges[hash_key_pos], end = key_ranges[hash_key_pos + 1];
                                     if (begin != end) {
                                         end = thrust::min(begin + 12, end);
                                         for (size_t p = begin; p < end; ++p) {
                                             uint64_t ref_position = index_start_pos + (uint64_t) pos_array[p] * step;
                                             uint64_t left_shift = read_position;
                                             if /*constexpr*/ (is_reverse_complement)
                                                 left_shift = read_len - k - read_position;
                                             if (ref_position < left_shift) continue;
                                             uint64_t match_pos = ref_position - left_shift;
                                             if (match_pos + read_len > ref_len) continue;

                                             if /*constexpr*/ (!is_reverse_complement) {
                                                 size_t mismatch_count = 0;
                                                 for (size_t o = 0; o < read_len; ++o) {
                                                     if (extract_base_unaligned_gpu(ref, match_pos + o) !=
                                                         N_reads_db[N_read_id * read_len + o]) {
                                                         mismatch_count++;
                                                         if (mismatch_count >= last_mismatch_count[i]) break;
                                                     }
                                                 }
                                                 if (mismatch_count < last_mismatch_count[i]) {
                                                     matches[i].strand_id = 0;
                                                     matches[i].pos = match_pos;
                                                     last_mismatch_count[i] = mismatch_count;
                                                     if (last_mismatch_count[i] <= target_mismatch_count)
                                                         return;
                                                 }
                                             } else {
                                                 size_t mismatch_count = 0;
                                                 for (size_t o = 0; o < read_len; ++o) {
                                                     if (extract_base_unaligned_gpu(ref, match_pos + o) !=
                                                         reverse_dna_base(N_reads_db[N_read_id * read_len + read_len - 1 - o])) {
                                                         mismatch_count++;
                                                         if (mismatch_count >= last_mismatch_count[i]) break;
                                                     }
                                                 }
                                                 if (mismatch_count < last_mismatch_count[i]) {
                                                     matches[i].strand_id = 1;
                                                     matches[i].pos = match_pos;
                                                     last_mismatch_count[i] = mismatch_count;
                                                     if (last_mismatch_count[i] <= target_mismatch_count)
                                                         return;
                                                 }
                                             }

                                         }
                                     }
                                 }
                             }
                         });
    }

    /// \param unmapping_reads_id reuse unmatch_reads_id memory
    void get_result(uint32_t* unmapping_reads_id, uint32_t & unmapping_reads_count,
                    uint32_t* & unmapping_N_reads_id, uint32_t & unmapping_N_reads_count,
                    std::vector<AlignmentRecord>& map_record, std::vector<uint8_t>& map_mis_cnt) {
        cudaSetDevice(device_id);
        cudaFree(reads_db);
        cudaFree(N_reads_db);
        cudaFree(ref);
        uint32_t * unmapping_reads_id_device;
        uint32_t * unmapping_N_reads_id_device;
        gpuErrorCheck(cudaMalloc((void**) &unmapping_reads_id_device, unmatch_reads_count * sizeof(uint32_t)));
        gpuErrorCheck(cudaMalloc((void**) &unmapping_N_reads_id_device, N_reads_count * sizeof(uint32_t)));
        cudaMemcpy(unmapping_reads_id_device, unmapping_reads_id, unmatch_reads_count * sizeof(uint32_t), cudaMemcpyHostToDevice);
        thrust::sequence(thrust::device, unmapping_N_reads_id_device, unmapping_N_reads_id_device + N_reads_count, 0);
        thrust::for_each(thrust::device, thrust::counting_iterator<uint32_t>(0), thrust::counting_iterator<uint32_t>(unmatch_reads_count),
                [matches = matches, unmapping_reads_id_device]__device__(uint32_t i) {
            matches[i].read_id = unmapping_reads_id_device[i];
        });
        uint32_t * indices;
        gpuErrorCheck(cudaMalloc((void**)&indices, (unmatch_reads_count + N_reads_count) * sizeof(uint32_t)));
        thrust::sequence(thrust::device, indices, indices + unmatch_reads_count + N_reads_count, 0);
        unmapping_reads_count =
                thrust::remove_if(thrust::device,
                                  thrust::make_zip_iterator(thrust::make_tuple(unmapping_reads_id_device, indices)),
                                  thrust::make_zip_iterator(thrust::make_tuple(unmapping_reads_id_device + unmatch_reads_count, indices + unmatch_reads_count)),
                                  [last_mismatch_count = last_mismatch_count, max_mismatch_count = param.max_mismatch_count] __device__ (const thrust::tuple<uint32_t, uint32_t>& item) {
                    return last_mismatch_count[thrust::get<1>(item)] <= max_mismatch_count;
                }) - thrust::make_zip_iterator(thrust::make_tuple(unmapping_reads_id_device, indices));
        unmapping_N_reads_count =
                thrust::remove_if(thrust::device,
                                  thrust::make_zip_iterator(thrust::make_tuple(unmapping_N_reads_id_device, indices + unmatch_reads_count)),
                                  thrust::make_zip_iterator(thrust::make_tuple(unmapping_N_reads_id_device + N_reads_count, indices + unmatch_reads_count + N_reads_count)),
                                  [last_mismatch_count = last_mismatch_count, max_mismatch_count = param.max_mismatch_count] __device__ (const thrust::tuple<uint32_t,uint32_t>& item) {
                    return last_mismatch_count[thrust::get<1>(item)] <= max_mismatch_count;
                }) - thrust::make_zip_iterator(thrust::make_tuple(unmapping_N_reads_id_device, indices + unmatch_reads_count));

        if (unmapping_N_reads_count > 0) {
            cudaMallocHost((void **) &unmapping_N_reads_id, unmapping_N_reads_count * sizeof(uint32_t));
            cudaMemcpy(unmapping_N_reads_id, unmapping_N_reads_id_device, unmapping_N_reads_count * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        }
        if (unmapping_reads_count > 0) {
            cudaMemcpy(unmapping_reads_id, unmapping_reads_id_device, unmapping_reads_count * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        }

        cudaFree(indices);
        cudaFree(unmapping_reads_id_device);
        cudaFree(unmapping_N_reads_id_device);

        size_t map_record_size =
                thrust::remove_if(thrust::device,
                                  thrust::make_zip_iterator(thrust::make_tuple(matches, last_mismatch_count)),
                                  thrust::make_zip_iterator(thrust::make_tuple(matches + unmatch_reads_count + N_reads_count, last_mismatch_count + unmatch_reads_count + N_reads_count)),
                                  [max_mismatch_count = param.max_mismatch_count] __device__ (const thrust::tuple<AlignmentRecord, uint8_t>& item) {
                    return thrust::get<1>(item) > max_mismatch_count;
                }) - thrust::make_zip_iterator(thrust::make_tuple(matches, last_mismatch_count));
        map_record.resize(map_record_size);
        map_mis_cnt.resize(map_record_size);
        cudaMemcpy(map_record.data(), matches, map_record_size * sizeof(AlignmentRecord), cudaMemcpyDeviceToHost);
        cudaMemcpy(map_mis_cnt.data(), last_mismatch_count, map_record_size * sizeof(uint8_t), cudaMemcpyDeviceToHost);
        cudaFree(matches);
        cudaFree(last_mismatch_count);
    }
};

template <size_t K>
__device__ __forceinline__ std::uint32_t maRushPrime1HashBinary(const uint64_t * binary_ref, size_t base_start) {
    char kmer_str[K];
#pragma unroll
    for (size_t i = 0; i < K; ++i) {
        kmer_str[i] = extract_base_unaligned_gpu(binary_ref, base_start + i);
    }
    std::uint64_t hash = K;
    char * kmer = kmer_str;
    for (std::uint32_t j = 0; j < K/4; ) {
        std::uint32_t k;
        memcpy(&k, kmer, 4);
        k += j++;
        hash ^= k;
        hash *= 171717;
        kmer += 4;
    }
    return (std::uint32_t)(hash);
}

template <size_t K>
__device__ __forceinline__ std::uint32_t maRushPrime1HashChar(const char *str) {
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

struct RefMatcher {
    struct MatchResult {
        uint64_t posSrc;
        uint64_t length;
        uint64_t posDest;

        __device__ MatchResult(): posSrc(0), length(0), posDest(0) {}

        __device__ MatchResult(uint64_t posSrc, uint64_t length, uint64_t posDest) : posSrc(posSrc), length(length), posDest(posDest) {}

        __device__ bool operator==(const MatchResult &rhs) const {
            return posSrc == rhs.posSrc &&
                   length == rhs.length &&
                   posDest == rhs.posDest;
        }

        __device__ bool operator!=(const MatchResult &rhs) const {
            return !(rhs == *this);
        }

        __device__ bool operator<(const MatchResult &rhs) const {
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

    static void reverse_complement(std::string& seq) {
        size_t L = seq.size();
#pragma omp parallel for
        for (size_t i = 0; i < L / 2; ++i) {
            char a = PgSAHelpers::reverseComplement(seq[i]);
            char b = PgSAHelpers::reverseComplement(seq[L - i - 1]);
            seq[i] = b;
            seq[L - i - 1] = a;
        }
        if (L % 2) seq[L / 2] = PgSAHelpers::reverseComplement(seq[L / 2]);
    }

//    static string getTotalMatchStat(uint64_t totalMatchLength, uint64_t destPgLength) {
//        return PgSAHelpers::toString(totalMatchLength) + " (" + PgSAHelpers::toString((totalMatchLength * 100.0) / destPgLength, 1) + "%)";
//    }

    static constexpr size_t src_step = 5;
    static constexpr size_t dest_step = 3;
    static constexpr size_t kmer_size = 36;
    static constexpr size_t target_match_len = 50;
    static constexpr size_t bucket_size_limit = 12;
    static constexpr double ref_match_result_ratio = 0.15;
    static constexpr double unmap_ref_match_result_ratio = 0.1;

private:
    size_t ref_len;
    int device_id;
    const Param& param;

    uint64_t * ref;
    uint32_t * key_ranges;
    uint32_t * pos_array;
    uint32_t hash_size_minus_one;

public:
    RefMatcher(size_t ref_len, int device_id, const Param& param): ref_len(ref_len), device_id(device_id), param(param) {}

    void init(const uint64_t * binary_ref, size_t binary_ref_size) {
        cudaSetDevice(device_id);
        size_t index_count = (ref_len - kmer_size + 1 + src_step - 1) / src_step;
        if (index_count > (1ULL << 30ULL)) { // ref_len require <= 5.3 Gbp, this is generally true under the block compression strategy
            printf("RefMatcher hash_size is limit to 4GB\n");
            std::exit(-1);
        }
        auto hash_size = get_hash_size(index_count);
        hash_size_minus_one = hash_size - 1;
        gpuErrorCheck(cudaMalloc((void**) &ref, binary_ref_size * sizeof(uint64_t)));
        cudaMemcpy(ref, binary_ref, binary_ref_size * sizeof(uint64_t), cudaMemcpyHostToDevice);

        uint8_t * counts;
        gpuErrorCheck(cudaMalloc((void**) &counts, hash_size + 2));
        thrust::uninitialized_fill(thrust::device, counts, counts + hash_size + 2, 0);
        thrust::for_each(thrust::device,
                         thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(index_count),
                         [ref = ref, hash_size_minus_one = hash_size_minus_one, counts]
                                 __device__ (size_t i) {
                             uint64_t pos = i * src_step;
                             uint32_t hash_pos = maRushPrime1HashBinary<kmer_size>(ref, pos) & hash_size_minus_one;
                             if (counts[hash_pos] < bucket_size_limit) ++counts[hash_pos];
                         });

        gpuErrorCheck(cudaMalloc((void**) &key_ranges, (hash_size + 2) * sizeof(uint32_t)));
        thrust::exclusive_scan(thrust::device, counts, counts + hash_size + 2, key_ranges, 0);

        uint32_t hash_count;
        cudaMemcpy(&hash_count, key_ranges + hash_size + 1, sizeof(uint32_t), cudaMemcpyDeviceToHost);
        gpuErrorCheck(cudaMalloc((void**) &pos_array, (hash_count + 2) * sizeof(uint32_t)));

        thrust::for_each(thrust::device,
                         thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(index_count),
                         [ref=ref, key_ranges=key_ranges, pos_array=pos_array, hash_size_minus_one = hash_size_minus_one, counts]
                                 __device__ (size_t i) {
                             uint64_t pos = i * src_step;
                             uint32_t hash_pos = maRushPrime1HashBinary<kmer_size>(ref, pos) & hash_size_minus_one;
                             if (counts[hash_pos]) {
                                 pos_array[key_ranges[hash_pos] + (--counts[hash_pos])] = i;
                             }
                         });

        cudaFree(counts);
    }

    void match_binary_dest(const uint64_t * dest, size_t dest_len,
                           MatchResult * matches, unsigned long long * matches_size_atomic,
                           size_t dest_index_count, bool dest_is_ref) {
        cudaSetDevice(device_id);
        constexpr size_t dest_segement_step = 4;
        size_t dest_segement_count = (dest_index_count + dest_segement_step - 1) / dest_segement_step;
        thrust::for_each(thrust::device,
                         thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(dest_segement_count /*dest_index_count*/),
                [dest, matches, matches_size_atomic, dest_len, dest_is_ref,
                 ref=ref, ref_len=ref_len, key_ranges=key_ranges, pos_array=pos_array, hash_size_minus_one=hash_size_minus_one]
                __device__(size_t i) {
                    for (size_t s = 0; s < dest_segement_step; ++s) {
                        uint64_t dest_pos = (i * dest_segement_step + s) * dest_step;
                        uint32_t hash_pos = maRushPrime1HashBinary<kmer_size>(dest, dest_pos) & hash_size_minus_one;
                        uint64_t begin = key_ranges[hash_pos], end = key_ranges[hash_pos + 1];
                        if (begin != end) {
                            for (size_t j = begin; j < end; ++j) {
                                uint64_t src_pos = pos_array[j] * src_step;
                                if (dest_is_ref && (dest_len - src_pos < dest_pos)) continue;

                                uint64_t dest_curr = dest_pos;
                                uint64_t src_curr = src_pos;
                                uint64_t p1 = dest_curr + kmer_size - 1;
                                uint64_t p2 = src_curr + kmer_size - 1;
                                while (++p1 != dest_len &&
                                       ++p2 != ref_len &&
                                       extract_base_unaligned_gpu(dest, p1) == extract_base_unaligned_gpu(ref, p2));

                                uint64_t right = p1;
                                p1 = dest_curr;
                                p2 = src_curr;
                                while (p1 != 0 &&
                                       p2 != 0 &&
                                       extract_base_unaligned_gpu(dest, --p1) == extract_base_unaligned_gpu(ref, --p2));

                                if (right - p1 > target_match_len) {
                                    if (extract_word_unaligned_gpu(dest, dest_curr, 32) !=
                                        extract_word_unaligned_gpu(ref, src_curr, 32))
                                        continue;
                                    dest_curr += 32;
                                    src_curr += 32;
                                    int compare_len = kmer_size - 32;
                                    while (compare_len && extract_base_unaligned_gpu(dest, dest_curr++) ==
                                                          extract_base_unaligned_gpu(ref, src_curr++)) {
                                        compare_len--;
                                    }
                                    if (compare_len == 0) {
                                        unsigned long long matches_idx = atomicAdd(matches_size_atomic, 1);
                                        /// assert(matches_idx < match_result_capacity)
                                        /// ignore the range check, thrust::system::error may be throwed in some extreme cases
                                        matches[matches_idx] = MatchResult(p2 + 1, right - p1 - 1, p1 + 1);
                                        goto finish;
                                        // break;
                                    }
                                }
                            }
                        }
                    }
                    finish:;
        });
    }

    void match_char_dest(const char * dest, size_t dest_len, MatchResult * matches, unsigned long long * matches_size_atomic, size_t dest_index_count) {
        cudaSetDevice(device_id);
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(dest_index_count),
                         [dest, matches, matches_size_atomic, dest_len,
                          ref=ref, ref_len=ref_len, key_ranges=key_ranges, pos_array=pos_array, hash_size_minus_one=hash_size_minus_one]
                          __device__(size_t i) {
                             uint64_t dest_pos = i * dest_step;
                             uint32_t hash_pos = maRushPrime1HashChar<kmer_size>(dest + dest_pos) & hash_size_minus_one;
                             uint64_t begin = key_ranges[hash_pos], end = key_ranges[hash_pos + 1];
                             if (begin != end) {
                                 for (size_t j = begin; j < end; ++j) {
                                     uint64_t src_pos = pos_array[j] * src_step;

                                     uint64_t dest_curr = dest_pos;
                                     uint64_t src_curr = src_pos;
                                     uint64_t p1 = dest_curr + kmer_size - 1;
                                     uint64_t p2 = src_curr + kmer_size - 1;
                                     while(++p1 != dest_len &&
                                           ++p2 != ref_len &&
                                           dest[p1] == extract_base_unaligned_gpu(ref, p2)) ;

                                     uint64_t right = p1;
                                     p1 = dest_curr;
                                     p2 = src_curr;
                                     while (p1 != 0 &&
                                            p2 != 0 &&
                                            dest[--p1] == extract_base_unaligned_gpu(ref, --p2));

                                     if (right - p1 > target_match_len) {
                                         int compare_len = kmer_size;
                                         while(compare_len && dest[dest_curr++] == extract_base_unaligned_gpu(ref, src_curr++)) {
                                             compare_len--;
                                         }
                                         if (compare_len == 0) {
                                             unsigned long long matches_idx = atomicAdd(matches_size_atomic, 1);
                                             /// assert(matches_idx < match_result_capacity)
                                             /// ignore the range check, thrust::system::error may be throwed in some extreme cases
                                             matches[matches_idx] = MatchResult(p2 + 1, right - p1 - 1, p1 + 1);
                                             break;
                                         }
                                     }
                                 }
                             }
                         });
    }

    // template<size_t read_unit_size>
    void match(RefRecord * dest_ref, std::string& MapOff, std::string& MapLen, bool dest_contain_N, bool dest_is_ref) {
        cudaSetDevice(device_id);
        // auto dest_ref_size = dest_ref->ref_string.size();
        MatchResult * matches, * matches_host;
        unsigned long long matches_size;
        unsigned long long * matches_size_atomic;
        size_t dest_index_count = (dest_ref->ref_string.size() - kmer_size + 1 + dest_step - 1) / dest_step ;
        size_t match_result_capacity;
        if (dest_is_ref) {
            match_result_capacity = dest_index_count * ref_match_result_ratio;
        } else {
            match_result_capacity = dest_index_count * unmap_ref_match_result_ratio;
        }
        // printf("match result capacity : %zu \n", match_result_capacity);
        gpuErrorCheck(cudaMalloc((void**) &matches, match_result_capacity * sizeof(MatchResult)));
        cudaMalloc((void**) &matches_size_atomic, sizeof(unsigned long long));
        cudaMemset(matches_size_atomic, 0, sizeof(unsigned long long));

        reverse_complement(dest_ref->ref_string);
        if (!dest_contain_N) {
            if (dest_ref->binary_ref_string == nullptr) {
                dest_ref->allocate_binary_ref();
            } else {
                dest_ref->reset_binary_ref();
            }
            dest_ref->compute_binary_ref();
            uint64_t * dest;
            gpuErrorCheck(cudaMalloc((void**) &dest, dest_ref->binary_ref_size * sizeof(uint64_t)));
            cudaMemcpy(dest, dest_ref->binary_ref_string, dest_ref->binary_ref_size * sizeof(uint64_t), cudaMemcpyHostToDevice);
            match_binary_dest(dest, dest_ref->ref_string.size(), matches, matches_size_atomic, dest_index_count, dest_is_ref);
            cudaFree(dest);
            dest_ref->erase_binary_ref();
            if (dest_is_ref) {
                cudaFree(ref);
                cudaFree(key_ranges);
                cudaFree(pos_array);
            }
        } else {
            char * dest;
            gpuErrorCheck(cudaMalloc((void**) &dest, dest_ref->ref_string.size()));
            cudaMemcpy(dest, dest_ref->ref_string.data(), dest_ref->ref_string.size(), cudaMemcpyHostToDevice);
            match_char_dest(dest, dest_ref->ref_string.size(), matches, matches_size_atomic, dest_index_count);
            cudaFree(dest);
        }

        cudaMemcpy(&matches_size, matches_size_atomic, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
        // printf("matches result size : %llu\n", matches_size);
        if (matches_size == 0) {
            cudaFree(matches);
            cudaFree(matches_size_atomic);
            reverse_complement(dest_ref->ref_string);
            return;
        }
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(matches_size),
                         [matches, dest_len=dest_ref->ref_string.size()] __device__(size_t i) {
                             matches[i].posDest = dest_len - (matches[i].posDest + matches[i].length);
                         });
        if (dest_is_ref) {
            thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(matches_size),
            [matches] __device__ (size_t i) {
                if (matches[i].posSrc > matches[i].posDest) {
                    uint64_t tmp = matches[i].posSrc;
                    matches[i].posSrc = matches[i].posDest;
                    matches[i].posDest = tmp;
                }
                uint64_t endPosSrc = (matches[i].posSrc + matches[i].length);
                if (endPosSrc > matches[i].posDest) {
                    uint64_t margin = (endPosSrc - matches[i].posDest + 1) / 2;
                    matches[i].length -= margin;
                    matches[i].posDest += margin;
                }
            });
        }

        thrust::sort(thrust::device, matches, matches + matches_size);
        matches_size = thrust::unique(thrust::device, matches, matches + matches_size) - matches;
        // printf("Unique exact matches: %llu\n", matches_size);
        if (matches_size == 0) {
            cudaFree(matches);
            cudaFree(matches_size_atomic);
            reverse_complement(dest_ref->ref_string);
            return;
        }
        cudaMallocHost((void**) &matches_host, matches_size * sizeof(MatchResult));
        cudaMemcpy(matches_host, matches, matches_size * sizeof(MatchResult), cudaMemcpyDeviceToHost);
        cudaFree(matches);
        cudaFree(matches_size_atomic);
        reverse_complement(dest_ref->ref_string);

        std::ostringstream MapOffDest;
        std::ostringstream MapLenDest;
        PgSAHelpers::writeUIntByteFrugal(MapLenDest, target_match_len);

        char * dest_ptr = &dest_ref->ref_string[0];
        uint64_t pos = 0, npos = 0, totalDestOverlap = 0, totalMatched = 0;
        bool is_ref_len_std = ref_len <= UINT32_MAX;
        for (size_t i = 0; i < matches_size; ++i) {
            auto & match = matches_host[i];
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
            dest_ref->ref_string[npos++] = '%'; // mark
            if (is_ref_len_std) {
                PgSAHelpers::writeValue<uint32_t>(MapOffDest, match.posSrc);
            } else {
                PgSAHelpers::writeValue<uint64_t>(MapOffDest, match.posSrc);
            }
            PgSAHelpers::writeUIntByteFrugal(MapLenDest, match.length - target_match_len);
            pos = match.posDest + match.length;
        }
        uint64_t length = dest_ref->ref_string.size() - pos;
        std::memmove(dest_ptr + npos, dest_ptr + pos, length);
        npos += length;
        dest_ref->ref_string.resize(npos);

        cudaFreeHost(matches_host);
        MapOff = MapOffDest.str();
        MapOffDest.clear();
        MapLen = MapLenDest.str();
        MapLenDest.clear();
//        printf("Final size of Pg: %zu (remove: %s ; %zu chars in overlapped dest symbol)\n",
//               npos, getTotalMatchStat(totalMatched, dest_ref_size).c_str(), totalDestOverlap);
    }

//    void finish_match() {
//        cudaSetDevice(device_id);
//        cudaFree(ref);
//        cudaFree(key_ranges);
//        cudaFree(pos_array);
//    }
};

// template<size_t read_unit_size>
void full_match_gpu(uint64_t * reads_db, uint32_t * prefix_reads_id, uint32_t * suffix_reads_id, uint32_t cnt,
                    uint32_t * next, uint32_t * prev, uint16_t * offset, const Param& param, int device_id) {
    cudaSetDevice(device_id);

    constexpr size_t alphabet_size = 4;
    thrust::device_vector<char> alphabet(alphabet_size);
    alphabet[0] = 'A'; alphabet[1] = 'C'; alphabet[2] = 'G'; alphabet[3] = 'T';

    auto read_len = param.read_len;
    auto read_unit_size = param.read_unit_size;
    for (int off = 1; off < read_len; ++off) {
        thrust::device_vector<uint32_t> suffix_sort_part(alphabet_size);
        thrust::lower_bound(thrust::device,
                            suffix_reads_id, suffix_reads_id + cnt,
                            alphabet.begin(), alphabet.end(),
                            suffix_sort_part.begin(), [=] __device__ (const uint32_t & id, char value) {
                    return extract_base_gpu(reads_db, id * read_unit_size, off - 1) < value;
                });

        uint32_t A_part_begin = suffix_sort_part[0];
        uint32_t C_part_begin = suffix_sort_part[1];
        uint32_t G_part_begin = suffix_sort_part[2];
        uint32_t T_part_begin = suffix_sort_part[3];
        {
            auto compare = [=] __device__ (uint32_t a, uint32_t b) {
                for (size_t i = off; i < read_len; i += 32) {
                    auto base_number = thrust::min(32ul, read_len - i);
                    auto a_word = extract_word_gpu(reads_db, a * read_unit_size, i, base_number);
                    auto b_word = extract_word_gpu(reads_db, b * read_unit_size, i, base_number);
                    if (a_word < b_word) return true;
                    else if (a_word > b_word) return false;
                }
                return false;
            };

            uint32_t * A_C_merge_suffix_id;
            gpuErrorCheck(cudaMalloc((void**)&A_C_merge_suffix_id, (G_part_begin - A_part_begin) * sizeof(uint32_t)));
            thrust::merge(thrust::device,
                          suffix_reads_id + A_part_begin, suffix_reads_id + C_part_begin,
                          suffix_reads_id + C_part_begin, suffix_reads_id + G_part_begin,
                          A_C_merge_suffix_id, compare);

            uint32_t * G_T_merge_suffix_id;
            gpuErrorCheck(cudaMalloc((void**)&G_T_merge_suffix_id, (cnt - G_part_begin) * sizeof(uint32_t)));
            thrust::merge(thrust::device,
                          suffix_reads_id + G_part_begin, suffix_reads_id + T_part_begin,
                          suffix_reads_id + T_part_begin, suffix_reads_id + cnt,
                          G_T_merge_suffix_id, compare);

            thrust::merge(thrust::device,
                          A_C_merge_suffix_id, A_C_merge_suffix_id + (G_part_begin - A_part_begin),
                          G_T_merge_suffix_id, G_T_merge_suffix_id + (cnt - G_part_begin),
                          suffix_reads_id, compare);

            cudaFree(A_C_merge_suffix_id);
            cudaFree(G_T_merge_suffix_id);
        }

        size_t block_prefix_len = read_len - off >= 10 ? 10 : read_len - off;
        size_t block_num = std::pow(alphabet_size, block_prefix_len);

        thrust::device_vector<uint32_t> blocks_id(block_num);
        thrust::sequence(thrust::device, blocks_id.begin(), blocks_id.end(), 0);
        thrust::device_vector<uint32_t> prefix_id_range(block_num + 1), suffix_id_range(block_num + 1);

        {
            thrust::lower_bound(thrust::device,
                                prefix_reads_id, prefix_reads_id + cnt,
                                blocks_id.begin(), blocks_id.end(),
                                prefix_id_range.begin(),
                                [=] __device__ (const uint32_t& id, uint32_t b_id) {
                                    return extract_word_gpu(reads_db, id * read_unit_size, 0, block_prefix_len) < b_id;
                                });
            thrust::lower_bound(thrust::device,
                                suffix_reads_id, suffix_reads_id + cnt,
                                blocks_id.begin(), blocks_id.end(),
                                suffix_id_range.begin(),
                                [=] __device__ (const uint32_t& id, uint32_t b_id) {
                                    return extract_word_gpu(reads_db, id * read_unit_size, off, block_prefix_len) < b_id;
                                });
        }
        prefix_id_range[block_num] = cnt;
        suffix_id_range[block_num] = cnt;

        uint32_t * prefix_id_range_ptr = thrust::raw_pointer_cast(prefix_id_range.data());
        uint32_t * suffix_id_range_ptr = thrust::raw_pointer_cast(suffix_id_range.data());
        thrust::for_each(thrust::device, blocks_id.begin(), blocks_id.end(), [=] __device__ (const uint32_t b_id) {
            uint32_t prefix_id_begin = prefix_id_range_ptr[b_id], prefix_id_end = prefix_id_range_ptr[b_id + 1];
            uint32_t suffix_id_begin = suffix_id_range_ptr[b_id], suffix_id_end = suffix_id_range_ptr[b_id + 1];
            uint32_t prefix_id = prefix_id_begin, suffix_id = suffix_id_begin;

            while (prefix_id < prefix_id_end && suffix_id < suffix_id_end) {
                uint32_t prefix_read_id = prefix_reads_id[prefix_id];
                uint32_t suffix_read_id = suffix_reads_id[suffix_id];

                int cmp = 0;
                size_t len = read_len - off;
                for (size_t i = block_prefix_len; i < len; i += 32) {
                    auto base_number = thrust::min(32ul, len - i);
                    auto w_prefix = extract_word_gpu(reads_db, prefix_read_id * read_unit_size, i, base_number);
                    auto w_suffix = extract_word_gpu(reads_db, suffix_read_id * read_unit_size, i + off, base_number);
                    if (w_prefix > w_suffix) {
                        cmp = 1; break;
                    } else if (w_prefix < w_suffix) {
                        cmp = -1; break;
                    }
                }

                if (cmp > 0) {        // prefix_read > suffix_read
                    suffix_id++;
                } else if (cmp < 0) { // prefix_read < suffix_read
                    prefix_id++;
                } else {              // cmp == 0
                    if (prefix_read_id == suffix_read_id) {
                        auto current_prefix_id = prefix_id;
                        int tmp = 1;
                        while (prefix_id + 1 < prefix_id_end) {
                            prefix_id++;
                            prefix_read_id = prefix_reads_id[prefix_id];

                            cmp = 0;
                            for (size_t i = block_prefix_len; i < len; i += 32) {
                                auto base_number = thrust::min(32ul, len - i);
                                auto w_prefix = extract_word_gpu(reads_db, prefix_read_id * read_unit_size, i, base_number);
                                auto w_suffix = extract_word_gpu(reads_db, suffix_read_id * read_unit_size, i + off, base_number);
                                if (w_prefix > w_suffix) {
                                    cmp = 1;
                                    break;
                                }
                            }
                            tmp = cmp;

                            if (tmp != 0) break;
                            if (prefix_read_id != suffix_read_id) break;
                            tmp = 1;
                        }

                        if (tmp != 0) {
                            prefix_id = current_prefix_id;
                            suffix_id++;
                            continue;
                        } else {
                            prefix_read_id = prefix_reads_id[prefix_id];
                            while (prefix_id > current_prefix_id) {
                                prefix_reads_id[prefix_id] = prefix_reads_id[prefix_id - 1];
                                prefix_id--;
                            }
                            prefix_reads_id[prefix_id] = prefix_read_id;
                        }
                    }

                    next[suffix_read_id] = prefix_read_id;
                    prev[prefix_read_id] = suffix_read_id;
                    prefix_id++;
                    suffix_id++;
                }
            }
        });

        thrust::for_each(thrust::device, suffix_reads_id, suffix_reads_id + cnt,
                         [=]__device__(uint32_t read_id){
                             if (next[read_id] != INVALID) offset[read_id] = off;
                         });

        auto * new_end1 = thrust::remove_if(thrust::device, prefix_reads_id, prefix_reads_id + cnt,
                                            [=]__device__(uint32_t read_id){ return prev[read_id] != INVALID; });
        thrust::remove_if(thrust::device, suffix_reads_id, suffix_reads_id + cnt,
                          [=]__device__(uint32_t read_id) { return next[read_id] != INVALID; });

        cnt = new_end1 - prefix_reads_id;
    }
}

// template<size_t read_unit_size>
void full_match_cpu(const char* N_reads_db_host, uint32_t * unmapping_N_reads_id,
                    uint32_t unmapping_N_reads_count, size_t N_reads_count,
                    const Param& param, RefRecord * ref_ptr) {
    __gnu_parallel::sort(unmapping_N_reads_id, unmapping_N_reads_id + unmapping_N_reads_count, [&](auto a, auto b) {
        int ret = std::memcmp(N_reads_db_host + a * param.read_len, N_reads_db_host + b * param.read_len, param.read_len);
        return ret < 0;
    });
    std::vector<uint32_t> prev(N_reads_count, INVALID), next(N_reads_count, INVALID);
    std::vector<uint16_t> offset(N_reads_count, INVALOFF);

    __gnu_parallel::for_each(thrust::counting_iterator<uint32_t>(1), thrust::counting_iterator<uint32_t>(unmapping_N_reads_count), [&](auto i) {
        auto a = unmapping_N_reads_id[i - 1], b = unmapping_N_reads_id[i];
        if (std::memcmp(N_reads_db_host + a * param.read_len, N_reads_db_host + b * param.read_len, param.read_len) == 0 ) {
            prev[unmapping_N_reads_id[i]] = unmapping_N_reads_id[i - 1];
            next[unmapping_N_reads_id[i - 1]] = unmapping_N_reads_id[i];
            offset[unmapping_N_reads_id[i - 1]] = 0;
        }
    });

    std::vector<uint32_t> read_prefix_sort_id;
    std::vector<uint32_t> read_last_suffix_sort_id;
    std::vector<uint32_t> read_suffix_sort_id;
    read_prefix_sort_id.reserve(unmapping_N_reads_count);
    read_last_suffix_sort_id.reserve(unmapping_N_reads_count);
    for (size_t i = 0; i < unmapping_N_reads_count; ++i) {
        auto id = unmapping_N_reads_id[i];
        if (prev[id] == INVALID) {
            read_prefix_sort_id.push_back(id);
        }
        if (next[id] == INVALID) {
            read_last_suffix_sort_id.push_back(id);
        }
    }
    read_prefix_sort_id.shrink_to_fit();
    read_last_suffix_sort_id.shrink_to_fit();
    read_suffix_sort_id.resize(read_last_suffix_sort_id.size());

    constexpr size_t alphabet_size = 5;
    constexpr std::array<char, alphabet_size> alphabet = { 'A', 'C', 'G', 'N', 'T'};
    int off = 1;

    while (true) {
        std::array<std::pair<uint32_t *, uint32_t *>, alphabet_size> suffix_sort_part;
        std::array<size_t, 256> count = {0};
        for (size_t i = 0; i < read_last_suffix_sort_id.size(); ++i) {
            count[N_reads_db_host[read_last_suffix_sort_id[i] * param.read_len + off - 1]]++;
        }
        suffix_sort_part[0] = std::make_pair(read_last_suffix_sort_id.data(),
                                             read_last_suffix_sort_id.data() + count['A']);
        suffix_sort_part[1] = std::make_pair(read_last_suffix_sort_id.data() + count['A'],
                                             read_last_suffix_sort_id.data() + count['A'] + count['C']);
        suffix_sort_part[2] = std::make_pair(read_last_suffix_sort_id.data() + count['A'] + count['C'],
                                             read_last_suffix_sort_id.data() + count['A'] + count['C'] + count['G']);
        suffix_sort_part[3] = std::make_pair(read_last_suffix_sort_id.data() + count['A'] + count['C'] + count['G'],
                                             read_last_suffix_sort_id.data() + count['A'] + count['C'] + count['G'] + count['N']);
        suffix_sort_part[4] = std::make_pair(read_last_suffix_sort_id.data() + count['A'] + count['C'] + count['G'] + count['N'],
                                             read_last_suffix_sort_id.data() + read_last_suffix_sort_id.size());

        __gnu_parallel::stable_multiway_merge (suffix_sort_part.begin(), suffix_sort_part.end(),
                                               read_suffix_sort_id.data(), read_suffix_sort_id.size(),
                                               [&](uint32_t a, uint32_t b) {
                                                   int cmp = std::memcmp(N_reads_db_host + a * param.read_len + off ,
                                                                         N_reads_db_host + b * param.read_len + off,
                                                                         param.read_len - off);
                                                   return cmp < 0;
                                               });

        size_t block_prefix_len = param.read_len - off >= 3 ? 3 : param.read_len - off;
        size_t block_num = std::pow(alphabet_size, block_prefix_len);
        std::vector<uint32_t> prefix_id_range(block_num + 1), suffix_id_range(block_num + 1);

#pragma omp parallel for schedule(dynamic, 1)
        for (size_t block_prefix_integer = 0; block_prefix_integer < block_num; ++block_prefix_integer) {
            std::string block_prefix (block_prefix_len, 0);
            size_t temp = block_prefix_integer;
            for (size_t i = 0; i < block_prefix_len; ++i) {
                block_prefix[block_prefix_len - i - 1] = alphabet[temp % alphabet_size];
                temp /= alphabet_size;
            }
            prefix_id_range[block_prefix_integer] =
                    std::lower_bound(read_prefix_sort_id.begin(),
                                     read_prefix_sort_id.end(),
                                     block_prefix, [&](uint32_t id, const std::string& pattern) {
                                return std::memcmp(N_reads_db_host + id * param.read_len, pattern.data(), block_prefix_len) < 0;
                            }) - read_prefix_sort_id.begin();
            suffix_id_range[block_prefix_integer] =
                    std::lower_bound(read_suffix_sort_id.begin(),
                                     read_suffix_sort_id.end(),
                                     block_prefix, [&](uint32_t id, const std::string& pattern) {
                                return std::memcmp(N_reads_db_host + id * param.read_len + off, pattern.data(), block_prefix_len) < 0;
                            }) - read_suffix_sort_id.begin();
        }
        prefix_id_range[block_num] = read_prefix_sort_id.size();
        suffix_id_range[block_num] = read_suffix_sort_id.size();

#pragma omp parallel for schedule(dynamic, 1)
        for (size_t b_id = 0; b_id < block_num; ++b_id) {
            uint32_t prefix_id_begin = prefix_id_range[b_id], prefix_id_end = prefix_id_range[b_id + 1];
            uint32_t suffix_id_begin = suffix_id_range[b_id], suffix_id_end = suffix_id_range[b_id + 1];
            uint32_t prefix_id = prefix_id_begin, suffix_id = suffix_id_begin;
            while (prefix_id < prefix_id_end && suffix_id < suffix_id_end) {
                uint32_t prefix_read_id = read_prefix_sort_id[prefix_id];
                uint32_t suffix_read_id = read_suffix_sort_id[suffix_id];
                uint64_t prefix_start_index = prefix_read_id * param.read_len;
                uint64_t suffix_start_index = suffix_read_id * param.read_len;

                int cmp = std::memcmp(N_reads_db_host + prefix_start_index + block_prefix_len,
                                      N_reads_db_host + suffix_start_index + off + block_prefix_len,
                                      param.read_len - off - block_prefix_len);

                if (cmp > 0) {        // prefix_read > suffix_read
                    suffix_id++;
                } else if (cmp < 0) { // prefix_read < suffix_read
                    prefix_id++;
                } else {
                    if (prefix_read_id == suffix_read_id) {
                        auto current_prefix_id = prefix_id;
                        int tmp = 1;
                        while (prefix_id + 1 < prefix_id_end) {
                            prefix_id++;
                            prefix_read_id = read_prefix_sort_id[prefix_id];
                            prefix_start_index = prefix_read_id * param.read_len;
                            tmp = std::memcmp(N_reads_db_host + prefix_start_index + block_prefix_len,
                                              N_reads_db_host + suffix_start_index + off + block_prefix_len,
                                              param.read_len - off - block_prefix_len);
                            if (tmp != 0) break;
                            if (prefix_read_id != suffix_read_id) break;
                            tmp = 1;
                        }

                        if (tmp != 0) {
                            prefix_id = current_prefix_id;
                            suffix_id++;
                            continue;
                        } else {
                            prefix_read_id = read_prefix_sort_id[prefix_id];
                            while (prefix_id > current_prefix_id) {
                                read_prefix_sort_id[prefix_id] = read_prefix_sort_id[prefix_id - 1];
                                prefix_id--;
                            }
                            read_prefix_sort_id[prefix_id] = prefix_read_id;
                        }
                    }
                    next[suffix_read_id] = prefix_read_id;
                    prev[prefix_read_id] = suffix_read_id;
                    prefix_id++;
                    suffix_id++;
                }
            }
        }

        for (uint32_t read_id : read_suffix_sort_id) {
            if (next[read_id] != INVALID) {
                offset[read_id] = off;
            }
        }

        read_prefix_sort_id.erase(std::remove_if(read_prefix_sort_id.begin(), read_prefix_sort_id.end(),
                                                 [&](uint32_t read_id) { return prev[read_id] != INVALID; }),
                                  read_prefix_sort_id.end());

        read_suffix_sort_id.erase(std::remove_if(read_suffix_sort_id.begin(), read_suffix_sort_id.end(),
                                                 [&](uint32_t read_id) { return next[read_id] != INVALID; }),
                                  read_suffix_sort_id.end());

        read_last_suffix_sort_id = read_suffix_sort_id;

        off++;
        if (off >= param.read_len) break;
    }
    ref_ptr->assembly(N_reads_db_host, unmapping_N_reads_id, unmapping_N_reads_count, next.data(), prev.data(), offset.data(), N_reads_count, param.read_len);
}

// template<size_t read_unit_size>
void block_compress(const uint64_t * reads_db_host,
                    std::vector<uint32_t> reads_id,
                    std::vector<char> N_reads_db_host,
                    std::vector<uint32_t> N_reads_id,
                    const Param& param,
                    uint8_t block_id,
                    int device_count) {
    size_t read_unit_size = param.read_unit_size;
    fs::path working_path = fs::path(param.working_parent_path) / fs::path(std::to_string(block_id));
    if (!fs::exists(working_path)) {
        fs::create_directories(working_path);
    } else {
        fs::remove_all(working_path);
        fs::create_directories(working_path);
    }
    std::ofstream output_file(fs::path(param.working_dir) / (std::to_string(block_id) + block_archive_name_suffix));
    uint32_t reads_count = reads_id.size() + N_reads_id.size();
    output_file.write(reinterpret_cast<const char*>(&reads_count), sizeof(uint32_t));
    cudaHostRegister(reads_id.data(), reads_id.size() * sizeof(uint32_t), cudaHostRegisterPortable);
    cudaHostRegister(N_reads_db_host.data(), N_reads_db_host.size() * sizeof(char), cudaHostRegisterPortable);
    int device_id;
    bool finish;
    uint32_t * next_host, * prev_host;
    uint16_t * offset_host;
    uint32_t * ref_reads_id_host, * unmatch_reads_id_host;
    uint32_t ref_reads_count, unmatch_reads_count;

    LOCK_START
    {
        cudaSetDevice(device_id);

        uint64_t * reads_db_device;
        uint32_t * reads_id_device;
        uint32_t * next;
        uint32_t * prev;
        uint16_t * offset;
        gpuErrorCheck(cudaMalloc((void**) &reads_db_device, reads_count * read_unit_size * sizeof(uint64_t)));
        gpuErrorCheck(cudaMalloc((void**) &reads_id_device, reads_id.size() * sizeof(uint32_t)));
        gpuErrorCheck(cudaMalloc((void**) &next, reads_count * sizeof(uint32_t)));
        gpuErrorCheck(cudaMalloc((void**) &prev, reads_count * sizeof(uint32_t)));
        gpuErrorCheck(cudaMalloc((void**) &offset, reads_count * sizeof(uint16_t)));
        cudaMemcpy(reads_db_device, reads_db_host, reads_count * read_unit_size * sizeof(uint64_t), cudaMemcpyHostToDevice);
        cudaMemcpy(reads_id_device, reads_id.data(), reads_id.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);
        thrust::uninitialized_fill(thrust::device, next, next + reads_count, INVALID);
        thrust::uninitialized_fill(thrust::device, prev, prev + reads_count, INVALID);
        thrust::uninitialized_fill(thrust::device, offset, offset + reads_count, INVALOFF);

        thrust::sort(thrust::device, reads_id_device, reads_id_device + reads_id.size(),
                     [reads_db_device, read_unit_size] __device__ (uint32_t a, uint32_t b) {
            for (size_t i = 0; i < read_unit_size; ++i) {
                if (reads_db_device[a * read_unit_size + i] < reads_db_device[b * read_unit_size + i]) return true;
                if (reads_db_device[a * read_unit_size + i] > reads_db_device[b * read_unit_size + i]) return false;
            }
            return false;
        });

        thrust::for_each(thrust::device,
                         thrust::counting_iterator<uint32_t>(0),
                         thrust::counting_iterator<uint32_t>(reads_id.size() - 1),
                [reads_db_device, reads_id_device, next, prev, offset, read_unit_size] __device__ (uint32_t i) {
            for (size_t j = 0; j < read_unit_size; ++j) {
                if (reads_db_device[reads_id_device[i] * read_unit_size + j] !=
                    reads_db_device[reads_id_device[i + 1] * read_unit_size + j]) return;
            }
            next[reads_id_device[i]] = reads_id_device[i + 1];
            prev[reads_id_device[i + 1]] = reads_id_device[i];
            offset[reads_id_device[i]] = 0;
        });

        {
            auto read_len = param.read_len;
            auto kmer_size = param.kmer_size;
            auto kmer_index_pos = param.kmer_index_pos;

            uint32_t * value_array;
            uint32_t * suffix_reads_id;
            uint32_t duplicate_count = thrust::count_if(thrust::device, next, next + reads_count, []__device__(uint32_t id) { return id != INVALID; });
            uint32_t prefix_reads_count, suffix_reads_count;
            prefix_reads_count = suffix_reads_count = reads_id.size() - duplicate_count;
            gpuErrorCheck(cudaMalloc((void**) &value_array, prefix_reads_count * sizeof(uint32_t)));
            gpuErrorCheck(cudaMalloc((void**) &suffix_reads_id, suffix_reads_count * sizeof(uint32_t)));

            thrust::copy_if(thrust::device, reads_id_device, reads_id_device + reads_id.size(), value_array,
                            [prev]__device__(uint32_t id) { return prev[id] == INVALID; });
            thrust::copy_if(thrust::device, reads_id_device, reads_id_device + reads_id.size(), suffix_reads_id,
                            [next]__device__(uint32_t id) { return next[id] == INVALID; });

            uint32_t * key_range;
            auto hash_size = get_hash_size(prefix_reads_count);
            auto hash_size_minus_one = hash_size - 1;
            gpuErrorCheck(cudaMalloc((void**) &key_range, (hash_size + 2) * sizeof(uint32_t)));
            thrust::uninitialized_fill(thrust::device, key_range, key_range + hash_size + 2, 0);

            uint32_t * hash_key_array;
            gpuErrorCheck(cudaMalloc((void**) &hash_key_array, prefix_reads_count * sizeof(uint32_t)));
            thrust::for_each(thrust::device,
                             thrust::counting_iterator<uint32_t>(0),
                             thrust::counting_iterator<uint32_t>(prefix_reads_count),
                             [reads_db_device, value_array, hash_key_array, kmer_index_pos, kmer_size, hash_size_minus_one, read_unit_size]
                             __device__(uint32_t i) {
                                 uint64_t kmer = extract_word_gpu(reads_db_device, value_array[i] * read_unit_size, kmer_index_pos, kmer_size);
                                 hash_key_array[i] = murmur_hash64(kmer) & hash_size_minus_one;
                             });
            thrust::stable_sort_by_key(thrust::device, hash_key_array, hash_key_array + prefix_reads_count, value_array);

            uint32_t * aux;
            gpuErrorCheck(cudaMalloc((void**) &aux, (prefix_reads_count + 1) * sizeof(uint32_t)));
            thrust::sequence(thrust::device, aux, aux + prefix_reads_count + 1, 0);
            thrust::pair<uint32_t*, uint32_t*> new_end =
                    thrust::unique_by_key(thrust::device, hash_key_array, hash_key_array + prefix_reads_count, aux);
            cudaMemcpy(new_end.second, &prefix_reads_count, sizeof(uint32_t), cudaMemcpyHostToDevice);
            uint32_t bucket_count = new_end.first - hash_key_array;
            thrust::for_each(thrust::device, thrust::counting_iterator<uint32_t>(0), thrust::counting_iterator<uint32_t>(bucket_count),
                             [key_range, hash_key_array, aux]__device__(uint32_t i) {
                                 key_range[hash_key_array[i]] = aux[i + 1] - aux[i];
                             });
            thrust::exclusive_scan(thrust::device, key_range, key_range + hash_size + 2, key_range, 0);
            cudaFree(hash_key_array);
            cudaFree(aux);

            for (int off = 1; off < param.max_off; ++off) {
                thrust::for_each(thrust::device, suffix_reads_id, suffix_reads_id + suffix_reads_count,
                                 [=] __device__ (uint32_t suffix_read_id) {
                                     uint64_t kmer = extract_word_gpu(reads_db_device, suffix_read_id * read_unit_size, kmer_index_pos + off, kmer_size);
                                     uint32_t kmer_hash = murmur_hash64(kmer) & hash_size_minus_one;
                                     uint32_t hash_value_start = key_range[kmer_hash], hash_value_end = key_range[kmer_hash + 1];
                                     size_t compare_len = read_len - off;
                                     uint32_t count = hash_value_end - hash_value_start;
                                     auto comp = [compare_len, reads_db_device, off, read_unit_size] __device__ (uint32_t prefix_idx, uint32_t suffix_idx) {
                                         for (size_t j = 0; j < compare_len; j += 32) {
                                             auto base_number = thrust::min(32ul, compare_len - j);
                                             auto prefix = extract_word_gpu(reads_db_device, prefix_idx * read_unit_size, j, base_number);
                                             auto suffix = extract_word_gpu(reads_db_device, suffix_idx * read_unit_size, j + off, base_number);
                                             if (prefix < suffix) return true;
                                             else if (prefix > suffix) return false;
                                         }
                                         return false;
                                     };
                                     while (count > 0) {
                                         uint32_t it = hash_value_start;
                                         uint32_t step = count / 2;
                                         it += step;
                                         if (comp(value_array[it], suffix_read_id)) {
                                             hash_value_start = ++it;
                                             count -= step + 1;
                                         } else {
                                             count = step;
                                         }
                                     }
                                     for (auto hash_value_idx = hash_value_start; hash_value_idx < hash_value_end; ++hash_value_idx) {
                                         uint32_t prefix_read_id = value_array[hash_value_idx];
                                         if ((prefix_read_id == suffix_read_id) || (prev[prefix_read_id] != INVALID)) {
                                             continue;
                                         }
                                         bool is_match = true;
                                         for (size_t j = 0; j < compare_len; j += 32) {
                                             auto base_number = thrust::min(32ul, compare_len - j);
                                             auto prefix = extract_word_gpu(reads_db_device, prefix_read_id * read_unit_size, j, base_number);
                                             auto suffix = extract_word_gpu(reads_db_device, suffix_read_id * read_unit_size, j + off, base_number);
                                             if (prefix != suffix) {
                                                 is_match = false;
                                                 break;
                                             }
                                         }
                                         if (is_match) {
                                             uint32_t old = atomicCAS(prev + prefix_read_id, INVALID, suffix_read_id);
                                             if (old == INVALID) {
                                                 next[suffix_read_id] = prefix_read_id;
                                                 offset[suffix_read_id] = off;
                                                 break;
                                             }
                                         } else {
                                             break;
                                         }
                                     }
                });

                suffix_reads_count = thrust::remove_if(thrust::device, suffix_reads_id, suffix_reads_id + suffix_reads_count,
                                  [next] __device__ (uint32_t id) { return next[id] != INVALID; }) - suffix_reads_id;
            }

            cudaFree(value_array);
            cudaFree(suffix_reads_id);
            cudaFree(key_range);
        }

        {
            uint8_t * res;
            gpuErrorCheck(cudaMalloc((void**)&res, reads_count * sizeof(uint8_t)));
            thrust::uninitialized_fill(thrust::device, res, res + reads_count, true);

            thrust::for_each(thrust::device, thrust::counting_iterator<uint32_t>(0), thrust::counting_iterator<uint32_t>(reads_count),
                     [next, prev, offset, res] __device__ (uint32_t read_id) {
                         if (next[read_id] != INVALID && prev[read_id] != INVALID) return;
                         if (next[read_id] != INVALID && offset[read_id] == 0) return;
                         if (prev[read_id] != INVALID && offset[prev[read_id]] == 0) return;
                         res[read_id] = false;
                     });

            thrust::for_each(thrust::device, thrust::counting_iterator<uint32_t>(0), thrust::counting_iterator<uint32_t>(reads_count),
                             [next, prev, res] __device__ (uint32_t read_id) {
                                 if (!res[read_id]) {
                                     if (next[read_id] != INVALID) {
                                         prev[next[read_id]] = INVALID;
                                         next[read_id] = INVALID;
                                     }
                                     if (prev[read_id] != INVALID) {
                                         next[prev[read_id]] = INVALID;
                                         prev[read_id] = INVALID;
                                     }
                                 }
                             });

            ref_reads_count = thrust::stable_partition(thrust::device, reads_id_device, reads_id_device + reads_id.size(),
                                             [res] __device__(uint32_t id) { return res[id]; }) - reads_id_device;
            unmatch_reads_count = reads_id.size() - ref_reads_count;
            cudaMallocHost((void**) &ref_reads_id_host, ref_reads_count * sizeof(uint32_t));
            cudaHostAlloc((void**) &unmatch_reads_id_host, unmatch_reads_count * sizeof(uint32_t), cudaHostAllocPortable);
            cudaMemcpy(ref_reads_id_host, reads_id_device, ref_reads_count * sizeof(uint32_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(unmatch_reads_id_host, reads_id_device + ref_reads_count, unmatch_reads_count * sizeof(uint32_t), cudaMemcpyDeviceToHost);
            cudaFree(res);

            uint32_t * prefix_reads_id, * suffix_reads_id;
            uint32_t prefix_reads_count, suffix_reads_count;
            prefix_reads_count = suffix_reads_count = thrust::count_if(thrust::device, reads_id_device, reads_id_device + ref_reads_count,
                                     [next] __device__ (uint32_t id) { return next[id] == INVALID; });
            gpuErrorCheck(cudaMalloc((void**) &prefix_reads_id, prefix_reads_count * sizeof(uint32_t)));
            gpuErrorCheck(cudaMalloc((void**) &suffix_reads_id, suffix_reads_count * sizeof(uint32_t)));
            thrust::copy_if(thrust::device, reads_id_device, reads_id_device + ref_reads_count, prefix_reads_id,
                            [prev] __device__ (uint32_t id) { return prev[id] == INVALID; });
            thrust::copy_if(thrust::device, reads_id_device, reads_id_device + ref_reads_count, suffix_reads_id,
                            [next] __device__ (uint32_t id) { return next[id] == INVALID; });

            full_match_gpu(reads_db_device, prefix_reads_id, suffix_reads_id, suffix_reads_count,
                                           next, prev, offset, param, device_id);
            cudaFree(prefix_reads_id);
            cudaFree(suffix_reads_id);
        }

        cudaMallocHost((void**) &next_host, reads_count * sizeof(uint32_t));
        cudaMallocHost((void**) &prev_host, reads_count * sizeof(uint32_t));
        cudaMallocHost((void**) &offset_host, reads_count * sizeof(uint16_t));
        cudaMemcpy(next_host, next, reads_count * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(prev_host, prev, reads_count * sizeof(uint32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(offset_host, offset, reads_count * sizeof(uint16_t), cudaMemcpyDeviceToHost);

        cudaFree(reads_db_device);
        cudaFree(reads_id_device);
        cudaFree(next);
        cudaFree(prev);
        cudaFree(offset);
    }
    LOCK_END
    cudaHostUnregister(reads_id.data());
    reads_id.clear();
    reads_id.shrink_to_fit();

    auto ref = std::make_unique<RefRecord>();
    ref->assembly(reads_db_host, ref_reads_id_host, ref_reads_count, next_host, prev_host, offset_host, reads_count, param.read_len, param.read_unit_size);
    ref->allocate_binary_ref();
    ref->compute_binary_ref();
    cudaFreeHost(next_host);
    cudaFreeHost(prev_host);
    cudaFreeHost(offset_host);
    cudaFreeHost(ref_reads_id_host);

    uint32_t * unmapping_N_reads_id;
    uint32_t unmapping_reads_count, unmapping_N_reads_count;
    std::vector<AlignmentRecord> map_record;
    std::vector<uint8_t> map_mis_cnt;
    LOCK_START
    {
        cudaSetDevice(device_id);
        ReadMapper read_mapper(device_id, ref->ref_string.size(), unmatch_reads_count, N_reads_id.size(), param);
        read_mapper.init(reads_db_host, N_reads_db_host.data(), unmatch_reads_id_host, ref->binary_ref_string, ref->binary_ref_size);
        read_mapper.template index_and_mapping<false>(Param::k1, Param::ref_index_step1, Param::bucket_limit, Param::read_index_step1, Param::target_mismatch_count1);
        read_mapper.template index_and_mapping<true>(Param::k1, Param::ref_index_step1, Param::bucket_limit, Param::read_index_step1, Param::target_mismatch_count1);
        read_mapper.template index_and_mapping<false>(Param::k2, Param::ref_index_step2, Param::bucket_limit, Param::read_index_step2, Param::target_mismatch_count2);
        read_mapper.template index_and_mapping<true>(Param::k2, Param::ref_index_step2, Param::bucket_limit, Param::read_index_step2, Param::target_mismatch_count2);
        read_mapper.get_result(unmatch_reads_id_host, unmapping_reads_count, unmapping_N_reads_id, unmapping_N_reads_count, map_record, map_mis_cnt);
    }
    LOCK_END
    cudaHostUnregister(N_reads_db_host.data());

    auto unmapping_ref = std::make_unique<RefRecord>();
    auto unmapping_N_ref = std::make_unique<RefRecord>();

    auto unmapping_ref_construct = std::async(std::launch::async, [&] {
        if (unmapping_reads_count == 0) return;
        uint64_t * reads_db_buffer;
        uint32_t * next_host, * prev_host;
        uint16_t * offset_host;
        LOCK_START
        {
            cudaSetDevice(device_id);
            uint64_t * reads_db;
            uint32_t * prefix_reads_id;
            uint32_t * suffix_reads_id;
            uint32_t * prev_device;
            uint32_t * next_device;
            uint16_t * offset_device;

            gpuErrorCheck(cudaMalloc((void**) &reads_db, unmapping_reads_count * read_unit_size * sizeof(uint64_t)));
            gpuErrorCheck(cudaMalloc((void**) &prefix_reads_id, unmapping_reads_count * sizeof(uint32_t)));
            gpuErrorCheck(cudaMalloc((void**) &suffix_reads_id, unmapping_reads_count * sizeof(uint32_t)));
            gpuErrorCheck(cudaMalloc((void**) &prev_device, unmapping_reads_count * sizeof(uint32_t)));
            gpuErrorCheck(cudaMalloc((void**) &next_device, unmapping_reads_count * sizeof(uint32_t)));
            gpuErrorCheck(cudaMalloc((void**) &offset_device, unmapping_reads_count * sizeof(uint16_t)));
            thrust::sequence(thrust::device, prefix_reads_id, prefix_reads_id + unmapping_reads_count, 0);
            thrust::sequence(thrust::device, suffix_reads_id, suffix_reads_id + unmapping_reads_count, 0);
            thrust::uninitialized_fill(thrust::device, prev_device, prev_device + unmapping_reads_count, INVALID);
            thrust::uninitialized_fill(thrust::device, next_device, next_device + unmapping_reads_count, INVALID);
            thrust::uninitialized_fill(thrust::device, offset_device, offset_device + unmapping_reads_count, INVALOFF);

            cudaMallocHost((void**) &reads_db_buffer, unmapping_reads_count * read_unit_size * sizeof(uint64_t));
#pragma omp parallel for
            for (size_t i = 0; i < unmapping_reads_count; ++i) {
#pragma omp simd
                for (size_t j = 0; j < read_unit_size; ++j) {
                    reads_db_buffer[i * read_unit_size + j] = reads_db_host[unmatch_reads_id_host[i] * read_unit_size + j];
                }
            }
            cudaMemcpy(reads_db, reads_db_buffer, unmapping_reads_count * read_unit_size * sizeof(uint64_t), cudaMemcpyHostToDevice);

            full_match_gpu(reads_db, prefix_reads_id, suffix_reads_id, unmapping_reads_count,
                                           next_device, prev_device, offset_device, param, device_id);

            cudaMallocHost((void**) &prev_host, unmapping_reads_count * sizeof(uint32_t));
            cudaMallocHost((void**) &next_host, unmapping_reads_count * sizeof(uint32_t));
            cudaMallocHost((void**) &offset_host, unmapping_reads_count * sizeof(uint16_t));
            cudaMemcpy(next_host, next_device, unmapping_reads_count * sizeof(uint32_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(prev_host, prev_device, unmapping_reads_count * sizeof(uint32_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(offset_host, offset_device, unmapping_reads_count * sizeof(uint16_t), cudaMemcpyDeviceToHost);

            cudaFree(reads_db);
            cudaFree(prefix_reads_id);
            cudaFree(suffix_reads_id);
            cudaFree(prev_device);
            cudaFree(next_device);
            cudaFree(offset_device);
        }
        LOCK_END

        std::vector<uint32_t> unmapping_ref_reads_id(unmapping_reads_count);
        std::iota(unmapping_ref_reads_id.begin(), unmapping_ref_reads_id.end(), 0);
        unmapping_ref->assembly(reads_db_buffer, unmapping_ref_reads_id.data(), unmapping_reads_count, next_host, prev_host, offset_host, unmapping_reads_count, param.read_len, param.read_unit_size);
        cudaFreeHost(reads_db_buffer);
        cudaFreeHost(next_host);
        cudaFreeHost(prev_host);
        cudaFreeHost(offset_host);

#pragma omp parallel for
        for (size_t i = 0; i < unmapping_ref->record.size(); ++i) {
            auto id = unmapping_ref->record[i].read_id;
            unmapping_ref->record[i].read_id = unmatch_reads_id_host[id];
        }
        cudaFreeHost(unmatch_reads_id_host);
    });

    auto unmapping_N_ref_construct = std::async(std::launch::async, [&] {
        if (unmapping_N_reads_count == 0) return ;
        full_match_cpu(N_reads_db_host.data(), unmapping_N_reads_id, unmapping_N_reads_count, N_reads_id.size(), param, unmapping_N_ref.get());
        cudaFreeHost(unmapping_N_reads_id);
    });

    std::vector<uint64_t> id_to_pos;
    if (param.is_preserve_order) id_to_pos.resize(reads_count);
    // std::string pe_flag;
    std::vector<uint32_t> pe_reads_order;
    if (!param.is_preserve_order && param.is_paired_end) pe_reads_order.reserve(reads_count); // pe_flag.reserve(reads_count);
    {
        auto & record = ref->record;
        auto record_size = record.size();
        record.insert(record.end(), map_record.begin(), map_record.end());
        map_record.clear(); map_record.shrink_to_fit();
        std::vector<uint8_t> mis_cnt(record.size(), 0);
        for (size_t i = 0; i < map_mis_cnt.size(); ++i) mis_cnt[record_size + i] = map_mis_cnt[i];
        map_mis_cnt.clear(); map_mis_cnt.shrink_to_fit();

        if (!param.is_preserve_order) {
            thrust::sort_by_key(thrust::omp::par, record.begin(), record.end(), mis_cnt.begin(),
                                [](const AlignmentRecord &a, const AlignmentRecord &b) {
                return a.pos < b.pos;
            });
        } else {
            thrust::sort_by_key(thrust::omp::par, record.begin(), record.end(), mis_cnt.begin(),
                                [&](const AlignmentRecord &a, const AlignmentRecord &b) {
                uint32_t a_id = a.read_id, b_id = b.read_id;
                if (a.is_contain_N) a_id = N_reads_id[a_id];
                if (b.is_contain_N) b_id = N_reads_id[b_id];
                return a_id < b_id;
            });
        }

        uint64_t last_pos = 0;
        size_t mismatch_base_context[256][256] = {0};
        std::string mismatch_base_context_stream;
        std::string mismatch_base_stream;
        std::string strand_id_stream;
        std::vector<uint16_t> reads_off_stream;
        std::vector<std::vector<uint16_t>> mismatch_offset_stream(param.max_mismatch_count + 1);

        strand_id_stream.reserve(record.size());
        if (!param.is_preserve_order) reads_off_stream.reserve(record.size());

        for (size_t i = 0; i < record.size(); ++i) {
            const uint32_t read_id = record[i].read_id;
            const uint64_t is_contain_N = record[i].is_contain_N;
            const uint64_t strand_id = record[i].strand_id;
            const uint64_t pos = record[i].pos;
            const uint8_t mismatch_count = mis_cnt[i];
            if (param.is_preserve_order) {
                if (is_contain_N) {
                    id_to_pos[N_reads_id[read_id]] = pos;
                } else {
                    id_to_pos[read_id] = pos;
                }
            } else {
                if (param.is_paired_end) {
                    if (is_contain_N) {
                        // pe_flag += (char) ((uint32_t)'0' + N_reads_id[read_id] % 2);
                        pe_reads_order.push_back(N_reads_id[read_id]);
                    } else {
                        // pe_flag += (char) ((uint32_t)'0' + read_id % 2);
                        pe_reads_order.push_back(read_id);
                    }
                }
                reads_off_stream.push_back(pos - last_pos);
                last_pos = pos;
            }
            if (strand_id == 0) strand_id_stream += "0"; else strand_id_stream += "1";

            if (mismatch_count > 0) {
                std::vector<uint16_t> mismatch_offset;
                if (!is_contain_N) {
                    for (size_t p = 0; p < param.read_len; ++p) {
                        uint64_t base_offset = (strand_id == 0) ? (p) : (param.read_len - 1 - p);
                        char read_base = extract_base_cpu(reads_db_host, read_id * read_unit_size, base_offset);
                        if (strand_id == 1) read_base = reverse_dna_base_cpu(read_base);
                        if (ref->ref_string[pos + p] != read_base) {
                            // mismatch_count++;
                            mismatch_base_context_stream += ref->ref_string[pos + p];
                            mismatch_base_stream += read_base;
                            mismatch_base_context[ref->ref_string[pos + p]][read_base]++;
                            mismatch_offset.push_back(p);
                        }
                    }
                } else {
                    for (size_t p = 0; p < param.read_len; ++p) {
                        uint64_t base_offset = (strand_id == 0) ? (read_id * param.read_len + p) : (read_id * param.read_len + param.read_len - 1 - p);
                        char read_base = N_reads_db_host[base_offset];
                        if (strand_id == 1) read_base = reverse_dna_base_cpu(read_base);
                        if (ref->ref_string[pos + p] != read_base) {
                            // mismatch_count++;
                            mismatch_base_context_stream += ref->ref_string[pos + p];
                            mismatch_base_stream += read_base;
                            mismatch_base_context[ref->ref_string[pos + p]][read_base]++;
                            mismatch_offset.push_back(p);
                        }
                    }
                }

                if (mismatch_offset[0] < (param.read_len - 1 - mismatch_offset[mismatch_count - 1])) { // left to right
                    mismatch_offset_stream[mismatch_count].push_back(0);
                    mismatch_offset_stream[mismatch_count].push_back(mismatch_offset[0]);
                    if (mismatch_count > 1) {
                        if (mismatch_count == 2) {
                            mismatch_offset_stream[mismatch_count].push_back(mismatch_offset[1] - mismatch_offset[0] - 1);
                        } else {
                            uint32_t min_off = std::numeric_limits<uint32_t>::max();
                            for (size_t m = 1; m < mismatch_count; ++m) {
                                if ((mismatch_offset[m] - mismatch_offset[m - 1]) < min_off) {
                                    min_off = mismatch_offset[m] - mismatch_offset[m - 1];
                                }
                            }
                            mismatch_offset_stream[mismatch_count].push_back(min_off - 1);
                            for (size_t m = 1; m < mismatch_count; ++m) {
                                mismatch_offset_stream[mismatch_count].push_back(mismatch_offset[m] - mismatch_offset[m - 1] - min_off);
                            }
                        }
                    }
                } else { // right to left
                    mismatch_offset_stream[mismatch_count].push_back(1);
                    mismatch_offset_stream[mismatch_count].push_back(param.read_len - 1 - mismatch_offset[mismatch_count - 1]);
                    if (mismatch_count > 1) {
                        if (mismatch_count == 2) {
                            mismatch_offset_stream[mismatch_count].push_back(mismatch_offset[1] - mismatch_offset[0] - 1);
                        } else {
                            uint32_t min_off = std::numeric_limits<uint32_t>::max();
                            for (size_t m = mismatch_count - 1; m > 0; --m) {
                                if ((mismatch_offset[m] - mismatch_offset[m - 1]) < min_off) {
                                    min_off = mismatch_offset[m] - mismatch_offset[m - 1];
                                }
                            }
                            mismatch_offset_stream[mismatch_count].push_back(min_off - 1);
                            for (size_t m = mismatch_count - 1; m > 0; --m) {
                                mismatch_offset_stream[mismatch_count].push_back(mismatch_offset[m] - mismatch_offset[m - 1] - min_off);
                            }
                        }
                    }
                }
            }
        }
        record.clear(); record.shrink_to_fit();

        auto mismatch_base_recode = [&] (char ctx, std::array<char, 4> dest) {
            std::array<size_t, 4> index = {0,1,2,3};
            std::sort(index.begin(), index.end(), [&](size_t a, size_t b) {
                return mismatch_base_context[ctx][dest[a]] > mismatch_base_context[ctx][dest[b]];
            });
            for (size_t i = 0; i < index.size(); ++i) {
                mismatch_base_context[ctx][dest[index[i]]] = i;
            }
        };
        mismatch_base_recode('A', {'C','G','N','T'});
        mismatch_base_recode('T', {'A','C','G','N'});
        mismatch_base_recode('C', {'A','G','N','T'});
        mismatch_base_recode('G', {'A','C','N','T'});

        auto mismatch_base_context_store = [&] (char ctx, std::array<char, 4> dest) {
            char value[4];
            for (char ch : dest) {
                value[mismatch_base_context[ctx][ch]] = ch;
            }
            output_file.write(value, 4);
        };
        mismatch_base_context_store('A', {'C','G','N','T'});
        mismatch_base_context_store('T', {'A','C','G','N'});
        mismatch_base_context_store('C', {'A','G','N','T'});
        mismatch_base_context_store('G', {'A','C','N','T'});

        std::ofstream mismatch_base_lzma (working_path / "mismatch_base.lzma");
        std::ofstream mismatch_count_lzma (working_path / "mismatch_count.lzma");
        std::ofstream read_off_lzma (working_path / "read_off.lzma");
#pragma omp parallel for
        for (size_t i = 0; i < mismatch_base_stream.size(); ++i) {
            mismatch_base_stream[i] = (char) mismatch_base_context[mismatch_base_context_stream[i]][mismatch_base_stream[i]];
        }
        mismatch_base_context_stream.clear(); mismatch_base_context_stream.shrink_to_fit();
        writeCompressed(mismatch_base_lzma, mismatch_base_stream, PPMD7_CODER, 3, 2, COMPRESSION_ESTIMATION_MIS_SYM);
        mismatch_base_stream.clear(); mismatch_base_stream.shrink_to_fit();
        writeCompressed(mismatch_count_lzma, reinterpret_cast<const char*>(mis_cnt.data()), mis_cnt.size(), PPMD7_CODER, 3, 2, 1/*COMPRESSION_ESTIMATION_MIS_CNT*/);
        mis_cnt.clear(); mis_cnt.shrink_to_fit();
        if (!param.is_preserve_order) {
            writeCompressed(read_off_lzma, reinterpret_cast<const char *>(reads_off_stream.data()), reads_off_stream.size() * 2, PPMD7_CODER, 3, 3, 1);
            reads_off_stream.clear(); reads_off_stream.shrink_to_fit();
        }
        {
            std::ofstream strand_id_txt(working_path / "strand_id.txt");
            strand_id_txt << strand_id_stream;
        }
        {
            bsc_compress((working_path / "strand_id.txt").c_str(), (working_path / "strand_id.bsc").c_str());
            fs::remove(working_path / "strand_id.txt");
            strand_id_stream.clear();
            strand_id_stream.shrink_to_fit();
        }
        for (size_t i = 1; i <= param.max_mismatch_count; ++i) {
            if (param.read_len <= 128) {
                mismatch_offset_compress<128>(mismatch_offset_stream[i], i, working_path / ("mismatch_off_" + std::to_string(i)));
            } else if (param.read_len <= 256) {
                mismatch_offset_compress<256>(mismatch_offset_stream[i], i, working_path / ("mismatch_off_" + std::to_string(i)));
            } else if (param.read_len <= 384) {
                mismatch_offset_compress<384>(mismatch_offset_stream[i], i, working_path / ("mismatch_off_" + std::to_string(i)));
            } else if (param.read_len <= 512) {
                mismatch_offset_compress<512>(mismatch_offset_stream[i], i, working_path / ("mismatch_off_" + std::to_string(i)));
            } else {
                mismatch_offset_compress<65536>(mismatch_offset_stream[i], i, working_path / ("mismatch_off_" + std::to_string(i)));
            }
        }
    }

    unmapping_ref_construct.get();
    unmapping_N_ref_construct.get();

    if (!param.is_preserve_order) {
        std::vector<uint16_t> reads_off_stream; // may be empty
        reads_off_stream.reserve(unmapping_ref->record.size());
        uint64_t last_pos = 0;
        for (size_t i = 0; i < unmapping_ref->record.size(); ++i) {
            uint32_t read_id = unmapping_ref->record[i].read_id ;
            if (param.is_paired_end) {
                // pe_flag += (char)((uint32_t)'0' + read_id % 2);
                pe_reads_order.push_back(read_id);
            }
            reads_off_stream.push_back(param.read_len - (unmapping_ref->record[i].pos - last_pos));
            last_pos = unmapping_ref->record[i].pos;
        }
        std::ofstream unmapping_reads_off_lzma(working_path / "unmapping_read_off.lzma");
        writeCompressed(unmapping_reads_off_lzma, reinterpret_cast<const char*>(reads_off_stream.data()), reads_off_stream.size() * 2, PPMD7_CODER, 3, 3, 1);
    } else {
        for (size_t i = 0; i < unmapping_ref->record.size(); ++i) {
            id_to_pos[unmapping_ref->record[i].read_id] = ref->ref_string.size() + unmapping_ref->record[i].pos;
        }
    }

    if (!param.is_preserve_order) {
        std::vector<uint16_t> N_reads_off_stream; // may be empty
        N_reads_off_stream.reserve(unmapping_N_ref->record.size());
        uint64_t last_pos = 0;
        for (size_t i = 0; i < unmapping_N_ref->record.size(); ++i) {
            uint32_t read_id = unmapping_N_ref->record[i].read_id;
            if (param.is_paired_end) {
                // pe_flag += (char)((uint32_t)'0' + N_reads_id[read_id] % 2);
                pe_reads_order.push_back(N_reads_id[read_id]);
            }
            N_reads_off_stream.push_back(param.read_len - (unmapping_N_ref->record[i].pos - last_pos));
            last_pos = unmapping_N_ref->record[i].pos;
        }
        std::ofstream N_reads_off_lzma(working_path / "unmapping_N_read_off.lzma");
        writeCompressed(N_reads_off_lzma, reinterpret_cast<const char*>(N_reads_off_stream.data()), N_reads_off_stream.size() * 2, PPMD7_CODER, 3, 3, 1);
    } else {
        for (size_t i = 0; i < unmapping_N_ref->record.size(); ++i) {
            id_to_pos[N_reads_id[unmapping_N_ref->record[i].read_id]] = ref->ref_string.size() + unmapping_ref->ref_string.size() + unmapping_N_ref->record[i].pos;
        }
    }

    auto pe_reads_order_compress_future = std::async(std::launch::async, [&] {
        if (!param.is_preserve_order && param.is_paired_end) {
            std::ofstream pe_reads_order_comp(working_path / "pe_order.comp");
            paired_end_reads_order_compress(pe_reads_order, pe_reads_order_comp, working_path, param.flzma2_level, param.flzma2_thread_num);
            pe_reads_order.clear();
            pe_reads_order.shrink_to_fit();

//            {
//                std::ofstream pe_flag_txt(working_path / "pe_flag.txt");
//                pe_flag_txt << pe_flag;
//            }
//            {
//                bsc_compress((working_path / "pe_flag.txt").c_str(), (working_path / "pe_flag.bsc").c_str());
//                fs::remove(working_path / "pe_flag.txt");
//                pe_flag.clear();
//                pe_flag.shrink_to_fit();
//            }
        }
    });

    uint64_t joinedRefLength = (ref->ref_string.size() + unmapping_ref->ref_string.size() + unmapping_N_ref->ref_string.size());
    uint8_t isJoinRefLengthStd = joinedRefLength <= UINT32_MAX;
    if (param.is_preserve_order) {
        output_file.write(reinterpret_cast<const char*>(&isJoinRefLengthStd), sizeof(uint8_t));
    }
    auto id_to_pos_compress_future = std::async(std::launch::async, [&] {
        if (param.is_preserve_order) {
            if (!param.is_paired_end) {
                std::ofstream id_pos_file(working_path / "id_pos.bin");
                if (isJoinRefLengthStd) {
                    std::vector<uint32_t> id_to_pos_32(id_to_pos.begin(), id_to_pos.end());
                    PgSAHelpers::writeArray(id_pos_file, (void*)id_to_pos_32.data(), id_to_pos_32.size() * 4);
                } else {
                    PgSAHelpers::writeArray(id_pos_file, (void*)id_to_pos.data(), id_to_pos.size() * 8);
                }
                lzma2::lzma2_compress((working_path / "id_pos.bin").c_str(), (working_path / "id_pos.comp").c_str(), param.flzma2_level, param.flzma2_thread_num);
            } else {
                std::ofstream id_pos_comp(working_path / "id_pos.comp");
                if (isJoinRefLengthStd) {
                    paired_end_id_to_pos_compression<uint32_t>(id_to_pos, id_pos_comp, working_path, param.flzma2_level, param.flzma2_thread_num);
                } else {
                    paired_end_id_to_pos_compression<uint64_t>(id_to_pos, id_pos_comp, working_path, param.flzma2_level, param.flzma2_thread_num);
                }
            }
            id_to_pos.clear();
            id_to_pos.shrink_to_fit();
        }
    });

    {
        uint8_t isRefLengthStd = ref->ref_string.size() <= UINT32_MAX;
        output_file.write(reinterpret_cast<const char*>(&isRefLengthStd), sizeof(uint8_t));

        std::string unmapRefMapOff, unmapRefMapLen;
        std::string unmapNRefMapOff, unmapNRefMapLen;
        std::string refMapOff, refMapLen;
        LOCK_START
        {
            cudaSetDevice(device_id);
            RefMatcher matcher(ref->ref_string.size(), device_id, param);
            matcher.init(ref->binary_ref_string, ref->binary_ref_size);
            if (!unmapping_ref->ref_string.empty()) {
                matcher.match(unmapping_ref.get(), unmapRefMapOff, unmapRefMapLen, false, false);
            }
            if (!unmapping_N_ref->ref_string.empty()) {
                matcher.match(unmapping_N_ref.get(), unmapNRefMapOff, unmapNRefMapLen, true, false);
            }
            matcher.match(ref.get(), refMapOff, refMapLen, false, true);
            // matcher.finish_match();
        }
        LOCK_END

        double estimated_ref_offset_ratio = simpleUintCompressionEstimate(ref->ref_string.size(), isRefLengthStd ? UINT32_MAX : UINT64_MAX);
        const int ref_offset_data_period_code = isRefLengthStd ? 2 : 3;
        auto create_map_file = [&]
                (const std::string &map_off, const std::string &map_len, const std::string &name) {
            std::ofstream map_off_lzma(working_path / (name + "_off.map.lzma"));
            std::ofstream map_len_lzma(working_path / (name + "_len.map.lzma"));
            writeCompressed(map_off_lzma, map_off.data(), map_off.size(), LZMA_CODER, 3, ref_offset_data_period_code, estimated_ref_offset_ratio);
            writeCompressed(map_len_lzma, map_len.data(), map_len.size(), LZMA_CODER, 3, 0, 1);
        };
        create_map_file(refMapOff, refMapLen, "ref");
        create_map_file(unmapRefMapOff, unmapRefMapLen, "unmap_ref");
        create_map_file(unmapNRefMapOff, unmapNRefMapLen,   "unmap_n_ref");
    }

    {
        uint64_t ref_size = ref->ref_string.size();
        uint64_t unmapping_ref_size = unmapping_ref->ref_string.size();
        uint64_t unmapping_N_ref_size = unmapping_N_ref->ref_string.size();
        output_file.write(reinterpret_cast<const char*>(&ref_size), sizeof(uint64_t));
        output_file.write(reinterpret_cast<const char*>(&unmapping_ref_size), sizeof(uint64_t));
        output_file.write(reinterpret_cast<const char*>(&unmapping_N_ref_size), sizeof(uint64_t));

        ref->ref_string.append(unmapping_ref->ref_string); unmapping_ref->ref_string.clear(); unmapping_ref->ref_string.shrink_to_fit();
        ref->ref_string.append(unmapping_N_ref->ref_string);  unmapping_N_ref->ref_string.clear(); unmapping_N_ref->ref_string.shrink_to_fit();

        size_t compLen = 0;
        char* var_len_encode_seq = Compress(compLen, ref->ref_string.data(), ref->ref_string.size(),
                                            VARLEN_DNA_CODER, 3, PgSAHelpers::VarLenDNACoder::getCoderParam
                                                    (PgSAHelpers::VarLenDNACoder::STATIC_CODES_CODER_PARAM,
                                                     PgSAHelpers::VarLenDNACoder::AG_EXTENDED_CODES_ID), 1);
        {
            std::ofstream ref_variable_length_encode_file(working_path / "ref_var_len_encode.bin");
            PgSAHelpers::writeArray(ref_variable_length_encode_file, (void*)var_len_encode_seq, compLen);
        }
        delete[] var_len_encode_seq;
        ref->ref_string.clear();
        ref->ref_string.shrink_to_fit();

        lzma2::lzma2_compress((working_path / "ref_var_len_encode.bin").c_str(), (working_path / "ref.lzma2").c_str(), param.flzma2_level, param.flzma2_thread_num);
    }

    auto read_file_and_output = [&](const fs::path &filename) {
        std::vector<uint8_t> data;
        read_vector_from_binary_file(data, filename);
        uint64_t size = data.size();
        output_file.write(reinterpret_cast<const char*>(&size), sizeof (uint64_t));
        output_file.write(reinterpret_cast<const char*>(data.data()), size);
    };

    pe_reads_order_compress_future.get();
    id_to_pos_compress_future.get();
    if (param.is_preserve_order) {
        read_file_and_output(working_path / "id_pos.comp");
    }
    read_file_and_output(working_path / "strand_id.bsc");
    read_file_and_output(working_path / "ref.lzma2");
    if (!param.is_preserve_order) {
        read_file_and_output(working_path / "read_off.lzma");
        read_file_and_output(working_path / "unmapping_N_read_off.lzma");
        read_file_and_output(working_path / "unmapping_read_off.lzma");
        if (param.is_paired_end) {
            // read_file_and_output(working_path / "pe_flag.bsc");
            read_file_and_output(working_path / "pe_order.comp");
        }
    }
    read_file_and_output(working_path / "ref_off.map.lzma");
    read_file_and_output(working_path / "ref_len.map.lzma");
    read_file_and_output(working_path / "unmap_ref_off.map.lzma");
    read_file_and_output(working_path / "unmap_ref_len.map.lzma");
    read_file_and_output(working_path / "unmap_n_ref_off.map.lzma");
    read_file_and_output(working_path / "unmap_n_ref_len.map.lzma");
    read_file_and_output(working_path / "mismatch_count.lzma");
    read_file_and_output(working_path / "mismatch_base.lzma");
    output_file.write(reinterpret_cast<const char*>(&param.max_mismatch_count), sizeof(uint8_t));
    for (size_t i = 1; i <= param.max_mismatch_count; ++i) {
        read_file_and_output(working_path / ("mismatch_off_" + std::to_string(i) + ".bin"));
    }

}

// template<size_t read_unit_size>
void process(const Param& param) {
    int device_count;
    gpuErrorCheck(cudaGetDeviceCount(&device_count));
    printf("GPU device count : %d\n", device_count);
    auto gpu_context_init_future = std::async(std::launch::async, [device_count]() {
#pragma omp parallel num_threads(device_count)
        {
            cudaSetDevice(omp_get_thread_num());
            cudaFree(nullptr);
        }
        printf("gpu context init finish\n");
    });

    uint64_t * reads_db_host;
    size_t reads_count = 0;
    uint8_t blocks_count = 0;
    size_t reads_db_host_capacity, fastq_bytes, block_bytes;
    if (!param.is_paired_end) {
        fastq_bytes = fs::file_size(param.f1_path);
        block_bytes = (size_t)((double) fastq_bytes * param.block_ratio);
        reads_db_host_capacity = (double) fastq_bytes * param.reads_data_ratio / param.read_len * param.read_unit_size * sizeof(uint64_t);
        printf("reads_db_host_capacity : %zu bytes\n", reads_db_host_capacity);
        cudaHostAlloc((void **) &reads_db_host, reads_db_host_capacity, cudaHostAllocPortable);
    } else {
        fastq_bytes = fs::file_size(param.f1_path);
        if (fastq_bytes != fs::file_size(param.f2_path)) throw std::runtime_error("paired-end file size isn't equal");
        block_bytes = (size_t)(2 * (double) fastq_bytes * param.block_ratio);
        reads_db_host_capacity = 2 * (double) fastq_bytes * param.reads_data_ratio / param.read_len * param.read_unit_size * sizeof(uint64_t);
        printf("reads_db_host_capacity : %zu bytes\n", reads_db_host_capacity);
        cudaHostAlloc((void **) &reads_db_host, reads_db_host_capacity, cudaHostAllocPortable);
    }

    std::vector<std::future<void>> block_compress_future;
    {
        cudaSetDevice(0);
        auto encode = [reads_db_host, read_unit_size = param.read_unit_size](char * buffer, uint8_t * contain_N_flags, const Param& param, size_t buffer_reads_count, size_t start_read_index) {
            cudaSetDevice(0);
            uint16_t read_len = param.read_len;
            uint64_t * binary_reads_buffer;
            cudaMalloc((void**) &binary_reads_buffer, buffer_reads_count * read_unit_size * sizeof(uint64_t));
            thrust::for_each(thrust::device,
                             thrust::counting_iterator<uint32_t>(0),
                             thrust::counting_iterator<uint32_t>(buffer_reads_count),
                             [buffer, contain_N_flags, binary_reads_buffer, read_len, read_unit_size = param.read_unit_size] __device__ (uint32_t i) {
                                 if (contain_N_flags[i]) return;
                                 uint64_t idx1 = 0, c;
                                 uint8_t idx2 = 0;
                                 for (size_t j = i * read_len; j < (i + 1) * read_len; ++j) {
                                     c = base_encode_table_mgpu[buffer[j]] & 0x03U;
                                     binary_reads_buffer[i * read_unit_size + idx1] |= (c << (62U - idx2));
                                     idx2 = (idx2 + 2U) & 0x3FU;
                                     if (idx2 == 0) {
                                         idx1++;
                                     }
                                 }
                             });
            cudaMemcpy(reads_db_host + start_read_index * read_unit_size, binary_reads_buffer, buffer_reads_count * read_unit_size * sizeof(uint64_t), cudaMemcpyDeviceToHost);
            cudaFree(binary_reads_buffer);
            cudaFree(buffer);
            cudaFree(contain_N_flags);
        };

        mio::mmap_source io_buffer, io_buffer_2;
        const size_t io_buffer_size = 50 * 1000 * 1000; // 50MB
        const size_t io_buffer_reads_size = io_buffer_size / 2;
        const size_t io_buffer_reads_count = io_buffer_reads_size / param.read_len;
        char * reads_db_buffer;
        cudaMallocHost((void**) &reads_db_buffer, io_buffer_reads_size);
        uint8_t * buffer_contain_N_flags;
        cudaMallocHost((void**) &buffer_contain_N_flags, io_buffer_reads_count);

        std::vector<std::future<void>> encode_futures;
        std::vector<uint32_t> reads_id, N_reads_id;
        std::vector<char> N_reads_db_host;
        size_t total_read_bytes = 0, last_total_read_bytes = 0, last_reads_count = 0;
        uint32_t global_reads_id = 0;

        while(total_read_bytes < fastq_bytes) {
            std::error_code error;
            if (!param.is_paired_end) {
                io_buffer.map(param.f1_path, total_read_bytes, std::min(io_buffer_size, fastq_bytes - total_read_bytes), error);
            } else {
                io_buffer.map(param.f1_path, total_read_bytes, std::min(io_buffer_size / 2, fastq_bytes - total_read_bytes), error);
                io_buffer_2.map(param.f2_path, total_read_bytes, std::min(io_buffer_size / 2, fastq_bytes - total_read_bytes), error);
            }
            auto buffer_end = io_buffer.size();
            total_read_bytes += buffer_end;
            if (total_read_bytes < fastq_bytes) {
                long first_line_len = 0;
                if (io_buffer[0] != '@') throw std::runtime_error("io_buffer doesn't start with @");
                while (io_buffer[first_line_len] != '\n') first_line_len++;
                auto new_buffer_end = buffer_end - first_line_len;
                int equal_count = 0;
                while(new_buffer_end > 0) {
                    if (io_buffer[new_buffer_end - 1] == '\n' && io_buffer[new_buffer_end] == '@') { // id or qual
                        for (size_t j = 0, k = new_buffer_end; j < first_line_len; j++, k++) {
                            if (io_buffer[j] == io_buffer[k]) {
                                equal_count++;
                                if (equal_count > 5) {
                                    goto get_end_finish;
                                }
                            } else {
                                equal_count = 0;
                            }
                        }
                    }
                    new_buffer_end--;
                }
                get_end_finish: ;
                total_read_bytes -= (buffer_end - new_buffer_end);
                buffer_end = new_buffer_end;
            }
            char base;
            long i = 0, j = 0;
            bool contain_N_flag = false;
            size_t buffer_reads_count = 0;

            while (i < buffer_end) {
                while(io_buffer[i++] != '\n');
                auto seq_start = i;
                while((base = io_buffer[i]) != '\n') {
                    if (base == 'N') contain_N_flag = true;
                    i++;
                }
                if (param.read_len != (i - seq_start)) throw std::runtime_error("read length isn't fix");

                if (contain_N_flag) {
                    N_reads_id.push_back(global_reads_id++);
                    N_reads_db_host.insert(N_reads_db_host.end(), io_buffer.data() + seq_start, io_buffer.data() + seq_start + param.read_len);
                } else {
                    reads_id.push_back(global_reads_id++);
                    std::memcpy(reads_db_buffer + buffer_reads_count * param.read_len, io_buffer.data() + seq_start, param.read_len);
                }
                buffer_contain_N_flags[buffer_reads_count++] = contain_N_flag;
                contain_N_flag = false;
                i++;
                while(io_buffer[i++] != '\n');
                while(io_buffer[i++] != '\n');

                if (param.is_paired_end) {
                    while(io_buffer_2[j++] != '\n');
                    seq_start = j;
                    while ((base = io_buffer_2[j]) != '\n') {
                        if (base == 'N') contain_N_flag = true;
                        j++;
                    }
                    if (param.read_len != (j - seq_start)) throw std::runtime_error("read length isn't fix");
                    if (contain_N_flag) {
                        N_reads_id.push_back(global_reads_id++);
                        N_reads_db_host.insert(N_reads_db_host.end(), io_buffer_2.data() + seq_start, io_buffer_2.data() + seq_start + param.read_len);
                    } else {
                        reads_id.push_back(global_reads_id++);
                        std::memcpy(reads_db_buffer + buffer_reads_count * param.read_len, io_buffer_2.data() + seq_start, param.read_len);
                    }
                    buffer_contain_N_flags[buffer_reads_count++] = contain_N_flag;
                    contain_N_flag = false;
                    j++;
                    while(io_buffer_2[j++] != '\n');
                    while(io_buffer_2[j++] != '\n');
                }
            }

            char * reads_db_buffer_device;
            uint8_t * buffer_contain_N_flags_device;
            cudaMalloc((void**) &reads_db_buffer_device, buffer_reads_count * param.read_len);
            cudaMalloc((void**) &buffer_contain_N_flags_device, buffer_reads_count * sizeof(uint8_t));
            cudaMemcpy(reads_db_buffer_device, reads_db_buffer, buffer_reads_count * param.read_len, cudaMemcpyHostToDevice);
            cudaMemcpy(buffer_contain_N_flags_device, buffer_contain_N_flags, buffer_reads_count * sizeof(uint8_t), cudaMemcpyHostToDevice);
            encode_futures.push_back(std::async(std::launch::async, encode, reads_db_buffer_device, buffer_contain_N_flags_device, param, buffer_reads_count, reads_count));
            reads_count += buffer_reads_count;
            if (reads_count * param.read_unit_size * sizeof(uint64_t) >= reads_db_host_capacity) {
                throw std::runtime_error("reads_db_host overflow");
            }

            size_t block_read_bytes = total_read_bytes - last_total_read_bytes;
            if (param.is_paired_end) block_read_bytes *= 2;
            if (block_read_bytes > block_bytes) {
                for (auto & f : encode_futures) f.get(); // force encode finish
                encode_futures.clear();
                block_compress_future.push_back(std::async(std::launch::async, block_compress,
                        reads_db_host + last_reads_count * param.read_unit_size, std::move(reads_id),
                        std::move(N_reads_db_host), std::move(N_reads_id), param, blocks_count++, device_count));
                if (blocks_count >= UINT8_MAX) {
                    throw std::runtime_error("too many blocks");
                }
                global_reads_id = 0;
                last_total_read_bytes = total_read_bytes;
                last_reads_count = reads_count;
            }
        }
        size_t block_read_bytes = total_read_bytes - last_total_read_bytes;
        if (block_read_bytes > 0) {
            for (auto & f : encode_futures) f.get();
            encode_futures.clear();
            block_compress_future.push_back(std::async(std::launch::async, block_compress,
                                                       reads_db_host + last_reads_count * param.read_unit_size, std::move(reads_id),
                                                       std::move(N_reads_db_host), std::move(N_reads_id), param, blocks_count++, device_count));
            if (blocks_count >= UINT8_MAX) {
                throw std::runtime_error("too many blocks");
            }
            global_reads_id = 0;
            last_total_read_bytes = total_read_bytes;
            last_reads_count = reads_count;
        }

        printf("reads count : %zu \n", reads_count);
        printf("block count : %u \n", blocks_count);
        cudaFreeHost(reads_db_buffer);
        cudaFreeHost(buffer_contain_N_flags);
    }

    fs::path output_path = fs::path(param.working_dir) / fs::path(param.output_name + archive_name_suffix);
    std::ofstream output_file(output_path);
    output_file.write(reinterpret_cast<const char*>(&param.is_preserve_order), sizeof(uint8_t));
    output_file.write(reinterpret_cast<const char*>(&param.is_paired_end), sizeof(uint8_t));
    output_file.write(reinterpret_cast<const char*>(&param.read_len), sizeof(uint16_t));
    output_file.write(reinterpret_cast<const char*>(&blocks_count), sizeof(uint8_t));

    printf("waiting block compress ...\n");
    for (size_t b_id = 0; b_id < block_compress_future.size(); ++b_id) {
        block_compress_future[b_id].get();
        std::vector<char> block_compressed_data;
        read_vector_from_binary_file(block_compressed_data, fs::path(param.working_dir) / (std::to_string(b_id) + block_archive_name_suffix));
        PgSAHelpers::writeArray(output_file, block_compressed_data.data(), block_compressed_data.size());
        fs::remove(fs::path(param.working_dir) / (std::to_string(b_id) + block_archive_name_suffix));
        printf("block %zu compress finish, compressed size : %zu\n", b_id, block_compressed_data.size());
    }
    cudaFreeHost(reads_db_host);

    output_file.flush();
    printf("archive size : %zu bytes \n", fs::file_size(output_path));
    int status = std::system(("rm -rf " + param.working_parent_path).c_str());
    if (status != 0) {
        printf("remove tmp directory fail\n");
    }
}

void compress(Param& param) {
    // size_t read_unit_size = 0;
    {
        size_t total_size = 0;
        std::ifstream input_fastq(param.f1_path);
        std::string line;
        std::getline(input_fastq, line);
        total_size += line.size();
        std::getline(input_fastq, line);
        total_size += line.size();
        param.read_len = line.size();
        param.read_unit_size = ceil<32>(param.read_len);
        std::getline(input_fastq, line);
        total_size += line.size();
        std::getline(input_fastq, line);
        total_size += line.size();
        printf("read length : %d\n", param.read_len);
        total_size += 4;
        param.reads_data_ratio = double(param.read_len + 1) / double(total_size);
    }

    param.max_off = param.read_len * 65 / 100;
    param.kmer_size = 32;
    if (param.read_len < 89) {
        param.kmer_size = 36 * param.read_len / 100;
    }
    param.kmer_index_pos = param.read_len - param.kmer_size - param.max_off;
    param.max_mismatch_count = param.read_len / 6;

    if (param.read_unit_size > 16) {
        printf("Read length must be less than 512");
        std::exit(0);
    }
    process(param);

//    switch(read_unit_size) {
//        case 1: // <= 32 bases
//            process<1>(param);
//            break;
//        case 2: // <= 64 bases
//            process<2>(param);
//            break;
//        case 3: // <= 96 bases
//            process<3>(param);
//            break;
//        case 4: // <= 128 bases
//            process<4>(param);
//            break;
//        case 5: // <= 160 bases
//            process<5>(param);
//            break;
//        case 6: // <= 192 bases
//            process<6>(param);
//            break;
//        case 7: // <= 224 bases
//            process<7>(param);
//            break;
//        case 8: // <= 256 bases
//            process<8>(param);
//            break;
//        case 9: // <= 288 bases
//            process<9>(param);
//            break;
//        case 10: // <= 320 bases
//            process<10>(param);
//            break;
//        case 11: // <= 352 bases
//            process<11>(param);
//            break;
//        case 12: // <= 384 bases
//            process<12>(param);
//            break;
//        case 13: // <= 416 bases
//            process<13>(param);
//            break;
//        case 14: // <= 448 bases
//            process<14>(param);
//            break;
//        case 15: // <= 480 bases
//            process<15>(param);
//            break;
//        case 16: // <= 512 bases
//            process<16>(param);
//            break;
//        default:
//            printf("Read length must be less than 512");
//            break;
//    }
}

