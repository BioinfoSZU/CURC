/**
Software License:
-----------------

Copyright (c) 2021 Zexuan Zhu<zhuzx@szu.edu.cn>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

PgRC:
-----------------

The PgRC is an in-memory algorithm for compressing the DNA stream of FASTQ
datasets, based on the idea of building an approximation of the shortest
common superstring over high-quality reads. CURC is based on the architecture
of PgRC and also uses parts of PgRC codes in backend encoding
(mainly variable-length encoding and pairing data encoding).

The PgRC copyright is as follows:

Copyright (c) 2020 Tomasz M. Kowalski, Szymon Grabowski All Rights Reserved.

See also the PgRC web site:
  https://github.com/kowallus/PgRC for more information.

*/

#include "preprocess.hpp"
#include <bits/stdc++.h>
#include <mio/mio.hpp>
#include <experimental/filesystem>
#include <lzma/VarLenDNACoder/LzmaLib.h>
#include "compress.cuh"
#include "constant.hpp"
#include "omp.h"

namespace fs = std::experimental::filesystem;

template<typename T>
static inline void read_vector_from_binary_file2(std::vector<T> &data, const fs::path &filename) {
    if (!fs::exists(filename)) throw std::runtime_error(filename.string() + " isn't exist");
    std::ifstream file;
    file.open(filename, std::ios::in | std::ios::binary);
    file.seekg(0, std::ios::end);
    auto size = file.tellg();
    file.seekg(0);
    data.resize(size / sizeof(T));
    file.read(reinterpret_cast<char *>(data.data()), size);
}

void process(const Param& param) {
    int device_count;
    get_device_count_and_init(device_count);
    printf("GPU device count : %d\n", device_count);

    uint64_t * reads_db_host;
    size_t reads_count = 0;
    size_t total_reads_count = 0;
    uint8_t blocks_count = 0;
    size_t reads_db_host_capacity, fastq_bytes, block_bytes;
    if (!param.is_paired_end) {
        fastq_bytes = fs::file_size(param.f1_path);
        block_bytes = (size_t)((double) fastq_bytes * param.block_ratio);
        reads_db_host_capacity = (double) (block_bytes + 50 * 1000 * 1000) * param.reads_data_ratio / param.read_len * param.read_unit_size * sizeof(uint64_t);
        gpu_malloc_host((void**)&reads_db_host, reads_db_host_capacity);
    } else {
        fastq_bytes = fs::file_size(param.f1_path);
        if (fastq_bytes != fs::file_size(param.f2_path)) throw std::runtime_error("paired-end file size isn't equal");
        block_bytes = (size_t)(2 * (double) fastq_bytes * param.block_ratio);
        reads_db_host_capacity = (double) (block_bytes + 50 * 1000 * 1000) * param.reads_data_ratio / param.read_len * param.read_unit_size * sizeof(uint64_t);
        gpu_malloc_host((void**)&reads_db_host, reads_db_host_capacity);
    }

    std::vector<std::future<void>> block_compress_future;
    {
        mio::mmap_source io_buffer, io_buffer_2;
        const size_t io_buffer_size = 50 * 1000 * 1000; // 50MB
        const size_t io_buffer_reads_size = io_buffer_size / 2;
        const size_t io_buffer_reads_count = io_buffer_reads_size / param.read_len;
        char * reads_db_buffer;
        gpu_malloc_host((void**) &reads_db_buffer, io_buffer_reads_size);

        std::vector<std::future<void>> encode_futures;
        std::vector<uint32_t> reads_id, N_reads_id;
        std::vector<char> N_reads_db_host;
        size_t total_read_bytes = 0, last_total_read_bytes = 0; //, last_reads_count = 0;
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

            uint64_t c = 0;
            uint8_t idx2 = 0;

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
                }
                std::memcpy(reads_db_buffer + buffer_reads_count * param.read_len, io_buffer.data() + seq_start, param.read_len);

                buffer_reads_count++;
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
                    }
                    std::memcpy(reads_db_buffer + buffer_reads_count * param.read_len, io_buffer_2.data() + seq_start, param.read_len);

                    buffer_reads_count++;
                    contain_N_flag = false;
                    j++;
                    while(io_buffer_2[j++] != '\n');
                    while(io_buffer_2[j++] != '\n');
                }
            }

            char * reads_db_buffer_device;

            gpu_malloc((void**) &reads_db_buffer_device, buffer_reads_count * param.read_len);
            gpu_mem_H2D(reads_db_buffer_device, reads_db_buffer, buffer_reads_count * param.read_len);
            encode_futures.emplace_back(std::async(std::launch::async, encode, reads_db_buffer_device, /*buffer_contain_N_flags_device,*/ param, buffer_reads_count, reads_count, reads_db_host, param.read_unit_size));
            reads_count += buffer_reads_count;
            if (reads_count * param.read_unit_size * sizeof(uint64_t) >= reads_db_host_capacity) {
                throw std::runtime_error("reads_db_host overflow");
            }

            size_t block_read_bytes = total_read_bytes - last_total_read_bytes;
            if (param.is_paired_end) block_read_bytes *= 2;
            if (block_read_bytes > block_bytes) {
                for (auto & f : encode_futures) f.get(); // force encode finish
                encode_futures.clear();
                block_compress_future.emplace_back(std::async(std::launch::async, block_compress,
                                                     reads_db_host,  // + last_reads_count * param.read_unit_size,
                                                     std::move(reads_id), std::move(N_reads_db_host), std::move(N_reads_id),
                                                     param, blocks_count++, device_count));

                if (blocks_count >= UINT8_MAX) {
                    throw std::runtime_error("too many blocks");
                }

                if (total_read_bytes < fastq_bytes) {
                    size_t unread_bytes = fastq_bytes - total_read_bytes;
                    if (param.is_paired_end) unread_bytes *= 2;
                    if (unread_bytes < block_bytes) {
                        reads_db_host_capacity = (double) unread_bytes * param.reads_data_ratio / param.read_len * param.read_unit_size * sizeof(uint64_t);
                    }
                    gpu_malloc_host((void **) &reads_db_host, reads_db_host_capacity);
                }

                global_reads_id = 0;
                last_total_read_bytes = total_read_bytes;

                total_reads_count += reads_count;
                reads_count = 0;
                // last_reads_count = reads_count;
            }
        }
        size_t block_read_bytes = total_read_bytes - last_total_read_bytes;
        if (block_read_bytes > 0) {
            for (auto & f : encode_futures) f.get();
            encode_futures.clear();
            block_compress_future.emplace_back(std::async(std::launch::async, block_compress,
                                                 reads_db_host,  // + last_reads_count * param.read_unit_size,
                                                 std::move(reads_id), std::move(N_reads_db_host), std::move(N_reads_id),
                                                 param, blocks_count++, device_count));

            if (blocks_count >= UINT8_MAX) {
                throw std::runtime_error("too many blocks");
            }
            global_reads_id = 0;
            last_total_read_bytes = total_read_bytes;

            total_reads_count += reads_count;
            reads_count = 0;
            // last_reads_count = reads_count;
        }

        printf("reads count : %zu \n", total_reads_count);
        printf("block count : %u \n", blocks_count);
        gpu_free_host(reads_db_buffer);
    }

    fs::path output_path = fs::path(param.working_dir) / fs::path(param.output_name + archive_name_suffix);
    std::ofstream output_file(output_path);
    output_file.write(reinterpret_cast<const char*>(&param.is_preserve_order), sizeof(uint8_t));
    output_file.write(reinterpret_cast<const char*>(&param.is_paired_end), sizeof(uint8_t));
    output_file.write(reinterpret_cast<const char*>(&param.read_len), sizeof(uint16_t));
    output_file.write(reinterpret_cast<const char*>(&blocks_count), sizeof(uint8_t));

    printf("waiting block compress ...\n");
    for (size_t b_id = 0; b_id < blocks_count; ++b_id) {
        block_compress_future[b_id].get();
        std::vector<char> block_compressed_data;
        read_vector_from_binary_file2(block_compressed_data, fs::path(param.working_dir) / (param.random_block_prefix_id + "_" + std::to_string(b_id) + block_archive_name_suffix));
        PgSAHelpers::writeArray(output_file, block_compressed_data.data(), block_compressed_data.size());
        fs::remove(fs::path(param.working_dir) / (param.random_block_prefix_id + "_" + std::to_string(b_id) + block_archive_name_suffix));
        printf("block %zu compress finish, compressed size : %zu\n", b_id, block_compressed_data.size());
    }

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
        param.read_unit_size = (param.read_len + 32 - 1) / 32;
        std::getline(input_fastq, line);
        total_size += line.size();
        std::getline(input_fastq, line);
        total_size += line.size();
        // printf("read length : %d\n", param.read_len);
        total_size += 4;
        param.reads_data_ratio = double(param.read_len + 1) / double(total_size);
    }

    param.max_off = param.read_len * 65 / 100;
    param.kmer_size = 32;
    if (param.read_len < 89) {
        param.kmer_size = 36 * param.read_len / 100;
    }
    param.kmer_index_pos = param.read_len - param.kmer_size - param.max_off;
    param.max_mismatch_count = param.read_len / 5;

    if (param.read_unit_size > 16) {
        printf("Read length must be less than 512");
        std::exit(0);
    }
    process(param);
}