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

#ifndef CURC_PARAM_HPP
#define CURC_PARAM_HPP
#include <string>

struct Param {
    // archive parameter
    std::string working_dir;
    std::string working_parent_path;
    std::string f1_path;
    std::string f2_path;
    std::string output_name;
    std::string random_block_prefix_id;
    double block_ratio;
    int flzma2_level;
    int flzma2_thread_num;
    uint8_t is_preserve_order;
    uint8_t is_paired_end;
    uint16_t read_len;
    size_t read_unit_size;
    double reads_data_ratio;
    uint64_t decode_buffer_size;

    // error-free match argument
    size_t max_off;
    size_t kmer_size;
    size_t kmer_index_pos;

    // read mapping argument
    uint8_t max_mismatch_count;
    static constexpr size_t bucket_limit = 12;
    static constexpr size_t k1 = /*23*/28, ref_index_step1 = 5, read_index_step1 = /*4*/2, target_mismatch_count1 = /*0*/1;
    static constexpr size_t k2 = 13, ref_index_step2 = 3, read_index_step2 = 2, target_mismatch_count2 = 2;
};

#endif //CURC_PARAM_HPP
