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

#ifndef CURC_COMPRESS_CUH
#define CURC_COMPRESS_CUH

#include "Param.hpp"

void block_compress(uint64_t * reads_db_host,
                    std::vector<uint32_t> reads_id,
                    std::vector<char> N_reads_db_host,
                    std::vector<uint32_t> N_reads_id,
                    const Param& param,
                    uint8_t block_id,
                    int device_count);

void encode (char * buffer, const Param& param, size_t buffer_reads_count, size_t start_read_index,
             uint64_t * reads_db_host, size_t read_unit_size);

void get_device_count_and_init(int & device_count);

void gpu_malloc_host(void** ptr, size_t size) ;

void gpu_malloc(void** ptr, size_t size, int device_id = 0) ;

void gpu_free_host(void * ptr) ;

void gpu_free(void * ptr, int device_id = 0) ;

void gpu_mem_H2D(void * dest, const void * src, size_t size, int device_id = 0);

void gpu_mem_D2H(void * dest, const void * src, size_t size, int device_id = 0);

#endif //CURC_COMPRESS_CUH
