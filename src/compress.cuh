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
