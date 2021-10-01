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
    double block_ratio;
    int flzma2_level;
    int flzma2_thread_num;
    uint8_t is_preserve_order;
    uint8_t is_paired_end;
    uint16_t read_len;
    double reads_data_ratio;

    // error-free match argument
    size_t max_off;
    size_t kmer_size;
    size_t kmer_index_pos;

    // read mapping argument
    uint8_t max_mismatch_count;
    static constexpr size_t bucket_limit = 12;
    static constexpr size_t k1 = 23, ref_index_step1 = 5, read_index_step1 = 4, target_mismatch_count1 = 0;
    static constexpr size_t k2 = 13, ref_index_step2 = 3, read_index_step2 = 2, target_mismatch_count2 = 2;
};

#endif //CURC_PARAM_HPP
