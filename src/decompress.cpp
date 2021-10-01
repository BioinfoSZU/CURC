#include "decompress.hpp"
#include <bits/stdc++.h>
#include <parallel/algorithm>
#include <experimental/filesystem>
#ifndef ABS
#define ABS std::abs
#endif
#include <fqzcomp/clr.cdr>
#include <fqzcomp/simple_model.h>
#include <lzma/VarLenDNACoder/LzmaLib.h>
#include "fast_lzma2_helper.hpp"

namespace fs = std::experimental::filesystem;

template<typename T>
std::vector<T> lzma_decompress_to_vector(const fs::path& path) {
    std::ifstream in(path);
    std::vector<T> ret;
    readCompressed(in, ret);
    return ret;
}

std::string lzma_decompress_to_string(const fs::path& path) {
    std::ifstream in(path);
    std::string ret;
    readCompressed(in, ret);
    return ret;
}

std::string extract_reference(const fs::path& working_path, const fs::path& path, size_t ref_size) {
    lzma2::lzma2_decompress(path.c_str(), (working_path / "ref_var_len_encode.bin").c_str());
    std::string var_len_encode_seq, ref_stream;
    {
        auto size = fs::file_size(working_path / "ref_var_len_encode.bin");
        var_len_encode_seq.resize(size);
        std::ifstream in(working_path / "ref_var_len_encode.bin");
        in.read(var_len_encode_seq.data(), size);
    }
    ref_stream.resize(ref_size);
    Uncompress((char *) ref_stream.data(), ref_stream.size(), var_len_encode_seq.data(), var_len_encode_seq.size(), VARLEN_DNA_CODER);
    return ref_stream;
}

std::string bsc_decompress(const fs::path& path) {
    std::string filename = path.stem().string();
    fs::path output_path(path.parent_path() / (filename + ".txt"));
    std::string bsc_command = std::string("./bsc") + " d "
                              + fs::absolute(path).c_str() + " "
                              + fs::absolute(output_path).c_str() + " "
                              + " > /dev/null";
    int status = std::system(bsc_command.c_str());
    if (status != 0) throw std::runtime_error("Error occurred during bsc decompress.");

    auto size = fs::file_size(output_path);
    std::string ret;
    ret.resize(size);
    std::ifstream output(output_path);
    output.read(ret.data(), size);
    return ret;
}

template<size_t max_read_len>
std::vector<uint16_t> mismatch_offset_decompress(const fs::path& path, size_t mismatch_count, size_t output_size) {
    std::vector<uint16_t> ret;
    ret.resize(output_size);

    std::vector<char> compressed_data(fs::file_size(path));
    std::ifstream input(path);
    input.read(compressed_data.data(), compressed_data.size());

    RangeCoder rc;
    rc.input(compressed_data.data());
    rc.StartDecode();

    SIMPLE_MODEL<2> flag;
    SIMPLE_MODEL<max_read_len / 2 + 1> first_value;
    SIMPLE_MODEL<max_read_len> value;

    size_t i = 0;
    while (i < ret.size()) {
        ret[i++] = flag.decodeSymbol(&rc);
        ret[i++] = first_value.decodeSymbol(&rc);
        if (mismatch_count > 1) {
            if (mismatch_count == 2) {
                ret[i++] = value.decodeSymbol(&rc);
            } else {
                ret[i++] = value.decodeSymbol(&rc);
                int j = mismatch_count - 1;
                while (j--) {
                    ret[i++] = value.decodeSymbol(&rc);
                }
            }
        }
    }

    rc.FinishDecode();
    return ret;
}

std::string restorePg(const std::string& destPg, std::string& srcPg,
                      std::istream& pgMapOffSrc, std::istream& pgMapLenSrc,
                      bool isPgLengthStd, bool srcIsDest) {
    std::string str;
    std::string & ret = srcIsDest ? srcPg : str;
    uint32_t target_match_len;
    PgSAHelpers::readUIntByteFrugal(pgMapLenSrc, target_match_len);
    uint64_t markPos = 0, posDest = 0;
    while((markPos = destPg.find('%', posDest)) != std::string::npos) {
        ret.append(destPg, posDest, markPos - posDest);
        posDest = markPos + 1;
        uint64_t matchSrcPos = 0;
        if (isPgLengthStd) {
            uint32_t tmp;
            PgSAHelpers::readValue<uint32_t>(pgMapOffSrc, tmp, false);
            matchSrcPos = tmp;
        } else {
            PgSAHelpers::readValue<uint64_t>(pgMapOffSrc, matchSrcPos, false);
        }
        uint64_t matchLength = 0;
        PgSAHelpers::readUIntByteFrugal(pgMapLenSrc, matchLength);
        matchLength += target_match_len;
        ret.append(PgSAHelpers::reverseComplement(srcPg.substr(matchSrcPos, matchLength)));
    }
    ret.append(destPg, posDest, destPg.length() - posDest);
    // printf("Restored Pg sequence of length: %zu\n", ret.length());
    return ret;
}

template<typename T>
void preserve_order_process(
        uint8_t max_mismatch_count, char mismatch_decoder_table[128][4],
        const std::string& hqPgSeq, const std::string& lqPgSeq, const std::string& nPgSeq,
        const std::string& mismatch_base_stream,
        const std::vector<uint8_t> & mismatch_count_stream,
        const std::vector<std::vector<uint16_t>> & mismatch_off_stream,
        const std::string& strand_id_stream,
        uint16_t read_len, bool is_paired_end, const fs::path& working_path,
        std::ofstream & o1, std::ofstream & o2) {
    std::vector<T> id_to_pos;
    if (!is_paired_end) {
        lzma2::lzma2_decompress((working_path / "id_pos.comp").c_str(), (working_path / "id_pos.bin").c_str());
        {
            auto size = fs::file_size(working_path / "id_pos.bin");
            id_to_pos.resize(size / sizeof(T));
            std::ifstream in(working_path / "id_pos.bin");
            in.read(reinterpret_cast<char *>(id_to_pos.data()), size);
        }
    } else {
        std::ifstream in(working_path / "id_pos.comp");
        auto prepare_file = [&](const std::string& filename) {
            uint64_t size;
            in.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
            std::vector<char> buffer(size);
            in.read(buffer.data(), size);
            std::ofstream out(working_path / filename);
            out.write(buffer.data(), size);
        };

        uint32_t total_reads_count;
        in.read(reinterpret_cast<char*>(&total_reads_count), 4);
        id_to_pos.reserve(total_reads_count);

        vector<uint8_t> offsetInUint16Flag;
        vector<uint8_t> offsetIsBaseFirstFlag;
        vector<uint16_t> offsetInUint16Value;
        vector<uint8_t> deltaInInt16Flag;
        vector<uint8_t> deltaIsBaseFirstFlag;
        vector<int16_t> deltaInInt16Value;
        vector<T> notBasePairPos;

        prepare_file("base_pair_pos.lzma2");
        lzma2::lzma2_decompress((working_path / "base_pair_pos.lzma2").c_str(), (working_path / "base_pair_pos.bin").c_str());
        {
            auto size = fs::file_size(working_path / "base_pair_pos.bin");
            id_to_pos.resize(size / sizeof(T));
            std::ifstream base_pair_pos_file(working_path / "base_pair_pos.bin");
            base_pair_pos_file.read(reinterpret_cast<char *>(id_to_pos.data()), size);
        }

        readCompressed(in, offsetInUint16Flag);
        readCompressed(in, offsetIsBaseFirstFlag);
        readCompressed(in, offsetInUint16Value);
        readCompressed(in, deltaInInt16Flag);
        readCompressed(in, deltaIsBaseFirstFlag);
        readCompressed(in, deltaInInt16Value);

        prepare_file("not_base_pair_pos.lzma2");
        lzma2::lzma2_decompress((working_path / "not_base_pair_pos.lzma2").c_str(), (working_path / "not_base_pair_pos.bin").c_str());
        {
            auto size = fs::file_size(working_path / "not_base_pair_pos.bin");
            notBasePairPos.resize(size / sizeof(T));
            std::ifstream not_base_pair_pos_file(working_path / "not_base_pair_pos.bin");
            not_base_pair_pos_file.read(reinterpret_cast<char *>(notBasePairPos.data()), size);
        }

        const uint32_t pairsCount = total_reads_count / 2;
        vector<uint32_t> bppRank;
        bppRank.reserve(pairsCount);
        for (uint32_t p = 0; p < pairsCount; p++)
            bppRank.push_back(p);
        __gnu_parallel::stable_sort(bppRank.begin(), bppRank.end(),
                                    [&](const uint32_t &idx1, const uint32_t &idx2) -> bool
                                    { return id_to_pos[idx1] < id_to_pos[idx2]; });

        id_to_pos.resize(total_reads_count);
        int64_t nbpPos = 0;
        int64_t offIdx = -1;
        int64_t delFlagIdx = -1;
        int64_t delIdx = -1;
        int64_t nbpPosIdx = -1;
        int64_t refPrev = 0;
        int64_t prev = 0;
        bool match = false;
        for (uint32_t i = 0; i < pairsCount; i++) {
            uint32_t p = bppRank[i];
            if (offsetInUint16Flag[i] == 1) {
                int64_t delta = offsetInUint16Value[++offIdx];
                if (offsetIsBaseFirstFlag[offIdx] == 0)
                    delta = -delta;
                nbpPos = id_to_pos[p] + delta;
            } else if (deltaInInt16Flag[++delFlagIdx]){
                int64_t delta = refPrev + deltaInInt16Value[++delIdx];
                refPrev = delta;
                prev = delta;
                if (deltaIsBaseFirstFlag[delIdx] == 0)
                    delta = -delta;
                nbpPos = id_to_pos[p] + delta;
                match = true;
            } else {
                nbpPos = notBasePairPos[++nbpPosIdx];
                int64_t delta = nbpPos - id_to_pos[p];
                if (delta < 0)
                    delta = -delta;
                if (!match || refPrev != prev)
                    refPrev = delta;
                match = false;
                prev = delta;
            }
            id_to_pos[pairsCount + p] = nbpPos;
        }
        /// 前面计算得到的 id_to_pos 左半部分是第一个文件的, 右半部分是第二个文件的, 这里需要进一步交叉开
        std::vector<T> right_part(pairsCount);
        for (size_t i = 0; i < pairsCount; ++i) {
            right_part[i] = id_to_pos[pairsCount + i];
            id_to_pos[pairsCount + i] = id_to_pos[i];
        }
        for (size_t i = 0; i < pairsCount; ++i) id_to_pos[i * 2] = id_to_pos[pairsCount + i];
        for (size_t i = 0; i < pairsCount; ++i) id_to_pos[i * 2 + 1] = right_part[i];
    }
    std::vector<size_t> mismatch_off_cur(129, 0);
    size_t mismatch_base_cur = 0;
    size_t mismatch_cnt_cur = 0;
    size_t strand_id_cur = 0;
    for (size_t i = 0; i < id_to_pos.size(); ++i) {
        uint64_t pos = id_to_pos[i];
        if (pos >= (hqPgSeq.size() + lqPgSeq.size())) { // n ref
            uint64_t match_pos = pos - (hqPgSeq.size() + lqPgSeq.size());
            if (!is_paired_end || i % 2 == 0) {
                o1 << nPgSeq.substr(match_pos, read_len) << "\n";
            } else {
                o2 << nPgSeq.substr(match_pos, read_len) << "\n";
            }
        } else if (pos >= hqPgSeq.size()) { // lq ref
            uint64_t match_pos = pos - hqPgSeq.size();
            if (!is_paired_end || i % 2 == 0) {
                o1 << lqPgSeq.substr(match_pos, read_len) << "\n";
            } else {
                o2 << lqPgSeq.substr(match_pos, read_len) << "\n";
            }
        } else {
            uint64_t match_pos = pos;
            uint8_t mismatch_cnt = mismatch_count_stream[mismatch_cnt_cur++];
            char strand_id = strand_id_stream[strand_id_cur++];
            std::string str = hqPgSeq.substr(match_pos, read_len);
            if (mismatch_cnt > 0) {
                for (size_t k = 0; k < mismatch_cnt; ++k) {
                    uint16_t mismatch_pos = mismatch_off_stream[mismatch_cnt][mismatch_off_cur[mismatch_cnt]++];
                    char origin_base = str[mismatch_pos];
                    char mismatch_base = mismatch_decoder_table[origin_base][mismatch_base_stream[mismatch_base_cur++]];
                    str[mismatch_pos] = mismatch_base;
                }
            }
            if (strand_id == '1') {
                PgSAHelpers::reverseComplementInPlace(str);
            }
            if (!is_paired_end || i % 2 == 0) {
                o1 << str << "\n";
            } else {
                o2 << str << "\n";
            }
        }
    }
}

void decompress_block(std::ifstream& input, std::ofstream & o1, std::ofstream & o2,
                      bool is_preserve_order, bool is_paired_end,
                      uint16_t read_len, const fs::path& working_path) {
    char mismatch_decoder_table[128][4];
    input.read(mismatch_decoder_table['A'], 4);
    input.read(mismatch_decoder_table['T'], 4);
    input.read(mismatch_decoder_table['C'], 4);
    input.read(mismatch_decoder_table['G'], 4);
    uint8_t isJoinPgLengthStd;
    if (is_preserve_order) {
        input.read(reinterpret_cast<char*>(&isJoinPgLengthStd), sizeof(uint8_t));
    }
    uint8_t isPgLengthStd;
    input.read(reinterpret_cast<char*>(&isPgLengthStd), sizeof(uint8_t));
    printf("isPgLengthStd : %u\n", isPgLengthStd);
    uint64_t HQ_ref_size, LQ_ref_size, N_ref_size;
    input.read(reinterpret_cast<char*>(&HQ_ref_size), sizeof(uint64_t));
    input.read(reinterpret_cast<char*>(&LQ_ref_size), sizeof(uint64_t));
    input.read(reinterpret_cast<char*>(&N_ref_size), sizeof(uint64_t));
    printf("HQ_ref_size : %zu; LQ_ref_size : %zu; N_ref_size : %zu\n", HQ_ref_size, LQ_ref_size, N_ref_size);
    auto prepare_file = [&](const std::string& filename) {
        uint64_t size;
        input.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
        std::vector<char> buffer(size);
        input.read(buffer.data(), size);
        std::ofstream out(working_path / filename);
        out.write(buffer.data(), size);
    };
    if (is_preserve_order) {
        prepare_file("id_pos.comp");
    }
    prepare_file("strand_id.bsc");
    prepare_file("ref.lzma2");
    if (!is_preserve_order) {
        prepare_file("read_off.lzma");
        prepare_file("N_read_off.lzma");
        prepare_file("LQ_read_off.lzma");
        if (is_paired_end) {
            prepare_file("pe_flag.bsc");
        }
    }
    prepare_file("hq_off.map.lzma");
    prepare_file("hq_len.map.lzma");
    prepare_file("lq_off.map.lzma");
    prepare_file("lq_len.map.lzma");
    prepare_file("n_off.map.lzma");
    prepare_file("n_len.map.lzma");
    prepare_file("mismatch_count.lzma");
    prepare_file("mismatch_base.lzma");
    uint8_t max_mismatch_count;
    input.read(reinterpret_cast<char*>(&max_mismatch_count), sizeof(uint8_t));
    printf("max mismatch count : %u\n", max_mismatch_count);
    for (size_t i = 1; i <= max_mismatch_count; ++i) {
        prepare_file("mismatch_off_" + std::to_string(i) + ".bin");
    }
    printf("extract archive finish\n");

    auto strand_id_stream = bsc_decompress(working_path / "strand_id.bsc");
    auto ref_stream = extract_reference(working_path, (working_path / "ref.lzma2"), HQ_ref_size + LQ_ref_size + N_ref_size);
    auto hqPgMapOff = lzma_decompress_to_string(working_path / "hq_off.map.lzma");
    auto hqPgMapLen = lzma_decompress_to_string(working_path / "hq_len.map.lzma");
    auto lqPgMapOff = lzma_decompress_to_string(working_path / "lq_off.map.lzma");
    auto lqPgMapLen = lzma_decompress_to_string(working_path / "lq_len.map.lzma");
    auto nPgMapOff = lzma_decompress_to_string(working_path / "n_off.map.lzma");
    auto nPgMapLen = lzma_decompress_to_string(working_path / "n_len.map.lzma");
    auto mismatch_count_stream = lzma_decompress_to_vector<uint8_t>(working_path / "mismatch_count.lzma");
    auto mismatch_base_stream = lzma_decompress_to_string(working_path / "mismatch_base.lzma");
    std::vector<size_t> mismatch_count_stat(max_mismatch_count + 1, 0);
    std::vector<std::vector<uint16_t>> mismatch_off_stream(max_mismatch_count + 1);
    for (size_t i = 0; i < mismatch_count_stream.size(); ++i) {
        mismatch_count_stat[mismatch_count_stream[i]]++;
    }
    for (size_t i = 1; i <= max_mismatch_count; ++i) {
        size_t output_size;
        if (i == 1) output_size = 2 * mismatch_count_stat[i];
        else if (i == 2) output_size = 3 * mismatch_count_stat[i];
        else output_size = (i + 2) * mismatch_count_stat[i];

        std::vector<uint16_t> ret;
        if (read_len <= 128) {
            ret = mismatch_offset_decompress<128>(working_path / ("mismatch_off_" + std::to_string(i) + ".bin"), i, output_size);
        } else if (read_len <= 256) {
            ret = mismatch_offset_decompress<256>(working_path / ("mismatch_off_" + std::to_string(i) + ".bin"), i, output_size);
        } else if (read_len <= 384) {
            ret = mismatch_offset_decompress<384>(working_path / ("mismatch_off_" + std::to_string(i) + ".bin"), i, output_size);
        } else if (read_len <= 512) {
            ret = mismatch_offset_decompress<512>(working_path / ("mismatch_off_" + std::to_string(i) + ".bin"), i, output_size);
        } else {
            ret = mismatch_offset_decompress<65536>(working_path / ("mismatch_off_" + std::to_string(i) + ".bin"), i, output_size);
        }

        mismatch_off_stream[i].reserve(i * mismatch_count_stat[i]);
        std::vector<uint16_t> tmp(i);
        size_t j = 0;
        while (j < ret.size()) {
            bool flag = ret[j++];
            if (!flag) {
                tmp[0] = ret[j++];
                if (i > 1) {
                    if (i == 2) {
                        tmp[1] = ret[j++] + tmp[0] + 1;
                    } else {
                        uint32_t min_off = ret[j++] + 1;
                        for (size_t k = 1; k < i; ++k) {
                            tmp[k] = ret[j++] + tmp[k - 1] + min_off;
                        }
                    }
                }
            } else {
                tmp[i - 1] = read_len - 1 - ret[j++];
                if (i > 1) {
                    if (i == 2) {
                        tmp[0] = tmp[1] - ret[j++] - 1;
                    } else {
                        uint32_t min_off = ret[j++] + 1;
                        for (size_t k = i - 1; k > 0; --k) {
                            tmp[k - 1] = tmp[k] - ret[j++] - min_off;
                        }
                    }
                }
            }
            mismatch_off_stream[i].insert(mismatch_off_stream[i].end(), tmp.begin(), tmp.end());
        }
    }

    std::string nPgSeq(ref_stream, HQ_ref_size + LQ_ref_size);
    ref_stream.resize(HQ_ref_size + LQ_ref_size);
    std::string lqPgSeq(ref_stream, HQ_ref_size);
    ref_stream.resize(HQ_ref_size);
    ref_stream.shrink_to_fit();
    std::string hqPgSeq = std::move(ref_stream);
    std::istringstream hqPgOffSrc(hqPgMapOff), hqPgLenSrc(hqPgMapLen);
    std::istringstream lqPgOffSrc(lqPgMapOff), lqPgLenSrc(lqPgMapLen);
    std::istringstream nPgOffSrc (nPgMapOff),  nPgLenSrc(nPgMapLen);
    std::string hqPgSrc;
    hqPgSeq = restorePg(hqPgSeq, hqPgSrc, hqPgOffSrc, hqPgLenSrc, isPgLengthStd, true);
    lqPgSeq = restorePg(lqPgSeq, hqPgSeq, lqPgOffSrc, lqPgLenSrc, isPgLengthStd, false);
    nPgSeq  = restorePg(nPgSeq,  hqPgSeq, nPgOffSrc, nPgLenSrc, isPgLengthStd, false);

    if (!is_preserve_order) {
        auto read_off_stream = lzma_decompress_to_vector<uint16_t>(working_path / "read_off.lzma");
        auto N_read_off_stream = lzma_decompress_to_vector<uint16_t>(working_path / "N_read_off.lzma");
        auto LQ_read_off_stream = lzma_decompress_to_vector<uint16_t>(working_path / "LQ_read_off.lzma");
        std::string pe_flag;
        if (is_paired_end) {
            pe_flag = bsc_decompress(working_path / "pe_flag.bsc");
        }

        uint64_t pos = 0;
        char flag;
        std::vector<size_t> mismatch_off_cur(max_mismatch_count + 1, 0);
        size_t mismatch_base_cur = 0;
        for (size_t i = 0; i < read_off_stream.size(); ++i) {
            uint8_t mis_cnt = mismatch_count_stream[i];
            uint16_t off = read_off_stream[i];
            pos += off;
            if (is_paired_end) flag = pe_flag[i];
            char strand_id = strand_id_stream[i];
            std::string str = hqPgSeq.substr(pos, read_len);
            if (mis_cnt != 0) {
                for (size_t k = 0; k < mis_cnt; ++k) {
                    uint16_t mismatch_pos = mismatch_off_stream[mis_cnt][mismatch_off_cur[mis_cnt]++];
                    char origin_base = str[mismatch_pos];
                    char mismatch_base = mismatch_decoder_table[origin_base][mismatch_base_stream[mismatch_base_cur++]];
                    str[mismatch_pos] = mismatch_base;
                }
            }
            if (strand_id == '0') {
                if (!is_paired_end || flag == '0') {
                    o1 << str << "\n";
                } else {
                    o2 << str << "\n";
                }
            } else {
                if (!is_paired_end || flag == '0') {
                    o1 << PgSAHelpers::reverseComplement(str) << "\n";
                } else {
                    o2 << PgSAHelpers::reverseComplement(str) << "\n";
                }
            }
        }
        pos = 0;
        for (size_t i = 0; i < LQ_read_off_stream.size(); ++i) {
            pos += (read_len - LQ_read_off_stream[i]); /// 注意偏移量的计算
            if (is_paired_end) flag = pe_flag[read_off_stream.size() + i];
            if (!is_paired_end || flag == '0') {
                o1 << lqPgSeq.substr(pos, read_len) << "\n";
            } else {
                o2 << lqPgSeq.substr(pos, read_len) << "\n";
            }
        }
        pos = 0;
        for (size_t i = 0; i < N_read_off_stream.size(); ++i) {
            pos += (read_len - N_read_off_stream[i]); /// 注意偏移量的计算
            if (is_paired_end) flag = pe_flag[read_off_stream.size() + LQ_read_off_stream.size() + i];
            if (!is_paired_end || flag == '0') {
                o1 << nPgSeq.substr(pos, read_len) << "\n";
            } else {
                o2 << nPgSeq.substr(pos, read_len) << "\n";
            }
        }
    } else {
        if (isJoinPgLengthStd) {
            preserve_order_process<uint32_t>(max_mismatch_count, mismatch_decoder_table, hqPgSeq, lqPgSeq, nPgSeq,
                                             mismatch_base_stream, mismatch_count_stream, mismatch_off_stream,
                                             strand_id_stream, read_len, is_paired_end, working_path, o1, o2);
        } else {
            preserve_order_process<uint64_t>(max_mismatch_count, mismatch_decoder_table, hqPgSeq, lqPgSeq, nPgSeq,
                                             mismatch_base_stream, mismatch_count_stream, mismatch_off_stream,
                                             strand_id_stream, read_len, is_paired_end, working_path, o1, o2);
        }
    }
}

void decompress(const Param& param) {
    std::ifstream input(param.f1_path);
    std::ofstream o1, o2;
    uint8_t is_preserve_order;
    uint8_t is_paired_end;
    uint16_t read_len;
    uint8_t block_count;

    input.read(reinterpret_cast<char*>(&is_preserve_order), sizeof(uint8_t));
    input.read(reinterpret_cast<char*>(&is_paired_end), sizeof(uint8_t));
    input.read(reinterpret_cast<char*>(&read_len), sizeof(uint16_t));
    input.read(reinterpret_cast<char*>(&block_count), sizeof(uint8_t));

    if (is_paired_end) {
        o1.open(fs::path(param.working_dir) / (param.output_name + "_1.seq"));
        o2.open(fs::path(param.working_dir) / (param.output_name + "_2.seq"));
    } else {
        o1.open(fs::path(param.working_dir) / (param.output_name + ".seq"));
    }

    for (size_t b_id = 0; b_id < block_count; ++b_id) {
        fs::path working_path = fs::path(param.working_parent_path) / (std::to_string(b_id));
        if (!fs::exists(working_path)) {
            fs::create_directories(working_path);
        } else {
            fs::remove_all(working_path);
            fs::create_directories(working_path);
        }
        decompress_block(input, o1, o2, is_preserve_order, is_paired_end, read_len, working_path);
    }

    int status = std::system(("rm -rf " + param.working_parent_path).c_str());
    if (status != 0) {
        printf("remove tmp directory fail\n");
    }
}
