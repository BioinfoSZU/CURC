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

#include "decompress.hpp"
#include <bits/stdc++.h>
#include <parallel/algorithm>
#include <thrust/system/omp/execution_policy.h>
#include <experimental/filesystem>
#ifndef ABS
#define ABS std::abs
#endif
#include <range_coder/clr.cdr>
#include <range_coder/simple_model.h>
#include <lzma/VarLenDNACoder/LzmaLib.h>
#include "fast_lzma2_helper.hpp"

namespace fs = std::experimental::filesystem;

template<typename T>
std::vector<T> decompress_to_vector(const fs::path& path) {
    std::ifstream in(path);
    std::vector<T> ret;
    readCompressed(in, ret);
    return ret;
}

std::string decompress_to_string(const fs::path& path) {
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

template<typename T>
void extract_id_to_pos(std::vector<T>& id_to_pos, const fs::path& working_path, bool is_paired_end, uint32_t total_reads_count) {
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

//        uint32_t total_reads_count;
//        in.read(reinterpret_cast<char*>(&total_reads_count), 4);
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

        std::vector<T> right_part(pairsCount);
        for (size_t i = 0; i < pairsCount; ++i) {
            right_part[i] = id_to_pos[pairsCount + i];
            id_to_pos[pairsCount + i] = id_to_pos[i];
        }
        for (size_t i = 0; i < pairsCount; ++i) id_to_pos[i * 2] = id_to_pos[pairsCount + i];
        for (size_t i = 0; i < pairsCount; ++i) id_to_pos[i * 2 + 1] = right_part[i];
    }
}

void extract_pe_pair_order(std::vector<uint32_t>& pe_pair_order, const fs::path& working_path, uint32_t reads_count) {
    std::ifstream in(working_path / "pe_order.comp");
    auto prepare_file = [&](const std::string& filename) {
        uint64_t size;
        in.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
        std::vector<char> buffer(size);
        in.read(buffer.data(), size);
        std::ofstream out(working_path / filename);
        out.write(buffer.data(), size);
    };
    std::vector<uint8_t> offsetPairFlag;
    std::vector<uint8_t> nonOffsetPairFlag;
    std::vector<uint8_t> offsetInUint8Flag;
    std::vector<uint8_t> offsetInUint8Value;
    std::vector<uint8_t> deltaInInt8Flag;
    std::vector<int8_t>  deltaInInt8Value;
    std::vector<uint32_t> offset;
    readCompressed<uint8_t>(in, offsetInUint8Flag);
    readCompressed<uint8_t>(in, offsetInUint8Value);
    readCompressed<uint8_t>(in, deltaInInt8Flag);
    readCompressed<uint8_t>(in, offsetPairFlag);
    readCompressed<uint8_t>(in, nonOffsetPairFlag);
    prepare_file("delta.lzma2");
    lzma2::lzma2_decompress((working_path / "delta.lzma2").c_str(), (working_path / "delta.bin").c_str());
    {
        auto size = fs::file_size(working_path / "delta.bin");
        deltaInInt8Value.resize(size / sizeof(int8_t));
        std::ifstream delta_file(working_path / "delta.bin");
        delta_file.read(reinterpret_cast<char *>(deltaInInt8Value.data()), size);
    }
    prepare_file("offset.lzma2");
    lzma2::lzma2_decompress((working_path / "offset.lzma2").c_str(), (working_path / "offset.bin").c_str());
    {
        auto size = fs::file_size(working_path / "offset.bin");
        offset.resize(size / sizeof(uint32_t));
        std::ifstream offset_file(working_path / "offset.bin");
        offset_file.read(reinterpret_cast<char *>(offset.data()), size);
    }

    pe_pair_order.resize(reads_count);
    std::vector<bool> visit(reads_count, false);
    int64_t pair_count = -1;
    int64_t offsetInUint8ValueIdx = -1;
    int64_t deltaInInt8FlagIdx = -1;
    int64_t deltaInInt8ValueIdx = -1;
    int64_t offsetIdx = -1;
    int64_t ref_prev_offset = 0;
    int64_t prev_offset = 0;
    int64_t pair_index_offset = 0;
    bool match = false;

    for (uint32_t i = 0; i < reads_count; ++i) {
        if (visit[i]) continue;
        if (offsetInUint8Flag[++pair_count]) {
            pair_index_offset = offsetInUint8Value[++offsetInUint8ValueIdx];
        } else if (deltaInInt8Flag[++deltaInInt8FlagIdx]) {
            pair_index_offset = ref_prev_offset + deltaInInt8Value[++deltaInInt8ValueIdx];
            ref_prev_offset = pair_index_offset;
            match = true;
            prev_offset = pair_index_offset;
        } else {
            pair_index_offset = offset[++offsetIdx];
            if (!match || ref_prev_offset != prev_offset) {
                ref_prev_offset = pair_index_offset;
            }
            match = false;
            prev_offset = pair_index_offset;
        }
        pe_pair_order[pair_count * 2] = i;
        pe_pair_order[pair_count * 2 + 1] = i + pair_index_offset;
        visit[i + pair_index_offset] = true;
    }
    int64_t offsetPairFlagIdx = -1;
    int64_t nonOffsetPairFlagIdx = -1;
    for (uint32_t i = 0; i < reads_count / 2; ++i) {
        bool is_swap = offsetInUint8Flag[i] ? offsetPairFlag[++offsetPairFlagIdx] : nonOffsetPairFlag[++nonOffsetPairFlagIdx];
        if (is_swap) {
            std::swap(pe_pair_order[i * 2], pe_pair_order[i * 2 + 1]);
        }
    }
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

std::string recoverRef(const std::string& destRef, std::string& srcRef,
                      std::istream& MapOffSrc, std::istream& MapLenSrc,
                      bool isRefLengthStd, bool srcIsDest) {
    std::string str;
    std::string & ret = srcIsDest ? srcRef : str;
    uint32_t target_match_len;
    PgSAHelpers::readUIntByteFrugal(MapLenSrc, target_match_len);
    uint64_t markPos = 0, posDest = 0;
    while((markPos = destRef.find('%', posDest)) != std::string::npos) {
        ret.append(destRef, posDest, markPos - posDest);
        posDest = markPos + 1;
        uint64_t matchSrcPos = 0;
        if (isRefLengthStd) {
            uint32_t tmp;
            PgSAHelpers::readValue<uint32_t>(MapOffSrc, tmp, false);
            matchSrcPos = tmp;
        } else {
            PgSAHelpers::readValue<uint64_t>(MapOffSrc, matchSrcPos, false);
        }
        uint64_t matchLength = 0;
        PgSAHelpers::readUIntByteFrugal(MapLenSrc, matchLength);
        matchLength += target_match_len;
        ret.append(PgSAHelpers::reverseComplement(srcRef.substr(matchSrcPos, matchLength)));
    }
    ret.append(destRef, posDest, destRef.length() - posDest);
    return ret;
}

template<typename T>
void preserve_order_process(
        char mismatch_decoder_table[128][4],
        const std::string& hqPgSeq, const std::string& lqPgSeq, const std::string& nPgSeq,
        const std::string& mismatch_base_stream,
        const std::vector<uint8_t> & mismatch_count_stream,
        const std::vector<std::vector<uint16_t>> & mismatch_off_stream,
        const std::string& strand_id_stream, const std::vector<T>& id_to_pos,
        uint16_t read_len, bool is_paired_end, const fs::path& working_path,
        std::ofstream & o1, std::ofstream & o2) {
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
    uint32_t reads_count;
    char mismatch_decoder_table[128][4];
    input.read(reinterpret_cast<char*>(&reads_count), sizeof(uint32_t));
    input.read(mismatch_decoder_table['A'], 4);
    input.read(mismatch_decoder_table['T'], 4);
    input.read(mismatch_decoder_table['C'], 4);
    input.read(mismatch_decoder_table['G'], 4);
    uint8_t isJoinRefLengthStd;
    if (is_preserve_order) {
        input.read(reinterpret_cast<char*>(&isJoinRefLengthStd), sizeof(uint8_t));
    }
    uint8_t isRefLengthStd;
    input.read(reinterpret_cast<char*>(&isRefLengthStd), sizeof(uint8_t));
    uint64_t ref_size, unmapping_ref_size, unmapping_N_ref_size;
    input.read(reinterpret_cast<char*>(&ref_size), sizeof(uint64_t));
    input.read(reinterpret_cast<char*>(&unmapping_ref_size), sizeof(uint64_t));
    input.read(reinterpret_cast<char*>(&unmapping_N_ref_size), sizeof(uint64_t));
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
        prepare_file("unmapping_N_read_off.lzma");
        prepare_file("unmapping_read_off.lzma");
        if (is_paired_end) {
            // prepare_file("pe_flag.bsc");
            prepare_file("pe_order.comp");
        }
    }
    prepare_file("ref_off.map.lzma");
    prepare_file("ref_len.map.lzma");
    prepare_file("unmap_ref_off.map.lzma");
    prepare_file("unmap_ref_len.map.lzma");
    prepare_file("unmap_n_ref_off.map.lzma");
    prepare_file("unmap_n_ref_len.map.lzma");
    prepare_file("mismatch_count.lzma");
    prepare_file("mismatch_base.lzma");
    uint8_t max_mismatch_count;
    input.read(reinterpret_cast<char*>(&max_mismatch_count), sizeof(uint8_t));
    for (size_t i = 1; i <= max_mismatch_count; ++i) {
        prepare_file("mismatch_off_" + std::to_string(i) + ".bin");
    }

    std::string strand_id_stream, ref_stream, /*pe_flag,*/ refMapOff, refMapLen, unmapRefMapOff, unmapRefMapLen, unmapNRefMapOff, unmapNRefMapLen, mismatch_base_stream;
    std::vector<uint8_t> mismatch_count_stream;
    std::vector<std::vector<uint16_t>> mismatch_off_stream(max_mismatch_count + 1);
    std::vector<uint16_t> read_off_stream, unmap_N_read_off_stream, unmap_read_off_stream;
    std::vector<uint32_t> id_to_pos_32;
    std::vector<uint64_t> id_to_pos_64;
    std::vector<uint32_t> pe_pair_order;
    std::vector<std::future<void>> futures;
    futures.emplace_back(std::async(std::launch::async, [&] {
        strand_id_stream = bsc_decompress(working_path / "strand_id.bsc");
    }));
    futures.emplace_back(std::async(std::launch::async, [&] {
        ref_stream = extract_reference(working_path, (working_path / "ref.lzma2"), ref_size + unmapping_ref_size + unmapping_N_ref_size);
    }));
    futures.emplace_back(std::async(std::launch::async, [&]  {
        refMapOff = decompress_to_string(working_path / "ref_off.map.lzma");
    }));
    futures.emplace_back(std::async(std::launch::async, [&] {
        refMapLen = decompress_to_string(working_path / "ref_len.map.lzma");
    }));
    futures.emplace_back(std::async(std::launch::async, [&] {
        unmapRefMapOff = decompress_to_string(working_path / "unmap_ref_off.map.lzma");
    }));
    futures.emplace_back(std::async(std::launch::async, [&] {
        unmapRefMapLen = decompress_to_string(working_path / "unmap_ref_len.map.lzma");
    }));
    futures.emplace_back(std::async(std::launch::async, [&] {
        unmapNRefMapOff = decompress_to_string(working_path / "unmap_n_ref_off.map.lzma");
    }));
    futures.emplace_back(std::async(std::launch::async, [&] {
        unmapNRefMapLen = decompress_to_string(working_path / "unmap_n_ref_len.map.lzma");
    }));
    futures.emplace_back(std::async(std::launch::async, [&] {
        mismatch_count_stream = decompress_to_vector<uint8_t>(working_path / "mismatch_count.lzma");
        std::vector<size_t> mismatch_count_stat(max_mismatch_count + 1, 0);
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
    }));
    futures.emplace_back(std::async(std::launch::async, [&] {
        mismatch_base_stream = decompress_to_string(working_path / "mismatch_base.lzma");
    }));
    futures.emplace_back(std::async(std::launch::async, [&] {
        if (!is_preserve_order) {
            read_off_stream = decompress_to_vector<uint16_t>(working_path / "read_off.lzma");
            unmap_N_read_off_stream = decompress_to_vector<uint16_t>(working_path / "unmapping_N_read_off.lzma");
            unmap_read_off_stream = decompress_to_vector<uint16_t>(working_path / "unmapping_read_off.lzma");
        }
    }));
    futures.emplace_back(std::async(std::launch::async, [&] {
        if (!is_preserve_order && is_paired_end) {
            // pe_flag = bsc_decompress(working_path / "pe_flag.bsc");
            extract_pe_pair_order(pe_pair_order, working_path, reads_count);
        }
    }));
    futures.emplace_back(std::async(std::launch::async, [&] {
        if (is_preserve_order) {
            if (isJoinRefLengthStd) {
                extract_id_to_pos(id_to_pos_32, working_path, is_paired_end, reads_count);
            } else {
                extract_id_to_pos(id_to_pos_64, working_path, is_paired_end, reads_count);
            }
        }
    }));
    for (auto& f : futures) f.get();

    std::string unmap_N_ref_seq(ref_stream, ref_size + unmapping_ref_size);
    ref_stream.resize(ref_size + unmapping_ref_size);
    std::string unmap_ref_seq(ref_stream, ref_size);
    ref_stream.resize(ref_size);
    ref_stream.shrink_to_fit();
    std::string ref_seq = std::move(ref_stream);
    std::istringstream refMapOffSrc(refMapOff), refMapLenSrc(refMapLen);
    std::istringstream unmapRefMapOffSrc(unmapRefMapOff), unmapRefMapLenSrc(unmapRefMapLen);
    std::istringstream unmapNRefMapOffSrc (unmapNRefMapOff),  unmapNRefMapLenSrc(unmapNRefMapLen);
    std::string ref_seq_src;
    ref_seq = recoverRef(ref_seq, ref_seq_src, refMapOffSrc, refMapLenSrc, isRefLengthStd, true);
    unmap_ref_seq = recoverRef(unmap_ref_seq, ref_seq, unmapRefMapOffSrc, unmapRefMapLenSrc, isRefLengthStd, false);
    unmap_N_ref_seq = recoverRef(unmap_N_ref_seq,  ref_seq, unmapNRefMapOffSrc, unmapNRefMapLenSrc, isRefLengthStd, false);

    if (!is_preserve_order) {
        if (!is_paired_end) {
            uint64_t pos = 0;
            std::vector<size_t> mismatch_off_cur(max_mismatch_count + 1, 0);
            size_t mismatch_base_cur = 0;
            for (size_t i = 0; i < read_off_stream.size(); ++i) {
                uint8_t mis_cnt = mismatch_count_stream[i];
                uint16_t off = read_off_stream[i];
                pos += off;
                char strand_id = strand_id_stream[i];
                std::string str = ref_seq.substr(pos, read_len);
                if (mis_cnt != 0) {
                    for (size_t k = 0; k < mis_cnt; ++k) {
                        uint16_t mismatch_pos = mismatch_off_stream[mis_cnt][mismatch_off_cur[mis_cnt]++];
                        char origin_base = str[mismatch_pos];
                        char mismatch_base = mismatch_decoder_table[origin_base][mismatch_base_stream[mismatch_base_cur++]];
                        str[mismatch_pos] = mismatch_base;
                    }
                }
                if (strand_id == '0') {
                    o1 << str << "\n";
                } else {
                    o1 << PgSAHelpers::reverseComplement(str) << "\n";
                }
            }
            pos = 0;
            for (size_t i = 0; i < unmap_read_off_stream.size(); ++i) {
                pos += (read_len - unmap_read_off_stream[i]);
                o1 << unmap_ref_seq.substr(pos, read_len) << "\n";
            }
            pos = 0;
            for (size_t i = 0; i < unmap_N_read_off_stream.size(); ++i) {
                pos += (read_len - unmap_N_read_off_stream[i]);
                o1 << unmap_N_ref_seq.substr(pos, read_len) << "\n";
            }
        } else {
            auto inclusive_scan = [read_len] (std::vector<uint64_t>& sum, const std::vector<uint16_t>& value, bool is_complement_off) {
                if (value.empty()) return;
                if (!is_complement_off) {
                    sum[0] = value[0];
                    for (size_t i = 1; i < sum.size(); ++i) sum[i] = sum[i - 1] + value[i];
                } else {
                    sum[0] = read_len - value[0];
                    for (size_t i = 1; i < sum.size(); ++i) sum[i] = sum[i - 1] + (read_len - value[i]);
                }
            };

            std::vector<uint64_t> read_pos(read_off_stream.size());
            inclusive_scan(read_pos, read_off_stream, false);
            read_off_stream.clear();
            read_off_stream.shrink_to_fit();

            std::vector<uint64_t> unmap_read_pos(unmap_read_off_stream.size());
            inclusive_scan(unmap_read_pos, unmap_read_off_stream, true);
            unmap_read_off_stream.clear();
            unmap_read_off_stream.shrink_to_fit();

            std::vector<uint64_t> unmap_N_read_pos(unmap_N_read_off_stream.size());
            inclusive_scan(unmap_N_read_pos, unmap_N_read_off_stream, true);
            unmap_N_read_off_stream.clear();
            unmap_N_read_off_stream.shrink_to_fit();

            std::vector<size_t> mismatch_off_cur(max_mismatch_count + 1, 0);
            size_t mismatch_base_cur = 0;
            std::vector<uint32_t> mismatch_off_start_pos(read_pos.size(), 0);
            std::vector<uint32_t> mismatch_base_start_pos(read_pos.size(), 0);
            for (size_t i = 0; i < read_pos.size(); ++i) {
                if (mismatch_count_stream[i] != 0) {
                    mismatch_off_start_pos[i] = mismatch_off_cur[mismatch_count_stream[i]];
                    mismatch_base_start_pos[i] = mismatch_base_cur;
                    mismatch_off_cur[mismatch_count_stream[i]] += mismatch_count_stream[i];
                    mismatch_base_cur += mismatch_count_stream[i];
                }
            }


            for (size_t idx = 0; idx < reads_count; ++idx) {
                uint32_t i = pe_pair_order[idx];
                if (i < read_pos.size()) {
                    uint8_t mis_cnt = mismatch_count_stream[i];
                    char strand_id = strand_id_stream[i];
                    std::string str = ref_seq.substr(read_pos[i], read_len);
                    if (mis_cnt != 0) {
                        for (size_t k = 0; k < mis_cnt; ++k) {
                            uint16_t mismatch_pos = mismatch_off_stream[mis_cnt][mismatch_off_start_pos[i] + k];
                            char origin_base = str[mismatch_pos];
                            char mismatch_base = mismatch_decoder_table[origin_base][mismatch_base_stream[mismatch_base_start_pos[i] + k]];
                            str[mismatch_pos] = mismatch_base;
                        }
                    }
                    if (strand_id == '0') {
                        if (idx % 2 == 0) {
                            o1 << str << "\n";
                        } else {
                            o2 << str << "\n";
                        }
                    } else {
                        if (idx % 2 == 0) {
                            o1 << PgSAHelpers::reverseComplement(str) << "\n";
                        } else {
                            o2 << PgSAHelpers::reverseComplement(str) << "\n";
                        }
                    }
                } else if (i < (read_pos.size() + unmap_read_pos.size())) {
                    i -= read_pos.size();
                    if (idx % 2 == 0) {
                        o1 << unmap_ref_seq.substr(unmap_read_pos[i], read_len) << "\n";
                    } else {
                        o2 << unmap_ref_seq.substr(unmap_read_pos[i], read_len) << "\n";
                    }
                } else {
                    i -= (read_pos.size() + unmap_read_pos.size());
                    if (idx % 2 == 0) {
                        o1 << unmap_N_ref_seq.substr(unmap_N_read_pos[i], read_len) << "\n";
                    } else {
                        o2 << unmap_N_ref_seq.substr(unmap_N_read_pos[i], read_len) << "\n";
                    }
                }
            }
        }
    } else {
        if (isJoinRefLengthStd) {
            preserve_order_process<uint32_t>(mismatch_decoder_table, ref_seq, unmap_ref_seq, unmap_N_ref_seq,
                                             mismatch_base_stream, mismatch_count_stream, mismatch_off_stream,
                                             strand_id_stream, id_to_pos_32, read_len, is_paired_end, working_path, o1, o2);
        } else {
            preserve_order_process<uint64_t>(mismatch_decoder_table, ref_seq, unmap_ref_seq, unmap_N_ref_seq,
                                             mismatch_base_stream, mismatch_count_stream, mismatch_off_stream,
                                             strand_id_stream, id_to_pos_64, read_len, is_paired_end, working_path, o1, o2);
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
