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

#include <bits/stdc++.h>
#include <experimental/filesystem>
#include "Param.hpp"
#include "constant.hpp"
#include "preprocess.hpp"
#include "decompress.hpp"
#include "../version.h"
#include "../third_party/cxxopts/cxxopts.hpp"
namespace fs = std::experimental::filesystem;

int main(int argc, char** argv) {
    cxxopts::Options options("curc", "CUDA Read Compressor " CURC_VERSION_STRING);
    options.add_options()
            ("working_dir", "working directory", cxxopts::value<std::string>()->default_value("."))
            ("c,compress", "compress file", cxxopts::value<bool>())
            ("d,decompress", "decompress archive", cxxopts::value<bool>())
            ("i,input", "input path (paired-end fastq paths are separated by commas)", cxxopts::value<std::vector<std::string>>())
            ("o,output", "output file name (compressed output file use .curc as extension, decompressed output file use .seq as extension)", cxxopts::value<std::string>())
            ("block_ratio", "ratio of block size", cxxopts::value<double>()->default_value("1"))
            ("flzma2_level", "fast-lzma2 compression level [1...10]", cxxopts::value<int>()->default_value("10"))
            ("flzma2_thread_num", "fast-lzma2 compression/decompression thread number", cxxopts::value<int>()->default_value("16"))
            ("preserve_order", "preserve order information", cxxopts::value<bool>()->default_value("false"))
            ("decode_buffer_size", "size of decompress buffer(MB)", cxxopts::value<uint64_t>()->default_value("1024"))
            ("v,version", "print version")
            ("h,help", "print usage");
    try {
        auto result = options.parse(argc, argv);
        if (result.count("version")) {
            std::cout << CURC_VERSION_STRING << "\n";
            std::exit(0);
        }
        if (result.count("help")    ||
            !result.count("input")  ||
            !result.count("output") ||
            (!result.count("compress") && !result.count("decompress"))) {
            std::cout << options.help() << "\n";
            std::exit(0);
        }
        if (result.count("compress") && result.count("decompress")) {
            std::cout << "compression and decompression options conflict\n";
            std::exit(0);
        }

        Param param;
        param.working_dir = result["working_dir"].as<std::string>();
        if (!fs::exists(param.working_dir)) {
            std::cout << "invalid working dir\n";
            std::exit(0);
        }

        auto working_parent_path_tmp = fs::path(param.working_dir) / "XXXXXX";
        param.working_parent_path = mkdtemp(&working_parent_path_tmp.string()[0]);
        if (!fs::exists(param.working_parent_path)) {
            fs::create_directories(param.working_parent_path);
        } else {
            fs::remove_all(param.working_parent_path);
            fs::create_directories(param.working_parent_path);
        }

        {
            std::string alphabet("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
            std::random_device rd;
            std::mt19937_64 generator(rd());
            std::shuffle(alphabet.begin(), alphabet.end(), generator);
            param.random_block_prefix_id = alphabet.substr(0, 10);
        }

        param.output_name = result["output"].as<std::string>();
        param.block_ratio = result["block_ratio"].as<double>();
        if (param.block_ratio <= 0 || param.block_ratio > 1) {
            std::cout << "block ratio is a float number in the range (0,1]\n";
            std::exit(0);
        }
        param.flzma2_level = result["flzma2_level"].as<int>();
        param.flzma2_thread_num = result["flzma2_thread_num"].as<int>();
        param.is_preserve_order = result["preserve_order"].as<bool>();
        if (result.count("compress")) {
            auto input = result["input"].as<std::vector<std::string>>();
            if (input.empty() || input.size() > 2) {
                std::cout << "invalid input path number\n";
                std::exit(0);
            }
            if (!std::experimental::filesystem::exists(input[0]) ||
               (fs::path(input[0]).extension() != ".fastq" && fs::path(input[0]).extension() != ".fq")) {
                std::cout << "please check the first fastq path (filename extension is .fastq or .fq)\n";
                std::exit(0);
            }
            param.f1_path = input[0];
            if (input.size() == 2) {
                if (!std::experimental::filesystem::exists(input[1]) ||
                   (fs::path(input[1]).extension() != ".fastq" && fs::path(input[1]).extension() != ".fq")) {
                    std::cout << "please check the second fastq path (filename extension is .fastq or .fq)\n";
                    std::exit(0);
                }
                param.f2_path = input[1];
            }
            param.is_paired_end = !param.f2_path.empty();
            compress(param);
        } else if (result.count("decompress")) {
            auto input = result["input"].as<std::vector<std::string>>();
            if (input.size() != 1) {
                std::cout << "invalid input path number\n";
                std::exit(0);
            }
            if (!std::experimental::filesystem::exists(input[0]) ||
                (fs::path(input[0]).extension() != archive_name_suffix)) {
                std::cout << "please check the archive path (filename extension is .curc)\n";
                std::exit(0);
            }
            param.f1_path = input[0];
            param.decode_buffer_size = result["decode_buffer_size"].as<uint64_t>();
            decompress(param);
        }
    } catch (cxxopts::OptionException& e) {
        std::cout << e.what() << "\n";
        std::exit(0);
    }
    return 0;
}
