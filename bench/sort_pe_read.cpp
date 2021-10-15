#include <bits/stdc++.h>
#include <cxxopts/cxxopts.hpp>
#include <parallel/algorithm>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
using namespace std;
int main(int argc, char** argv) {
    cxxopts::Options options("sort_read", "sort pair sequence file in pair order");
    options.add_options()
            ("input",  "input file path", cxxopts::value<std::vector<std::string>>())
            ("output", "output file path", cxxopts::value<std::vector<std::string>>())
            ("h,help", "print usage");
    try {
        auto result = options.parse(argc, argv);
        if (result.count("help") || !result.count("input") || !result.count("output")) {
            std::cout << options.help() << "\n";
            std::exit(0);
        }
        auto input_path = result["input"].as<std::vector<std::string>>();
        auto output_path = result["output"].as<std::vector<std::string>>();
        if (input_path.size() != 2 || output_path.size() != 2) {
            std::cout << "use pair sequence input/output path\n";
            std::exit(0);
        }
        if (!std::experimental::filesystem::exists(input_path[0]) ||
            !std::experimental::filesystem::exists(input_path[1])) {
            std::cout << "please check the input path\n";
            std::exit(0);
        }
        std::ifstream in1(input_path[0]), in2(input_path[1]);
        std::ofstream out1(output_path[0]), out2(output_path[1]);
        std::string line1, line2;
        std::vector<std::pair<std::string,std::string>> data;
        while (std::getline(in1, line1)) {
            std::getline(in2, line2);
            data.emplace_back(line1, line2);
        }
        std::vector<uint32_t> id(data.size());
        std::iota(id.begin(), id.end(), 0);
        __gnu_parallel::sort(id.begin(), id.end(), [&] (auto a, auto b) {
            return data[a] < data[b];
        });
        for (auto i : id) {
            out1 << data[i].first << "\n";
            out2 << data[i].second << "\n";
        }
    } catch (cxxopts::OptionException& e) {
        std::cout << e.what() << "\n";
        std::exit(0);
    }
}

