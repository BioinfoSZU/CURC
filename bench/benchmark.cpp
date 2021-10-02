#include "../third_party/cxxopts/cxxopts.hpp"
#include <bits/stdc++.h>
#include <experimental/filesystem>
using namespace std;
namespace fs = std::experimental::filesystem;

std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

double get_second(const std::string & time) {
    auto tmp = split(time, ':');
    if (tmp.size() == 2) { // m:s
        return std::stod(tmp[0]) * 60 + std::stod(tmp[1]);
    } else {               // h:m:s
        return std::stod(tmp[0]) * 3600 + std::stod(tmp[1]) * 60 + std::stod(tmp[2]);
    }
}

int main(int argc,char ** argv) {
    cxxopts::Options options("benchmark", "Fastq compression benchmark");
    options.add_options()
            ("working_dir", "working directory", cxxopts::value<std::string>()->default_value("."))
            ("root_dir", "root dir of dataset", cxxopts::value<std::string>())
            ("se_db_path", "single-end database file path", cxxopts::value<std::string>())
            ("pe_db_path", "paired-end database file path", cxxopts::value<std::string>())
            ("mode", "se/pe", cxxopts::value<std::string>())
            ("preserve_order", "enable preserve order", cxxopts::value<bool>()->default_value("false"))
            ("all", "compress the entire fastq file", cxxopts::value<bool>()->default_value("false"))
            ("program_name", "spring/pgrc/gpu", cxxopts::value<std::string>())
            ("program_path", "", cxxopts::value<std::string>())
            ("thread_num", "thread number in compression", cxxopts::value<int>())
            ("gpu_ids", "gpu device id list", cxxopts::value<std::vector<int>>())
            ("h,help", "print usage");
    try {
        auto result = options.parse(argc, argv);
        if (result.count("help")          ||
            !result.count("root_dir")     ||
            !result.count("mode")         ||
            !result.count("program_name") ||
            !result.count("program_path") ||
            !result.count("thread_num")) {
            std::cout << options.help() << "\n";
            std::exit(0);
        }
        auto working_dir = result["working_dir"].as<std::string>();
        if (!fs::exists(working_dir)) {
            printf("please check working dir\n");
            std::exit(0);
        }
        auto root_dir = result["root_dir"].as<std::string>();
        if (!fs::exists(root_dir)) {
            printf("please check root dir \n");
            std::exit(0);
        }
        auto mode = result["mode"].as<std::string>();
        if (mode != "se" && mode != "pe") {
            printf("please check mode [se/pe]\n");
            std::exit(0);
        }
        if (mode == "se" && !result.count("se_db_path")) {
            printf("please give se_db_path\n");
            std::exit(0);
        }
        if (mode == "pe" && !result.count("pe_db_path")) {
            printf("please give pe_db_path\n");
            std::exit(0);
        }
        auto program_name = result["program_name"].as<std::string>();
        if (program_name != "spring" && program_name != "pgrc" && program_name != "gpu") {
            printf("please check program name [spring/pgrc/gpu]\n");
            std::exit(0);
        }
        auto program_path = result["program_path"].as<std::string>();
        if (!fs::exists(program_path)) {
            printf("please check program path\n");
            std::exit(0);
        }
        auto enable_compress_all = result["all"].as<bool>();
        auto preserve_order = result["preserve_order"].as<bool>();
        if (enable_compress_all) {
            if (program_name != "spring") {
                printf("only spring support compress entire fastq file\n");
                std::exit(0);
            }
            if (!preserve_order) {
                printf("compress entire fastq only support preserve order mode\n");
                std::exit(0);
            }
        }
        auto thread_num = result["thread_num"].as<int>();
        auto max_thread_num = std::thread::hardware_concurrency();
        if (thread_num == 0 || thread_num > max_thread_num) thread_num = (int) max_thread_num;

        std::vector<int> device_ids;
        if (program_name == "gpu") {
            if (!result.count("gpu_ids")) {
                printf("gpu must provide gpu device id list\n");
                std::exit(0);
            }
            device_ids = result["gpu_ids"].as<std::vector<int>>();
            if (device_ids.empty()) {
                printf("parse device_ids error\n");
                std::exit(0);
            }
        }

        std::string output_file_name = program_name;
        if (enable_compress_all) {
            output_file_name += "_fastq";
        } else {
            output_file_name += "_base";
        }
        output_file_name += ("_" + mode);
        if (preserve_order) {
            output_file_name += "_order";
        } else {
            output_file_name += "_reorder";
        }
        fs::path working_path = output_file_name;
        if (!fs::exists(working_path)) {
            fs::create_directories(working_path);
        }
        std::ofstream csv(output_file_name + ".csv");
        csv << "dataset archive_size(byte) comp_time(s) comp_memory(kb) comp_cpu_usage decomp_time(s) decomp_memory(kb) decomp_cpu_usage \n";

        std::string comp_time_cmd = "/usr/bin/time -f '%P %M %E' -o comp_time.log ";
        std::string decomp_time_cmd = "/usr/bin/time -f '%P %M %E' -o decomp_time.log ";
        std::ifstream db_file, db_path;
        std::string f;
        if(mode == "se") {
            db_file.open(result["se_db_path"].as<std::string>());
        } else {
            db_file.open(result["pe_db_path"].as<std::string>());
        }

        int status;
        while (std::getline(db_file, f)) {
            fs::path comp_output_path = working_path / (f + "_comp.log");
            fs::path decomp_output_path = working_path / (f + "_decomp.log");
            fs::path f1_path, f2_path;
            if (mode == "se") {
                f1_path = fs::path(root_dir) / (f + ".fastq");
            } else {
                f1_path = fs::path(root_dir) / (f + "_1.fastq");
                f2_path = fs::path(root_dir) / (f + "_2.fastq");
            }
            double block_ratio = 0;
            if (program_name == "gpu") {
                std::string block_ratio_str;
                std::getline(db_file, block_ratio_str);
                block_ratio = std::stod(block_ratio_str);
            }
            std::string comp_cmd;
            if (program_name == "gpu") {
                comp_cmd = "CUDA_VISIBLE_DEVICES=" + std::to_string(device_ids[0]);
                for (size_t i = 1; i < device_ids.size(); ++i) {
                    comp_cmd += "," + std::to_string(device_ids[i]);
                }
                comp_cmd += " " + comp_time_cmd;
            } else {
                comp_cmd = comp_time_cmd;
            }
            std::string decomp_cmd = decomp_time_cmd;
            size_t archive_size;

            if (program_name == "spring") {
                if (mode == "se") {
                    decomp_cmd += program_path + " -d -t " + std::to_string(thread_num) + " -i archive.spring -o archive.fastq";
                } else {
                    decomp_cmd += program_path + " -d -t " + std::to_string(thread_num) + " -i archive.spring -o archive_1.fastq archive_2.fastq";
                }
                if (enable_compress_all) { // only preserve order
                    comp_cmd += program_path + " -c -t " + std::to_string(thread_num) + " -i " + f1_path.c_str() + " ";
                    if (mode == "pe") comp_cmd += f2_path.c_str() ;
                    comp_cmd += " -o archive.spring";
                } else {
                    if (preserve_order) {
                        comp_cmd += program_path + " -c --no-ids --no-quality -t " + std::to_string(thread_num) + " -i " + f1_path.c_str() + " ";
                        if (mode == "pe") comp_cmd += f2_path.c_str() ;
                        comp_cmd += " -o archive.spring";
                    } else {
                        comp_cmd += program_path + " -c --no-ids --no-quality -r -t " + std::to_string(thread_num) + " -i " + f1_path.c_str() + " ";
                        if (mode == "pe") comp_cmd += f2_path.c_str() ;
                        comp_cmd += " -o archive.spring";
                    }
                }
                comp_cmd += " > " + comp_output_path.string() ;
                decomp_cmd += " > " + decomp_output_path.string();

                printf("%s\n", comp_cmd.c_str());
                status = std::system(comp_cmd.c_str());
                if (status != 0) {
                    printf("spring compress error\n");
                    std::exit(0);
                }
                archive_size = fs::file_size("archive.spring");

                printf("%s\n", decomp_cmd.c_str());
                status = std::system(decomp_cmd.c_str());
                if (status != 0) {
                    printf("spring decompress error\n");
                    std::exit(0);
                }

                std::system("rm archive.spring");
                if (mode == "se") {
                    std::system("rm archive.fastq");
                } else {
                    std::system("rm archive_1.fastq && rm archive_2.fastq");
                }
            } else if (program_name == "pgrc") {
                decomp_cmd += program_path + " -t " + std::to_string(thread_num) + " -d archive.pgrc > " + decomp_output_path.c_str();
                comp_cmd += program_path + " -c 3 -t " + std::to_string(thread_num) ;
                if (preserve_order) comp_cmd += " -o ";
                comp_cmd += " -i " + f1_path.string() + " ";
                if (mode == "pe") comp_cmd += f2_path.string() + " ";
                comp_cmd += "archive.pgrc > " + comp_output_path.string();

                printf("%s\n", comp_cmd.c_str());
                status = std::system(comp_cmd.c_str());
                if (status != 0) {
                    printf("pgrc compress error\n");
                    std::exit(0);
                }
                archive_size = fs::file_size("archive.pgrc");

                printf("%s\n", decomp_cmd.c_str());
                status = std::system(decomp_cmd.c_str());
                if (status != 0) {
                    printf("pgrc decompress error\n");
                    std::exit(0);
                }

                std::system("rm archive.pgrc");
                if (mode == "se") {
                    std::system("rm archive.pgrc_out");
                } else {
                    std::system("rm archive.pgrc_out_1 && rm archive.pgrc_out_2");
                }
            } else if (program_name == "gpu") {
                decomp_cmd += program_path + " -d -i archive.curc -o archive > " + decomp_output_path.c_str();
                comp_cmd += program_path + " -c -i " + f1_path.string();
                if (mode == "pe") comp_cmd += "," + f2_path.string() + " ";
                if (preserve_order) comp_cmd += " --preserve_order ";
                if (block_ratio == 0) {
                    printf("block ratio is empty\n");
                    std::exit(0);
                }
                comp_cmd += " --block_ratio " + std::to_string(block_ratio);
                comp_cmd += " -o archive > " + comp_output_path.string();

                printf("%s\n", comp_cmd.c_str());
                status = std::system(comp_cmd.c_str());
                if (status != 0) {
                    printf("gpu compress error\n");
                    // std::exit(0);
                }
                archive_size = fs::file_size("archive.curc");

                printf("%s\n", decomp_cmd.c_str());
                status = std::system(decomp_cmd.c_str());
                if (status != 0) {
                    printf("gpu decompress error\n");
                    // std::exit(0);
                }

                std::system("rm archive.curc");
                if (mode == "se") {
                    std::system("rm archive.seq");
                } else {
                    std::system("rm archive_1.seq && rm archive_2.seq");
                }
            }
            std::ifstream comp_time_log("comp_time.log");
            std::ifstream decomp_time_log("decomp_time.log");
            std::string comp_cpu_usage, comp_memory, comp_time;
            std::string decomp_cpu_usage, decomp_memory, decomp_time;
            comp_time_log >> comp_cpu_usage >> comp_memory >> comp_time;
            decomp_time_log >> decomp_cpu_usage >> decomp_memory >> decomp_time;

            csv << f << " " << archive_size << " "
                << get_second(comp_time) << " " << comp_memory << " " << comp_cpu_usage << " "
                << get_second(decomp_time) << " " << decomp_memory << " " << decomp_cpu_usage << " \n" << std::flush;
        }

    } catch (cxxopts::OptionException& e) {
        std::cout << e.what() << "\n";
        std::exit(0);
    }
}

