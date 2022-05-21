## CURC - A CUDA-based reference-free read compressor
A GPU-accelerated reference-free compressor for high-throughput sequencing reads of FASTQ files.

The current implementation of CURC only supports constant-length reads (<= 512 bases).

### Installation on Linux
CURC is designed to run on Linux operating systems with NVIDIA GPUs.
We have fully tested the program on NVIDIA GPU architectures with compute capability >= 3.7 such as Maxwell, Pascal, Volta, Turing 
(check compute capability of your NVIDIA GPU on https://developer.nvidia.com/cuda-gpus#compute). 
To build and run CURC on old architectures, you can use the `-DCURC_DISABLE_ARCH_CHECK=ON` option for `cmake` that disables checking GPU architecture.
The instructions below can be used to create the executable file in the build directory. 

#### Compiling requirement
CURC should be compiled with GCC version 7.3 or later (supporting C++17 standard) and CUDA Toolkit version 10.1 or later. 
1. To check the GCC version, please use `gcc --version`. If the GCC version is older than 7.3, you can use the following command to install GCC 7.
- On Ubuntu
```bash
sudo add-apt-repository ppa:jonathonf/gcc
sudo apt-get update
sudo apt-get install gcc-7 g++-7
```
- On CentOS
```bash
yum install centos-release-scl
yum install devtoolset-7-gcc-c++
scl enable devtoolset-7 bash # optional step (if you want to set GCC 7 as default compiler in bash)
```
- On anaconda3 without root privileges
```bash
conda create -n gcc7 
conda activate gcc7
conda install gcc_linux-64=7.3.0 gxx_linux-64=7.3.0
cd <anaconda3_root>/envs/gcc7/bin/  
ln -s x86_64-conda_cos6-linux-gnu-gcc gcc 
ln -s x86_64-conda_cos6-linux-gnu-g++ g++ 
export CC=<anaconda3_root>/envs/gcc7/bin/gcc  
export CXX=<anaconda3_root>/envs/gcc7/bin/g++ 
```

2. To check the CUDA version, please use `cat <cuda_path>/version.txt`(eg `cat /usr/local/cuda/version.txt`) or `nvcc --version`.
If CUDA isn't installed or with version older than 10.1, you can download and install CUDA from https://developer.nvidia.com/cuda-toolkit-archive. 

#### Build
CURC uses cmake as the build system, and you can check cmake version using `cmake --version`.
CMake uses `nvcc_path` to detect the CUDA toolkit settings. 
The `nvcc_path` is `<cuda_path>/bin/nvcc` (eg `/usr/local/cuda/bin/nvcc`). 
Make sure that CUDA version >= 10.1 before building, and you can check CUDA version using `cat <cuda_path>/version.txt`). 

Choose one of the following commands to build the source code based on your system environment.
- If cmake with version 3.19 or later has been installed:
```bash
git clone https://github.com/BioinfoSZU/CURC.git
cd CURC 
mkdir build
cd build
export CC=<gcc_path>  # eg export CC=/usr/bin/gcc-7
export CXX=<g++_path> # eg export CXX=/usr/bin/g++-7
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER=<nvcc_path> .. 
make
```

- Otherwise:
```bash
git clone https://github.com/BioinfoSZU/CURC.git
cd CURC
mkdir build
cd build
wget https://github.com/Kitware/CMake/releases/download/v3.21.3/cmake-3.21.3-linux-x86_64.tar.gz
tar -xzf cmake-3.21.3-linux-x86_64.tar.gz
export CC=<gcc_path>  # eg export CC=/usr/bin/gcc-7
export CXX=<g++_path> # eg export CXX=/usr/bin/g++-7
./cmake-3.21.3-linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER=<nvcc_path> .. 
make
```

Note that: if some requirements are not met during the CMake build, and you reinstall the required software, 
please delete `CMakeCache.txt` (if present) from the build directory first and re-run the building command.

### Usage
Run the CURC executable in the build directory with the options below:
```text
CUDA Read Compressor v1.0.0
Usage:
  curc [OPTION...]

      --working_dir arg         working directory (default: .)
  -c, --compress                compress file
  -d, --decompress              decompress archive
  -i, --input arg               input path (paired-end fastq paths are 
                                separated by commas)
  -o, --output arg              output file name (compressed output file 
                                use .curc as extension, decompressed output 
                                file use .seq as extension)
      --block_ratio arg         ratio of block size (default: 1)
      --flzma2_level arg        fast-lzma2 compression level [1...10] 
                                (default: 10)
      --flzma2_thread_num arg   fast-lzma2 compression/decompression thread 
                                number (default: 16)
      --preserve_order          preserve order information
      --decode_buffer_size arg  size of decompress buffer(MB) (default: 
                                1024)
  -v, --version                 print version
  -h, --help                    print usage
```

### GPU warmup and selection
GPU initialization could be slow on some systems. To avoid the impact of GPU initialization for CURC, 
you can execute `nvidia-smi -l 10` in the background to warm up GPU from idle before running CURC.
Another way is to enable persistence mode using `nvidia-smi -i <target gpu> -pm ENABLED`. 

If there are multiple GPUs in the system and some devices are occupied by other compute-intensive tasks, 
you can use `CUDA_VISIBLE_DEVICES` to make only those idle devices visible to CURC. For example:
```bash
CUDA_VISIBLE_DEVICES=2 ./curc <args>      # specify gpu device 2 that CURC uses
```

### Block ratio
CURC processes FASTQ file block by block. The block ratio is defined as the single block size divided by the size of the entire FASTQ file.
The default value of block_ratio is 1, which means compressing the entire FASTQ in one block.
If the GPU has sufficient memory (eg, 16 GB), the block size can be set to a large value (eg, 50 GB), otherwise, 
the block size should be kept small to avoid the out-of-memory issue (around 20 GB is reasonable).

### Example Usage
- To compress single-end FASTQ reads in order non-preserving mode with one block
```bash
./curc -c -i in.fastq --block_ratio 1 -o single_end_archive # compressed output is single_end_archive.curc
```

- To compress single-end FASTQ reads in order preserving mode with one block
```bash
./curc -c -i in.fastq --block_ratio 1 --preserve_order -o single_end_archive # compressed output is single_end_archive.curc
```

- To compress paired-end FASTQ reads in order non-preserving mode with two size-equal blocks
```bash
./curc -c -i in_1.fastq,in_2.fastq --block_ratio 0.5 -o paired_end_archive # compressed output is paired_end_archive.curc
```

- To compress paired-end FASTQ reads in order preserving mode with two size-equal blocks
```bash
./curc -c -i in_1.fastq,in_2.fastq --block_ratio 0.5 --preserve_order -o paired_end_archive # compressed output is paired_end_archive.curc
```

- To decompress single-end compressed archive
```bash
./curc -d -i single_end_archive.curc -o out   # decompressed output is out.seq 
```

- To decompress paired-end compressed archive
```bash
./curc -d -i paired_end_archive.curc -o out   # decompressed output is out_1.seq and out_2.seq
```

### Sample data
|   dataset    |  download command | direct download link   |
|   :----:     |      :----:       |        :----:          |
| SRR635193_1  |`wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR635/SRR635193/SRR635193_1.fastq.gz` | [SRR635193_1.fastq.gz](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR635/SRR635193/SRR635193_1.fastq.gz) |
| SRR445724    |`wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR445/SRR445724/SRR445724.fastq.gz`   | [SRR445724.fastq.gz](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR445/SRR445724/SRR445724.fastq.gz)   |
| SRR065390_1  |`wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065390/SRR065390_1.fastq.gz` | [SRR065390_1.fastq.gz](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065390/SRR065390_1.fastq.gz) |
| SRR065390_2  |`wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065390/SRR065390_2.fastq.gz` | [SRR065390_2.fastq.gz](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065390/SRR065390_2.fastq.gz) |

### Publications
[Xie, Shaohui et al. “CURC: A CUDA-based reference-free read compressor.” Bioinformatics (Oxford, England), btac333. 17 May. 2022, doi:10.1093/bioinformatics/btac333](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btac333/6586792)

### Credits
CURC is based on the architecture of PgRC and also uses parts of PgRC codes in backend encoding.

[Kowalski, S. and Grabowski, S. (2020) PgRC: pseudogenome-based read compressor. Bioinformatics, 36, 2082–2089.](https://academic.oup.com/bioinformatics/article/36/7/2082/5670526)

[Kowalski T, Grabowski SP, Engineering the compression of sequencing reads, BioRxiv2020, https://www.biorxiv.org/content/10.1101/2020.05.01.071720v1.full.](https://www.biorxiv.org/content/10.1101/2020.05.01.071720v1)
