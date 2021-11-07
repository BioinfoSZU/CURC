## CURC - A CUDA-based reference-free read compressor
A GPU-accelerated reference-free compressor for high-throughput sequencing reads of FASTQ files.
### Installation on Linux
CURC currently only supports running on Linux operating systems with Nvidia GPUs.
We have fully tested on the NVIDIA GPU architectures with compute capability >= 3.7 such as Maxwell, Pascal, Volta, Turing 
(check compute capability of your NVIDIA GPU on https://developer.nvidia.com/cuda-gpus#compute). 
If you plan to build and run CURC on old architectures, you can use the `-DCURC_DISABLE_ARCH_CHECK=ON` option for `cmake` that disables checking GPU architecture.
The instructions below will create the executable in the build directory. 

#### Compiler and CUDA requirement
CURC requires GCC version >= 7.3 (support C++17 standard) and CUDA Toolkit version >= 10.1. 
1. Checking GCC version using `gcc --version`. If the GCC version is older than 7.3, use the following command to install GCC7.
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

2. Checking CUDA version using `cat <cuda_path>/version.txt`(eg `cat /usr/local/cuda/version.txt`) or `nvcc --version`.
If CUDA isn't installed or with version older than 10.1, you can download and install CUDA from https://developer.nvidia.com/cuda-toolkit-archive. 

#### Build
CURC uses cmake as the build system, and you can check cmake version using `cmake --version`.
CMake uses `nvcc_path` to detect the CUDA toolkit settings. 
The `nvcc_path` is `<cuda_path>/bin/nvcc` (eg `/usr/local/cuda/bin/nvcc`). 
Make sure that CUDA version >= 10.1 before building, and you can check CUDA version using `cat <cuda_path>/version.txt`). 

Choose one of the following commands to build source based on your system environment.
- with cmake installed and version at least 3.19:
```bash
git clone https://github.com/BioinfoSZU/CURC.git
cd CURC 
mkdir build
cd build
export CC=<gcc_path>  # eg /usr/bin/gcc-7
export CXX=<g++_path> # eg /usr/bin/g++-7 
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER=<nvcc_path> .. 
make
```

- with cmake not installed or cmake version older than 3.19:
```bash
git clone https://github.com/BioinfoSZU/CURC.git
cd CURC
mkdir build
cd build
wget https://github.com/Kitware/CMake/releases/download/v3.21.3/cmake-3.21.3-linux-x86_64.tar.gz
tar -xzf cmake-3.21.3-linux-x86_64.tar.gz
export CC=<gcc_path>  # eg /usr/bin/gcc-7
export CXX=<g++_path> # eg /usr/bin/g++-7 
./cmake-3.21.3-linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER=<nvcc_path> .. 
make
```

### Usage
Run the CURC executable in the build directory with the options below:
```text
CUDA Read Compressor v1.0.0
Usage:
  curc [OPTION...]

      --working_dir arg        working directory (default: .)
  -c, --compress               compress file
  -d, --decompress             decompress archive
  -i, --input arg              input path (paired-end fastq paths are 
                               separated by commas)
  -o, --output arg             output file name (compressed output file use 
                               .curc as extension, decompressed output file 
                               use .seq as extension)
      --block_ratio arg        ratio of block size (default: 1)
      --flzma2_level arg       fast-lzma2 compression level [1...10] 
                               (default: 10)
      --flzma2_thread_num arg  fast-lzma2 compression/decompression thread 
                               number (default: 16)
      --preserve_order         preserve order information
  -h, --help                   print usage
```

### GPU warmup
GPU initialization may be slow on some systems. 
To avoid the impact of GPU initialization for CURC, 
you can execute `nvidia-smi -l 10` in the background to warm up GPU from idle before running CURC.
Another way is to enable persistence mode using `nvidia-smi -i <target gpu> -pm ENABLED`. 

### Block ratio
CURC processes FASTQ file block by block. The block ratio is the single block size divided by the size of the entire FASTQ file.
The default value of block_ratio is 1, which compresses the entire FASTQ in one block.
If the GPU has enough memory (eg, 16 GB), the block size can be set to large enough (eg, 50 GB).
If the GPU memory is small, the block size needs to be reduced to avoid out-of-memory (10-20 GB is reasonable).

### Example Usage
- compress single-end FASTQ reads in order non-preserving mode with one block
```bash
./curc -c -i in.fastq --block_ratio 1 -o single_end_archive # compressed output is single_end_archive.curc
```

- compress single-end FASTQ reads in order preserving mode with one block
```bash
./curc -c -i in.fastq --block_ratio 1 --preserve_order -o single_end_archive # compressed output is single_end_archive.curc
```

- compress paired-end FASTQ reads in order non-preserving mode with two size-equal blocks
```bash
./curc -c -i in_1.fastq,in_2.fastq --block_ratio 0.5 -o paired_end_archive # compressed output is paired_end_archive.curc
```

- compress paired-end FASTQ reads in order preserving mode with two size-equal blocks
```bash
./curc -c -i in_1.fastq,in_2.fastq --block_ratio 0.5 --preserve_order -o paired_end_archive # compressed output is paired_end_archive.curc
```

- decompress single-end compressed archive
```bash
./curc -d -i single_end_archive.curc -o out   # decompressed output is out.seq 
```

- decompress paired-end compressed archive
```bash
./curc -d -i paired_end_archive.curc -o out   # decompressed output is out_1.seq and out_2.seq
```

### Sample data
[SRR554369_1_sample.fastq.gz](data/SRR554369_1_sample.fastq.gz)
