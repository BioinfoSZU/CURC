## CURC - CUDA Read Sequence Compressor

### Software Requirement
- GCC >= 7.3 (support C++17 standard)
- CMake 3.19
- CUDA Toolkit >= 10.1
- OpenMP

### Platform Requirement
- GPU compute capability >= 6.0

### Operator System
- Linux

### Build
```bash
$ git clone https://github.com/junior-2016/CURC.git && cd CURC 
$ mkdir build && cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER=<nvcc_compiler_path> ..
$ make -j<thread_num>
```

### Command line argument
```bash 
CUDA Read Compressor v1.0.0
Usage:
  curc [OPTION...]

      --working_dir arg        working directory (default: .)
  -c, --compress               compress file
  -d, --decompress             decompress archive
  -i, --input arg              input file path (paired-end fastq paths are 
                               separated by commas)
  -o, --output arg             output file name
      --block_ratio arg        ratio of block size (default: 1)
      --flzma2_level arg       fast-lzma2 compression level [1...10] 
                               (default: 10)
      --flzma2_thread_num arg  fast-lzma2 compression/decompression thread 
                               number (default: 16)
      --preserve_order         preserve order information
  -h, --help                   print usage
```

### Usage
Single-end order non-preserving mode
```bash 
./curc -c -i SRR519094_1.fastq --block_ratio 1 -o SRR519094_1 # compress output is SRR519094_1.curc
./curc -d -i SRR519094_1.curc -o SRR519094_1   # decompress output is SRR519094_1.seq
```
Single-end order preserving mode
```bash 
./curc -c -i SRR519094_1.fastq --block_ratio 1 --preserve_order -o SRR519094_1 # compress output is SRR519094_1.curc
./curc -d -i SRR519094_1.curc -o SRR519094_1   # decompress output is SRR519094_1.seq
```
Paired-end order non-preserving mode
```bash 
./curc -c -i SRR519094_1.fastq,SRR519094_2.fastq --block_ratio 0.5 -o SRR519094 # compress output is SRR519094.curc
./curc -d -i SRR519094.curc -o SRR519094   # decompress output is SRR519094_1.seq and SRR519094_2.seq
```
Paired-end order preserving mode
```bash 
./curc -c -i SRR519094_1.fastq,SRR519094_2.fastq --block_ratio 0.5 --preserve_order -o SRR519094 # compress output is SRR519094.curc
./curc -d -i SRR519094.curc -o SRR519094   # decompress output is SRR519094_1.seq and SRR519094_2.seq
```

### Sample data
[sample.fastq.gz](data/sample.fastq.gz)

### Block ratio
The block ratio is the ratio of the single block size to the size of the entire Fastq file.
The default value of block_ratio is 1, which compresses the entire FASTQ in one block.
If the GPU memory is large enough (eg, Tesla P100-PCI-E 16G), the block size can be set to 50-60 GB. 
If the GPU memory is small, the block size needs to be reduced to avoid out-of-memory exceptions (can be set to about 20 GB).

### Bench result
