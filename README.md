## CURC - CUDA Read Sequence Compressor

### Software Requirement
- GCC >= 7.3
- CMake 3.19
- CUDA Toolkit >= 10.1

### Platform Requirement
- GPU compute capability >= 6.0

### Operator System
- Linux

### Build
```bash
$ git clone https://github.com/junior-2016/CURC.git && cd CURC 
$ mkdir build && cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_COMPILER=<nvcc_path> ..
$ make -j<thread_num>
```

### Usage
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

### Example
```bash 

```