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
Single-end order non-preserve mode
```bash 
./curc -c -i ERP001775_1.fastq --block_ratio 0.5 -o ERP001775_1 # compress output is ERP001775_1.curc
./curc -d -i ERP001775_1.curc -o ERP001775_1                    # decompress output is ERP001775_1.seq
```
Single-end order preserve mode
```bash 
./curc -c -i ERP001775_1.fastq --block_ratio 0.5 -o ERP001775_1 --preserve_order # compress output is ERP001775_1.curc
./curc -d -i ERP001775_1.curc -o ERP001775_1                                     # decompress output is ERP001775_1.seq
```
Paired-end order non-preserve mode
```bash 
./curc -c -i ERP001775_1.fastq,ERP001775_2.fastq --block_ratio 0.25 -o ERP001775 # compress output is ERP001775.curc
./curc -d -i ERP001775.curc -o ERP001775   # decompress output is ERP001775_1.seq and ERP001775_2.seq
```
Paired-end order preserve mode
```bash 
./curc -c -i ERP001775_1.fastq,ERP001775_2.fastq --block_ratio 0.25 -o ERP001775 --preserve_order # compress output is ERP001775.curc
./curc -d -i ERP001775.curc -o ERP001775   # decompress output is ERP001775_1.seq and ERP001775_2.seq
```

### Block ratio
The block ratio is the ratio of the single block size to the size of the entire Fastq file.
If the GPU memory is large enough (eg, Tesla P100-PCI-E 16G), the block size can be set to 50-60 GB. 
If the GPU memory is small, the block size needs to be reduced to avoid out-of-memory exceptions (can be set to about 20 GB).

### Bench result
