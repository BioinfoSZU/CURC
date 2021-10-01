## Benchmark script

### Usage
```bash 
Fastq compression benchmark
Usage:
  benchmark [OPTION...]

      --working_dir arg   working directory (default: .)
      --root_dir arg      root dir of dataset
      --se_db_path arg    single-end database file path
      --pe_db_path arg    paired-end database file path
      --mode arg          se/pe
      --preserve_order    enable preserve order
      --all               compress the entire fastq file
      --program_name arg  spring/pgrc/gpu
      --program_path arg  
      --thread_num arg    thread number in compression
  -h, --help              print usage
```

### FASTQ database file example
- single-end db file example
```text 
SRR554369_1 // each line is a single-end file name
ERR038698_2
....
 
```

- paired-end db file example
```text 
SRR554369 // each line implies two file (SRR554369_1.fastq and SRR554369_2.fastq)
ERR966765
.... 

```

### Example
#### Spring 
- single-end reorder
```bash 
 ./bench --root_dir /run/media/junior/Data/gnome_data_2 --mode se --program_name spring --program_path /home/junior/project/genome_compression/Spring/spring --thread_num 16  --se_db_path ../bench/se_db.txt
```
- single-end preserve order
```bash 
 ./bench --root_dir /run/media/junior/Data/gnome_data_2 --mode se --program_name spring --program_path /home/junior/project/genome_compression/Spring/spring --thread_num 16 --preserve_order --se_db_path ../bench/se_db.txt
```
- paired-end reorder
```bash 
./bench --root_dir /run/media/junior/Data/gnome_data_2 --mode pe --program_name spring --program_path /home/junior/project/genome_compression/Spring/spring --thread_num 16  --pe_db_path ../bench/pe_db.txt 
```
- paired-end preserve order
```bash 
./bench --root_dir /run/media/junior/Data/gnome_data_2 --mode pe --program_name spring --program_path /home/junior/project/genome_compression/Spring/spring --thread_num 16 --preserve_order --pe_db_path ../bench/pe_db.txt
```

#### PgRC
- single-end reorder
```bash 
./bench --root_dir /run/media/junior/Data/gnome_data_2 --mode se --program_name pgrc --program_path /home/junior/project/genome_compression/PgRC/cmake-build-release/PgRC-dev --thread_num 16  --se_db_path ../bench/se_db.txt 
```

- paired-end reorder
```bash 
./bench --root_dir /run/media/junior/Data/gnome_data_2 --mode pe --program_name pgrc --program_path /home/junior/project/genome_compression/PgRC/cmake-build-release/PgRC-dev --thread_num 16  --pe_db_path ../bench/pe_db.txt 
```