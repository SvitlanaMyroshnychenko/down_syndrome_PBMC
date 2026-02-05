@echo off

echo Trimming SRR11856162...
java -jar E:\Programs\Trimmomatic\trimmomatic-0.40.jar PE -threads 4 ^
data\fastq-raw\SRR11856162_1.fastq.gz ^
data\fastq-raw\SRR11856162_2.fastq.gz ^
data\fastq-trimmed\SRR11856162_1_paired.fastq.gz ^
data\fastq-trimmed\SRR11856162_1_unpaired.fastq.gz ^
data\fastq-trimmed\SRR11856162_2_paired.fastq.gz ^
data\fastq-trimmed\SRR11856162_2_unpaired.fastq.gz ^
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 ^
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

echo Trimming SRR11856164...
java -jar E:\Programs\Trimmomatic\trimmomatic-0.40.jar PE -threads 4 ^
data\fastq-raw\SRR11856164_1.fastq.gz ^
data\fastq-raw\SRR11856164_2.fastq.gz ^
data\fastq-trimmed\SRR11856164_1_paired.fastq.gz ^
data\fastq-trimmed\SRR11856164_1_unpaired.fastq.gz ^
data\fastq-trimmed\SRR11856164_2_paired.fastq.gz ^
data\fastq-trimmed\SRR11856164_2_unpaired.fastq.gz ^
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 ^
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

echo Trimming SRR11856165...
java -jar E:\Programs\Trimmomatic\trimmomatic-0.40.jar PE -threads 4 ^
data\fastq-raw\SRR11856165_1.fastq.gz ^
data\fastq-raw\SRR11856165_2.fastq.gz ^
data\fastq-trimmed\SRR11856165_1_paired.fastq.gz ^
data\fastq-trimmed\SRR11856165_1_unpaired.fastq.gz ^
data\fastq-trimmed\SRR11856165_2_paired.fastq.gz ^
data\fastq-trimmed\SRR11856165_2_unpaired.fastq.gz ^
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 ^
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

echo Trimming SRR11856166...
java -jar E:\Programs\Trimmomatic\trimmomatic-0.40.jar PE -threads 4 ^
data\fastq-raw\SRR11856166_1.fastq.gz ^
data\fastq-raw\SRR11856166_2.fastq.gz ^
data\fastq-trimmed\SRR11856166_1_paired.fastq.gz ^
data\fastq-trimmed\SRR11856166_1_unpaired.fastq.gz ^
data\fastq-trimmed\SRR11856166_2_paired.fastq.gz ^
data\fastq-trimmed\SRR11856166_2_unpaired.fastq.gz ^
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 ^
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

echo Trimming finished
pause
