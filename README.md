# 2022 Summer URAP - 신유리(생공20)

# 1주차 (7/25~7/29)

## RNA-seq  (HISAT2 → featureCounts → edgeR)

### 1)  HISAT2

- FASTQ to SAM
- FASTQ 파일의 Read들이 Reference genome에 mapping 완료되면 각 Read별로 Reference genome에서의 염색체 번호 및 위치가 기록되는데, 이를 Sam 파일이라고 함.

```r
hisat2 --max-intronlen 2000 -p 3 -x hello -1 BokF2_val_1.fq.gz -2 BokF2_val_2.fq.gz -S BokF2.sam
```

### 2) samtools

- SAM to BAM
- sorting 및 indexing

```r
samtools view -b BokF2.sam -o BokF2.bam  #용량을 줄이기 위해 Bam 파일로 변환
samtools sort BokF2.bam -@ 3 -o BokF2_sorted.bam  #염색체와 위치(coordinate)순으로 정렬
```

### 3) featureCounts

- 

```r
featureCounts -a wolfiporia.gff -o wolfiporia.out -g ID -T 3 -t gene *sorted.bam
```

### 4) edgeR

- MDS Plot
    - Group별로 전체적인 발현량 패턴에 유의미한 차이가 있는지 확인
    - 전체적으로 서열들의 sequencing이 잘 되었는지 확인
- Nomalization

## IGV (Intergrative Genomics View)

IGV를 설치하고, 간단한 사용법을 배웠다. wolfiporia cocos reference genome에 대해 align된 read를 시각화했다.

- 통합적인 Genome 데이터셋을 시각화해주는 그래픽 기반 프로그램
- 다양한 포맷의 데이터 로드 가능(array-based, NGS, annotation data)
- 사용한 파일
    - Wolfiporia reference sequence file (복령 참조 염기서열)
    - BokM1 bam file (input file, IGV에서 보고자 하는 sample의 파일)
    - Wolfiporia.gff 파일
    

## FunGAP

Lynux 환경에서 conda를 사용하여 FunGAP을 설치했다. GitHub의 READme 및 <FunGAP: Fungal Genome Annotation Pipeline using evidence-based gene model evaluation> 가이드를 읽고 **FunGAP의 사용 방법 및 기능**에 대해 알아보았다. 

![Untitled](2022%20Summer%20URAP%20-%20%E1%84%89%E1%85%B5%E1%86%AB%E1%84%8B%E1%85%B2%E1%84%85%E1%85%B5(%E1%84%89%E1%85%A2%E1%86%BC%E1%84%80%E1%85%A9%E1%86%BC20)%20baf0f04a85934cdc948cbb4b2cda662e/Untitled.png)

<실행 코드>

```python
python fungap.py \
  --output_dir  FUNGAP_wolfiporia \
  --trans_read_1 ../FUNGAP_wolfiporia/BokM_1.fastq \
  --trans_read_2 ../FUNGAP_wolfiporia/BokM_2.fastq \
  --genome_assembly ../FUNGAP_wolfiporia/wolfiporia_filtered.fasta \
  --augustus_species  phanerochaete_chrysosporium \
  --sister_proteome ./sister_orgs/prot_db.faa \
  --num_cores 3 \
  --busco fungi_odb10
```

## Rosalind
