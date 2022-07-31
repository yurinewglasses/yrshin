# 2022 Summer URAP - 신유리(생공20)




## 1주차

### RNA-seq (HISAT2 -> featureCounts -> edgeR)
1) Lynux 환경에서 간단한 명령어 사용법을 학습하고, 필요한 프로그램(HISAT2, FastQC, samtools 등)을 설치했다.
2) HISAT2를 통해 mapping (FASTQ to SAM)


~~~R
hisat2 --max-intronlen 2000 -p 3 -x hello -1 BokF2_val_1.fq.gz -2 BokF2_val_2.fq.gz -S BokF2.sam
~~~
* FASTQ 파일의 Read들이 Reference genome에 mapping 완료되면 각 Read별로 Reference genome에서의 염색체 번호 및 위치가 기록되는데, 이를 Sam 파일이라고 한다.

3) samtools를 이용해 sorting 및 indexing (SAM to BAM)

~~~R
samtools view -b BokF2.sam -o BokF2.bam  #용량을 줄이기 위해 Bam 파일로 변환
samtools sort BokF2.bam -@ 3 -o BokF2_sorted.bam  #염색체와 위치(coordinate)순으로 정렬
~~~

4) featureCounts를 통해 

~~~R
featureCounts -a wolfiporia.gff -o wolfiporia.out -g ID -T 3 -t gene *sorted.bam
~~~

5) edgeR 실행
MDS Plot : 전체적으로 서열들의 sequencing이 잘 되었는지 확인하기 위함, Group 별로 전체적인 발현량 패턴에 차이가 있는지 확인하기 위함

---------------------

### 7/28(목)
1) **IGV**를 설치하고, 간단한 사용법을 배웠다. wolfiporia cocos reference genome에 대해 align된 read를 시각화했다. 

* IGV(Intergrative Genomics View)
  - 통합적인 Genome 데이터셋을 시각화해주는 그래픽 기반 프로그램
  - 다양한 포맷의 데이터 로드 가능(array-based, NGS, annotation data)
  - Wolfiporia reference sequence file (참조 염기서열) / BokM1 bam file (input file, IGV에서 보고자 하는 sample의 파일) /  

3) Lynux 환경에서 conda를 사용하여 **FunGAP을 설치**했다.
4) GitHub의 READme와 참고논문 <FunGAP: Fungal Genome Annotation Pipeline using evidence-based gene model evaluation>을 읽고 **FunGAP의 사용 방법 및 기능**에 대해 알아보았다.



![image](https://user-images.githubusercontent.com/110142232/181458499-53ebe964-d9e3-4ad1-b3a7-096122da1067.png)

#### 1단계
**1) Repeat masking**
: 트랜스포존 반복서열과 같은 genomic regions는 잘못된 alignment를 만들고 gene prediction을 방해하기 때문에 Repeat masking 단계가 필요하다. 

**2) Assembly fo mRNA reads**


### 7/29(금)

~~~python
python fungap.py \
  --output_dir  FUNGAP_wolfiporia \
  --trans_read_1 ../FUNGAP_wolfiporia/BokM_1.fastq \
  --trans_read_2 ../FUNGAP_wolfiporia/BokM_2.fastq \
  --genome_assembly ../FUNGAP_wolfiporia/wolfiporia_filtered.fasta \
  --augustus_species  phanerochaete_chrysosporium \
  --sister_proteome ./sister_orgs/prot_db.faa \
  --num_cores 3 \
  --busco fungi_odb10
  ~~~
