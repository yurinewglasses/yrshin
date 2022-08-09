# 2022 Summer URAP - 신유리(생공20)

# 1주차 (7/25~7/29)

- 복령(wolfiporia)의 발달 단계에 따른 차등 발현 유전자(Differential Expressesd Gene) 확인

    <초기 파일>
    - Bok(F1,F2,F3,M1,M2,S1,S2,S3) fastq file
    - wolfiporia_cocos_reference fasta file

## RNA-seq analysis  (HISAT2 → featureCounts → edgeR)

![image](https://user-images.githubusercontent.com/110142232/182016478-331519f0-e5b0-42cf-8039-e3149c53997f.png)


### 1)  HISAT2

- Raw read들을 wolfiporia Reference genome에 mapping
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

- mapping된 read들을 정량화

```r
featureCounts -a wolfiporia.gff -o wolfiporia.out -g ID -T 3 -t gene *sorted.bam
```

### 4) edgeR

- MDS Plot
    - Group별로 전체적인 발현량 패턴에 유의미한 차이가 있는지 확인
    - 전체적으로 서열들의 sequencing이 잘 되었는지 확인
![wolfiporia_mds](https://user-images.githubusercontent.com/110142232/182065930-a71ce901-4f05-47e0-81cb-6e831339ca79.png)


- Data filtering
    - read counts가 너무 적은 유전자는 filtering하여 제외함
    - read counts가 너무 적은 경우, differential expression 분석 시 더욱 민감하게 적용되므로, 신뢰도가 떨어질 수 있음 (read count 10개-50개 차이보다 1000개-5000개 차이가 더 유의미함)


- Normalization
    -  mapping된 read의 개수로 발현량을 정의하기에는 sample별로 시퀀싱 데이터 크기가 다를 수 있고, 유전자나 transcript의 길이에 따라 mapping된 read의 수도 다르기 때문에 객관적인 값으로 보기는 힘들다
    - 이러한 오차를 줄여 객관적인 값을 보여주도록 유전자의 길이, library size를 모두 고려하여 normalization을 진행한다.
    - FPKM, RPKM, TMM(Trimmed mean of M values : edgeR에서 사용)

- BCV plot

![wolfiporia](https://user-images.githubusercontent.com/110142232/182072866-509e24aa-9c9e-41bd-996a-9499608c8eff.png)


- exactTest(Differential Expression)
    
   - BokM-BokF
    
    ![슬라이드1](https://user-images.githubusercontent.com/110142232/182081395-9b408689-b9e7-4162-b589-a33484e220d7.PNG)


    
   - BokF-BokS
        
    ![슬라이드2](https://user-images.githubusercontent.com/110142232/182081425-ff04ced7-cf4c-4e12-9037-9b682f4b9671.PNG)


   - BokM-BokS
       
    ![슬라이드3](https://user-images.githubusercontent.com/110142232/182081454-560c80d8-8f58-46c1-86d4-6b14665ffd9f.PNG)
    
    ![슬라이드4](https://user-images.githubusercontent.com/110142232/182081577-be95c0d4-5cf1-47e0-8feb-ca843fca444b.PNG)

       
       
      

       










## IGV (Intergrative Genomics View)

IGV를 설치하고, 간단한 사용법을 배웠다. wolfiporia cocos reference genome에 대해 align된 read를 시각화했다.

- 통합적인 Genome 데이터셋을 시각화해주는 그래픽 기반 프로그램
- 다양한 포맷의 데이터 로드 가능(array-based, NGS, annotation data)
- 사용한 파일
    - Wolfiporia reference sequence file (복령 참조 염기서열)
    - BokM1 bam file (input file, IGV에서 보고자 하는 sample의 파일)
    - Wolfiporia.gff 파일
    

## FunGAP

Lynux 환경에서 FunGAP을 설치하고, **FunGAP의 사용 방법 및 기능**에 대해 알아보았다. 
- GitHub의 READme
- FunGAP: Fungal Genome Annotation Pipeline using evidence-based gene model evaluation 

FunGAP은 다양한 유전자 예측 프로그램을 통합하고 유전자 예측 모델을 평가하여 진균류 유전체를 해독하는 시스템이다. 

![image](https://user-images.githubusercontent.com/110142232/182183462-77e424dd-322a-4ffd-bb18-49e44e5ce4a2.png)


**1단계: Preprocessing of input data**

- **RepeatMasker** : 조립된 유전체 서열 내 반복 서열을 찾아 masking하여 반복 서열이 유전자 정렬 및 예측을 방해하는 것을 방지
- **Hisat / Trinity** : mRNA 서열분석 reads를 정렬 및 조립하여 contig를 생성




**2단계: Gene prediction**

여러 유전자 예측 프로그램을 사용하여 유전자 모델을 생성하는 단계

- Augustus
- Maker
- Braker




**3단계: Gene model evaluation and filtration**

유전자 모델을 평가하여 최종 유전자 모델을 선별하는 단계이다. 

유전자 모델의 평가는 BLASTp, BUSCO 및 InterProScan으로부터 얻은 점수를 방정식에 따라 계산된 evidence score가 높은 유전자 모델을 최종 유전자 모델로 선별하는 방식이다. 

$$
Evidence score = (BLAST score * coverage) +BUSCO score + ∑Pfam scores
$$

- coverage : 서열 길이 보정값 → blast 결과 출력에서 query coverage
- Pfam scores : InterProScan의 XML 결과에서 hmmer3-match 점수


- BLASTp
    - 주어진 서열 데이터베이스로부터 유사한 서열을 검색하는 프로그램
    - 예측된 유전자 서열이 데이터베이스에 존재하는지를 검색하여 유전자 모델 평가에 사용될 수 있음
    - download_sister_orgs.py 스크립트를 통해 주어진 생물체의 분류군에 대해 NCBI 단백질 서열을 다운로드 할 수 있으며 이는 BLASTp의 데이터베이스로 이용됨
    - 서열의 길이가 길수록 더 높은 bit score를 얻으므로 이를 보정하기 위해 서열 길이 coverage를 곱하여 최종 BLAST evidence 점수를 계산함


- BUSCO
    - 모든 진균 유전체에 보존된 단일 카피 ortholog에 대해 hidden Markov 모델을 데이터베이스 형태로 제공함
    - 예측 유전자가 BUSCO 데이터베이스 모델에 정렬될 경우 다른 그렇지 않은 예측 유전자에 비해 실제 유전자일 가능성이 높다고 판단 가능
    
    
- InterProScan
    - Pfam 도메인 검색
    - 수동으로 큐레이트된 단백질 군의 데이터베이스 제공
    - Pfam 도메인으로 주석 처리된 유전자 모델이 실제 유전자일 가능성이 더 높다는 가정 하에 유전자 모델 평가에 사용
    - Pfam에 대한 점수는 InterProScan의 XML 결과에서 hmmer3-match 점수에 의해 직접 제공됨
    
    
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

![fungap_out](https://user-images.githubusercontent.com/110142232/182061268-d1f34ae7-ca90-4547-9818-83498e0e673f.jpg)

![fungap_out_prot_len_dist](https://user-images.githubusercontent.com/110142232/182065429-bd6e90bc-5dc6-4c80-915f-fd701b342583.png)
![fungap_out_trans_len_dist](https://user-images.githubusercontent.com/110142232/182065446-c58ed160-3a66-4bf8-90bb-b887688ced8c.png)
![SmartSelectImage_2022-08-01-12-14-50](https://user-images.githubusercontent.com/110142232/182065153-d816b928-5d8e-4539-a8fe-50e6ad231747.png)

----------------------------------------------------------------



# 2주차 (8/1~8/5)

목표 : Wild type과 p53 knock-out zebrafish에 1-naphthol을 처리했을 때 유전자 발현량 차이를 비교한다. 

---

group은 다음과 같이 4가지로 나눌 수 있다.

1) Wild type / No 1-naphthol treatment (31~35)

2) Wild type / Treated with 1-naphthol (36~40)

3) p53 knock out / No 1-naphthol treatment (21~25)

4) p53 knock out / Treated with 1-naphthol (26~30)

DEG 분석을 통해 1&2 group, 3&4 group을 비교하여 1-naphthol 처리에 따른 발현량 차이를 분석하는 것을 목표로 삼았다.

SRA 데이터를 제공한 논문과 RNA-seq 데이터 분석 과정의 차이를 비교하며 진행했다. 자세한 논문 내용은 아래 링크를 통해 확인할 수 있다. 

논문 : **TP53 Modulates Oxidative Stress in Gata1+ Erythroid Cells**

[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5312256/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5312256/)


![슬라이드6](https://user-images.githubusercontent.com/110142232/183565957-47016bdc-bc91-4b66-8b0c-d7c589f3180d.PNG)


Reference sequence의 경우 참고 논문에서는 UCSC danRer7(=NCBI Zv9)을 사용했다. 이 서열은 2010년 버전의 것이므로, 더 유의미한 결과를 얻어내기 위해서는 가장 최신화된 서열인 danRer11을 사용할 필요가 있었다. 하지만 분석 과정에서 생겨나는 코딩상의 많은 오류를 해결하기 위해서 우선 참고 논문과 같은 danRer7을 사용한 후, 분석 과정의 신뢰성이 보장된다면(즉, 동일한 데이터를 사용한 논문과 비교적 유사한 결과가 나온다면) 같은 방식으로 danRer11 버전을 사용하기로 방향을 잡았다.


## QC & Pre-processing

Trim Galore을 통해 20개의 sample 데이터(SRR5125621~SRR5125640) 전처리를 수행했다. 

```bash
prefetch SRR5125623 && fastq-dump --split-files SRR5125621
trim_galore --paired -j 4 -o QC SRR5125623_1.fastq SRR5125623_2.fastq
```

## Read alignment

Hisat2 및 samtools를 통해 20개 데이터의 alignment를 진행한다.  사용한 reference genome은 NCBI Zv9이다. max-introlen 값은 zebrafish 관련 유사 연구를 참고하여 설정하고자 했다. 하지만 연구마다 그 범위가 너무 다양했으므로, 10,000/20,000/50,000/100,000/500,000 총 4번의 hisat2를 진행하여 alignment rate가 가장 높게 나타나는 길이를 택했다.

```bash
hisat2 --max-intronlen 50000 -p 3 -x index2 -1 /espeon/analysis1/yrshin/SRA_naphthol/QC/SRR5125623_1_val_1.fq -2 /espeon/analysis1/yrshin/SRA_naphthol/QC/SRR5125623_2_val_2.fq 2> Zv9_SRR5125623hisat2.log | samtools view -@ 3 -bSF4 - | samtools sort -@ 3 - -o Zv9_SRR5125623_sorted.bam
```

## Expression quantification

처음 quantification을 수행할 때는 all.count 데이터를 만들도록 설계된 R pipeline code를 이용하여 진행했으며, 추후 신뢰할만한 결과가 나왔는지 확인하기 위해 FeatureCount를 통해 비교 및 검증하였다.  

- R pipeline code
- FeatureCount

## DE analysis

### MDS Plot

- group간 근접성(prximity)를 시각화하기 위해 MDS plot을 그렸다. 
- 빨간색 표지된 sample이 1-naphthol 처리를 하지 않은 group이고, 하늘색 표지된 sample이 1-naphthol 처리한 group이다. 
- 같은 group끼리 거리가 가까운 경우도 있지만, 몇 개의 sample은 같은 group임에도 동떨어져 있는 모습이 보인다. 
- 개선점 : sample 이름을 너무 길게 설정하여 겹쳐보이는 부분이 많고 정확하게 식별하기 어렵다. 조금 더 간결하게 sample 이름을 설정할 필요가 있다.

![MDS plot](https://user-images.githubusercontent.com/110142232/183566506-0182e8db-83b8-4fba-959c-1c19a96ec111.png)

- grouping이 뚜렷하게 되지 않는 원인이 sample의 문제인지, MA Plot을 그리는 과정에서 코드의 문제인지 파악하는 과정이 필요했다. 데이터를 제공한 논문의 GEO data(FeatureCount 직후)를 다운받아 동일한 코드를 실행해봤다.

![exon_zebrafish_p53_degZv9_p53_mds](https://user-images.githubusercontent.com/110142232/183568211-a1e44c65-19da-490b-86d0-dbfd65924148.png)

- 비교적 뚜렷하게 group이 나뉜 것을 확인했다.


### Heatmap

- Wild type zebrafish에 1-naphthol을 처리했을 때의 발현 차이를 시각화했다.
- Whole genes를 대상으로 그린 것보다 DEG를 대상으로 그린 Heatmap에서 더 뚜렷한 group별 차이가 드러났다.

![슬라이드8](https://user-images.githubusercontent.com/110142232/183565411-a5a4f293-66f6-4577-8dfb-2b17d2a0ba40.PNG)

![슬라이드9](https://user-images.githubusercontent.com/110142232/183565417-9b490993-4090-410d-a152-20d6d2f15b3e.PNG)

![슬라이드10](https://user-images.githubusercontent.com/110142232/183565430-72384f23-5115-49ab-b5c2-ed38c6976fc1.PNG)

![슬라이드11](https://user-images.githubusercontent.com/110142232/183565440-65d5f1ef-689b-4cbd-9a08-1da9ee0e4714.PNG)


뚜렷한 색 차이가 드러나지는 않지만 전반적인 경향성을 보았을 때 1-Naphthol을 처리하지 않은 group과 처리한 group 사이 발현 차이가 드러난다. 특히 No naphthol group에서 up-regulated DEG 범위가 유사하다고 판단했다. 

![exon_zebrafish_p53_degZv9_p53_heatmap_deg](https://user-images.githubusercontent.com/110142232/183571442-aa589109-3a9a-4cc6-b701-70cff1b369b6.png)

검증의 일환으로 논문 GEO expression set data를 이용한 Heatmap(p53)도 그렸다. 


### DEG table

R pipeline을 통해 얻어낸 DEG table이다. 아래 표는 logFC 값이 큰 순서대로 상위 15개의 DEG를 내림차순 정렬한 것이다. 전체 table은 github에 추가로 기재했다. 

1) p53 knock out zebrafish

![SmartSelectImage_2022-08-09-01-03-16](https://user-images.githubusercontent.com/110142232/183565635-dca0fb78-fdef-4618-91fa-05586abf1763.png)

2) Wild type zebrafish

![SmartSelectImage_2022-08-09-01-15-35](https://user-images.githubusercontent.com/110142232/183565592-0ba42a43-ec82-46bd-a405-9a309e2e4495.png)



