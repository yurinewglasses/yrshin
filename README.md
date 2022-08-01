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

