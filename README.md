# DNA Variant Calling pipeline
 
![image](https://github.com/Juan-Jeffery/DNA_Variant_Calling_pipeline/blob/main/img/DNA_pipeline.png)

## Package version
Make sure to verify and install the required packages by reviewing the `.yml` file, excluding MutsigCV.
### WES_preprocessing_SNP (QC + Alinment + Variant_Calling(Mutect,varscan,muse))
```bash
Name                    Version                   Build  Channel
bwa                       0.7.17               h5bf99c6_8    bioconda 
picard                    2.18.29                       0    bioconda
samtools                  1.6                  hb116620_7    bioconda 
------------------------------------------------------------------------
fastqc                    0.12.1               hdfd78af_0    bioconda 
multiqc                   1.17                     pypi_0    pypi
trimmomatic               0.39                 hdfd78af_2    bioconda
qualimap                  2.2.2a                        1    bioconda
------------------------------------------------------------------------
gatk4                     4.1.0.0                       0    bioconda 
muse                      1.0.rc               h2e03b76_5    bioconda
varscan                   2.4.6                hdfd78af_0    bioconda
```
### WES_Indel (Variant_Calling(Pindel))
```bash
Name                    Version                   Build  Channel
pindel                    0.2.5b9             h84372a0_10    bioconda
```
### WES_CNV (Variant_Calling(CNVkit))
```bash
Name                    Version                   Build  Channel
cnvkit                    0.9.10             pyhdfd78af_0    bioconda
gistic2                   2.0.23                        0    hcc
```
### WES_Annotation (Annotation(vep, vcf2maf, SigProFiler))
```bash
Name                    Version                   Build  Channel
igv-reports               1.12.0             pyh7cba7a3_0    bioconda
sigmut                    1.0                  hdfd78af_2    bioconda
samtools                  1.10                 h2e538c0_3    bioconda
vcf2maf                   1.6.21               hdfd78af_0    bioconda
ensembl-vep               103.1           pl5262h4a94de4_2    bioconda 
```




