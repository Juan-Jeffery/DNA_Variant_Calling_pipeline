# DNA Variant Calling pipeline
 
![image](https://github.com/Juan-Jeffery/DNA_Variant_Calling_pipeline/blob/main/img/DNA_pipeline.png)

## Package version
---
### WES_preprocessing_SNP (QC + Alinment + Variant_Calling(Mutect,varscan,muse))
'''bash
# Name                    Version                   Build  Channel
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
'''
