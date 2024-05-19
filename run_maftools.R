rm(list=ls())#to remove all objects from the current R environment
setwd("/Users/")
library(maftools)
library(readr)

#將轉完的maf檔轉成tsv檔讀取
tsv <- read_tsv('all_samples_2callers_vep103_maftools.maf')
tsv_L <-read_tsv('all_samples_2callers_vep103.maf')
#合併
tsv <- cbind(tsv,tsv_L$"HGVSp")
names(tsv)[ncol(tsv)] <- "HGVSp"
maf <- read.maf(maf=tsv)

#plotmafSummary
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#oncoplot
oncoplot(maf = maf, top = 20,showTumorSampleBarcodes = T) 
#Transition and Transversions
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = titv)
#lollipopPlot
lollipopPlot( maf = maf,  gene = c('gene_name'),  
              AACol = 'HGVSp_Short',  
              showMutationRate = TRUE,  
              showDomainLabel = FALSE) #　避免重疊
#tmb
tmb(maf, captureSize = 50, logScale = TRUE)
#mutload
mutload = tcgaCompare(maf = maf, cohortName = 'Example-CHOL', logscale = TRUE, capture_size = 50)


###########加入臨床資料###########
clinical = readxl::read_xlsx("clinical.xlsx")
maf=read.maf(maf = tsv, clinicalData = clinical)

#oncoplot
pdf("oncoplot.pdf",20,10)
#One can use any colors, here in this example color palette from RColorBrewer package is used
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del')
print(vc_cols)
#oncoplot(maf = maf_80, colors = vc_cols, top = 20)
oncoplot(maf = tsv, 
         colors = vc_cols, top = 20, sortByAnnotation = T,
         annotationOrder = c("0", "1"), 
         #fontSize = 0.8, 
         legendFontSize = 1.5, anno_height =2,#legend_height = 1,
         annotationFontSize = 1.5, titleFontSize = 1.5) 
#SampleNamefontSize = 1, showTumorSampleBarcodes = TRUE, draw_titv = TRUE)
dev.off()
