##2020.02.21
## Okamoto wheat mRNA-seq with udonko
## KALLISTO-SlEUTH pairwise DEG detection & PCA,VC plots

library(sleuth)
s2c <- read.table("meta(udonko).txt", header=TRUE, stringsAsFactors=FALSE)
head(s2c)
#sample infec geno                         path            run
#1 udon_01    no null 2_kallisto.pe/01_Non_Nu1_S1_ 01_Non_Nu1_S1_
#2 udon_02    no null 2_kallisto.pe/02_Non_Nu2_S2_ 02_Non_Nu2_S2_
#3 udon_03    no null 2_kallisto.pe/03_Non_Nu3_S3_ 03_Non_Nu3_S3_
#4 udon_04    no   ox 2_kallisto.pe/04_Non_Ox1_S4_ 04_Non_Ox1_S4_
#5 udon_05    no   ox 2_kallisto.pe/05_Non_Ox2_S5_ 05_Non_Ox2_S5_
#6 udon_06    no   ox 2_kallisto.pe/06_Non_Ox3_S6_ 06_Non_Ox3_S6_

##### PCA whole samples
so <- sleuth_prep(s2c)
plot_pca(so, text_labels = T, color_by = "geno")

##### 1.Whole-sample analysis (udonko-effect upon PYL-OX; LRT)
so <- sleuth_fit(so, ~infec + geno + infec:geno, 'full')
so <- sleuth_fit(so, ~infec + geno, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
results_table_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
table(results_table_lrt[,"qval"] < 0.05)
#FALSE  TRUE 
#27465 30662 #from PE data ...Too many

##Anyway, make the output table
#collect DEG ID
sig_ids <- results_table_lrt$target_id[which(results_table_lrt$qval < 0.05)]
sig_table <- results_table_lrt[results_table_lrt$target_id %in% sig_ids,][1:3]
head(sig_table)
#target_id         pval         qval
#1 TRAESCS6D02G065100.1 4.089755e-19 2.377252e-14
#2 TRAESCS7A02G339900.1 2.144114e-16 4.154364e-12
#3 TRAESCS4B02G055300.1 1.522033e-16 4.154364e-12
#4 TRAESCS5D02G217200.1 8.763454e-16 1.273483e-11
#5 TRAESCS2B02G516700.1 2.100946e-15 2.442434e-11

#Extract TPM data
tpm_matrix <- sleuth_to_matrix(so, 'obs_raw', 'tpm')
d <- tpm_matrix #for sure
d <- cbind(rownames(d), data.frame(d, row.names=NULL))
colnames(d)[1] <- "target_id"
tpm <- d[d$target_id %in% sig_ids,]
dim(tpm)
#[1] 30662    13

#Merging & outputs
out <- merge(tpm, sig_table, by="target_id", sort=T)
dim(out)
#[1] 30662    15
write.csv(out, file="0_DEG_mixed.effect.csv", row.names=FALSE)


##### 2.Do pairwise anlaysis (LRT+Wald) repeat x4
smpl <- s2c[c(1:6),]
smpl
#sample infec geno                         path            run
#1 udon_01    no null 2_kallisto.pe/01_Non_Nu1_S1_ 01_Non_Nu1_S1_
#2 udon_02    no null 2_kallisto.pe/02_Non_Nu2_S2_ 02_Non_Nu2_S2_
#3 udon_03    no null 2_kallisto.pe/03_Non_Nu3_S3_ 03_Non_Nu3_S3_
#4 udon_04    no   ox 2_kallisto.pe/04_Non_Ox1_S4_ 04_Non_Ox1_S4_
#5 udon_05    no   ox 2_kallisto.pe/05_Non_Ox2_S5_ 05_Non_Ox2_S5_
#6 udon_06    no   ox 2_kallisto.pe/06_Non_Ox3_S6_ 06_Non_Ox3_S6_

## Generate sleuth object with liklihood ratio test (LRT)
so <- sleuth_prep(smpl, ~ geno)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
results_table_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
table(results_table_lrt[,"qval"] < 0.05)
#FALSE  TRUE 
#33496 24100  #1
#14885 44890  #2
#14779 44213 #3
#26890 32168 #4


#Optional, output PCA plot to check the data repeatablity
plot_pca(so, color_by = "geno")

########### Do the test again with Wald test (WT)
models(so)
#[  full  ]
#formula:  ~geno 
#data modeled:  obs_counts 
#transform sync'ed:  TRUE 
#coefficients:
#(Intercept)
#genoox
#[  reduced  ]
#formula:  ~1 
#data modeled:  obs_counts 
#transform sync'ed:  TRUE 
#coefficients:
#  (Intercept)

so <- sleuth_wt(so, "genoox")
results_table_wt <- sleuth_results(so, 'genoox')

################ Prepare output data
#Extract common DEG from two LRT & WT analyses
d5.lrt.sig_ids <- results_table_lrt$target_id[which(results_table_lrt$qval < 0.05)]
d5.wt.sig_ids <- results_table_wt$target_id[which(results_table_wt$qval < 0.05)]
shared_ids <- d5.wt.sig_ids[d5.wt.sig_ids %in% d5.lrt.sig_ids]
shared_results <- results_table_wt[results_table_wt$target_id %in% shared_ids,]
out <- shared_results[,c(1:4)]
out$log2FC <- log(exp(out$b),2)   #<<- Add Log2FC column
out <- subset(out, select=-c(b))  #<<- Remove b column
out <- subset(out, out$log2FC >= 1 | out$log2FC <= -1) #<<- select by FC
head(out)
#             target_id pval qval     log2FC
#1 TRAESCS2B02G516700.1    0    0  -5.045720
#2 TRAESCS3B02G181300.1    0    0  -7.813321
#3 TRAESCS3B02G002500.1    0    0  -5.937973
#4 TRAESCS4D02G357400.1    0    0 -13.302269
#5  TRAESCSU02G263800.1    0    0   3.930788
#6 TRAESCS2D02G033700.1    0    0 -13.262597
shared_ids <- as.vector(out$target_id)  #<<-prepare new id list

#extract raw TPM data
tpm_matrix <- sleuth_to_matrix(so, 'obs_raw', 'tpm')
d <- tpm_matrix #for sure
d <- cbind(rownames(d), data.frame(d, row.names=NULL))
colnames(d)[1] <- "target_id"
tpm <- d[d$target_id %in% shared_ids,]
head(tpm)
#              target_id    udon_01   udon_02   udon_03     udon_04     udon_05
#1      synTaPYL1a_GFP.1   0.192933  0.000000  0.000000 599.3718696 540.0187906
#3  TRAESCS1A02G000200.1  21.082491  9.056873  8.929441   5.5781065   3.0852960
#4  TRAESCS1A02G000300.1   4.163074  1.661490  1.369813   0.8915954   0.6285683
#5  TRAESCS1A02G000400.1 105.590076 45.115719 27.451366  17.5930123  13.7870934
#10 TRAESCS1A02G000900.1  44.528381 22.782112 14.971474  10.4312024   7.1084031
#11 TRAESCS1A02G001000.1   3.593321  1.727277  1.717704   1.0074722   0.7047051
#udon_06
#1  546.1632851
#3    3.1503855
#4    0.8455380
#5   13.2998150
#10   6.5608455

#merge two tables 
tmp <- merge(tpm, out, by="target_id", sort=T)
head(tmp)
#             target_id    udon_01   udon_02   udon_03     udon_04     udon_05
#1     synTaPYL1a_GFP.1   0.192933  0.000000  0.000000 599.3718696 540.0187906
#2 TRAESCS1A02G000200.1  21.082491  9.056873  8.929441   5.5781065   3.0852960
#3 TRAESCS1A02G000300.1   4.163074  1.661490  1.369813   0.8915954   0.6285683
#4 TRAESCS1A02G000400.1 105.590076 45.115719 27.451366  17.5930123  13.7870934
#5 TRAESCS1A02G000900.1  44.528381 22.782112 14.971474  10.4312024   7.1084031
#6 TRAESCS1A02G001000.1   3.593321  1.727277  1.717704   1.0074722   0.7047051
#udon_06         pval         qval    log2FC
#1 546.1632851 6.833632e-26 2.540929e-24 14.154461
#2   3.1503855 4.641212e-04 2.012310e-03 -1.616758
#3   0.8455380 4.196066e-03 1.404933e-02 -1.410684
#4  13.2998150 2.081919e-03 7.604221e-03 -1.723979
#5   6.5608455 6.892059e-04 2.860113e-03 -1.598406
#6   0.8832189 3.892803e-04 1.719137e-03 -1.304420

#preprare "whole" output (just repeat the steps)
fdr <- results_table_wt[,c(1:4)]
fdr$log2FC <- -log(exp(fdr$b),2)   #<<- Add Log2FC column
fdr <- subset(fdr, select=-c(b))   #<<- Remove b column
all <- merge(d, fdr, by="target_id", sort=T)
head(all)

#write-down
write.csv(tmp, file="sleuth-DEG.csv", row.names=FALSE)
write.csv(all, file="sleuth-whole.csv", row.names=FALSE) #optional

##############[Optional] Volcano-plot
#collect the data & roundup numbers
vc <- all[,c("target_id","qval","log2FC")]
vc$qval <- -log(vc$qval, 10)
vc$qval <- round(vc$qval, digits=1)
vc$log2FC <- round(vc$log2FC, digits=1)
head(vc)
#            target_id qval log2FC
#1     synTaPYL1a_GFP.1 23.6  -14.2
#2 TRAESCS1A02G000100.1   NA     NA
#3 TRAESCS1A02G000200.1  2.7    1.6
#4 TRAESCS1A02G000300.1  1.9    1.4
#5 TRAESCS1A02G000400.1  2.1    1.7
#6 TRAESCS1A02G000500.1   NA     NA

#Generate & cound DEGs
v1 <- subset(vc, vc$log2FC >= 1 & vc$qval > sqrt(2)) ## UP-DEG
v2 <- subset(vc, vc$log2FC <= -1 & vc$qval > sqrt(2)) ## DN-DEG
deg_ids <- c(as.vector(v1$target_id), as.vector(v2$target_id))
`%notin%` <- Negate(`%in%`) #temporary designates "not-in" func
v3 <- vc[vc$target_id %notin% deg_ids,] ## non-DEG

# Counting DEGs
up.deg <- length(v1[,1])
dn.deg <- length(v2[,1])

# Reduce spots
v1 <- unique(v1)
v2 <- unique(v2)
v3 <- unique(v3)

# Drawing plot with solid limits
plot(v3$log2FC,v3$qval, pch=20, cex=0.4, col="grey70", ylim=c(0,350), xlim=c(-15,15), xlab="Log2FC", ylab="-log10FDR")
points(v1$log2FC,v1$qval, pch=20, cex=0.5, col="#ED0422")
points(v2$log2FC,v2$qval, pch=20, cex=0.5, col="#3498DB")
abline(v=0, lty=2)

# Add numbers of DEG on the plot
text(10,340, up.deg, pos=4, cex=1.4, col="#ED0422")
text(10,310, dn.deg, pos=4, cex=1.4, col="#3498DB")
