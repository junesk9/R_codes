##2020.02.07
## ALRC wheat mRNA-seq
## KALLISTO-SlEUTH pairwise DEG detection & PCA,VC plots

library(sleuth)
s2c <- read.table("meta(se).txt", header=TRUE, stringsAsFactors=FALSE)

#### Select rows to be analysis
smpl <- s2c[c(2:6),]
smpl
#       sample cond Time                     path        run
#      sample cond Time                     path        run
#2 matsu-ta02 cool    0 3_kallisto.se/matsu-ta02 matsu-ta02
#3 matsu-ta03 cool    0 3_kallisto.se/matsu-ta03 matsu-ta03
#4 matsu-ta04 cool  1D3 3_kallisto.se/matsu-ta04 matsu-ta04
#5 matsu-ta05 cool  1D3 3_kallisto.se/matsu-ta05 matsu-ta05
#6 matsu-ta06 cool  1D3 3_kallisto.se/matsu-ta06 matsu-ta06

######## Generate sleuth object with liklihood ratio test (LRT)
so <- sleuth_prep(smpl, ~ Time)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
results_table_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
table(results_table_lrt[,"qval"] < 0.05)
#
#FALSE  TRUE 
#67130  1794 

#Optional, output PCA plot to check the data repeatablity
plot_pca(so, color_by = "Time")

########### Do the test again with Wald test (WT)
models(so)
#[  full  ]
#formula:  ~Time 
#data modeled:  obs_counts 
#transform sync'ed:  TRUE 
#coefficients:
#(Intercept)
#Time1D3
#[  reduced  ]
#formula:  ~1 
#data modeled:  obs_counts 
#transform sync'ed:  TRUE 
#coefficients:
#  (Intercept)

so <- sleuth_wt(so, "Time1D3")
results_table_wt <- sleuth_results(so, 'Time1D3')

################ Prepare output data
#Extract common DEG from two LRT & WT analyses
d5.lrt.sig_ids <- results_table_lrt$target_id[which(results_table_lrt$qval < 0.05)]
d5.wt.sig_ids <- results_table_wt$target_id[which(results_table_wt$qval < 0.05)]
shared_ids <- d5.wt.sig_ids[d5.wt.sig_ids %in% d5.lrt.sig_ids]
shared_results <- results_table_wt[results_table_wt$target_id %in% shared_ids,]
out <- shared_results[,c(1:4)]
out$log2FC <- -log(exp(out$b),2)   #<<- Add Log2FC column
out <- subset(out, select=-c(b))  #<<- Remove b column
out <- subset(out, out$log2FC >= 1 | out$log2FC <= -1) #<<- select by FC
head(out)
#             target_id          pval          qval     log2FC
#1 TraesCS4A02G162200.1  0.000000e+00  0.000000e+00  -4.527090
#2 TraesCS4D02G207500.1 2.489079e-240 8.577864e-236  -3.616269
#3 TraesCS3B02G135400.1 1.601526e-237 3.679451e-233 -10.595830
#4 TraesCS4A02G097900.2 1.650201e-219 2.843461e-215  -4.233299
#5 TraesCS3D02G118200.2 4.733792e-191 6.525437e-187 -10.216819
#6 TraesCS4B02G206700.1 8.324538e-166 9.562674e-162  -3.398584
shared_ids <- as.vector(out$target_id)  #<<-prepare new id list

#extract raw TPM data
tpm_matrix <- sleuth_to_matrix(so, 'obs_raw', 'tpm')
d <- tpm_matrix #for sure
d <- cbind(rownames(d), data.frame(d, row.names=NULL))
colnames(d)[1] <- "target_id"
tpm <- d[d$target_id %in% shared_ids,]
head(tpm)
#               target_id matsu.ta02 matsu.ta03 matsu.ta04 matsu.ta05 matsu.ta06
#534 TraesCS1A02G046800.1   1.399190   1.074606   3.230287   4.248172   3.532629
#662 TraesCS1A02G058400.1  11.172633  14.749862  44.213661  31.611015  39.869216
#663 TraesCS1A02G058400.2   2.415274   2.763716  11.629040   9.894266  12.424415
#664 TraesCS1A02G058400.3   0.000000   0.000000   1.686642   3.219718   1.845818
#893 TraesCS1A02G076900.1   4.402836   3.110870  40.815573  47.283302  45.454829
#894 TraesCS1A02G076900.2   0.000000   0.000000   5.442433   7.497449   5.826613

#merge two tables 
tmp <- merge(tpm, out, by="target_id", sort=T)
head(tmp)
#target_id matsu.ta02 matsu.ta03 matsu.ta04 matsu.ta05 matsu.ta06
#1 TraesCS1A02G046800.1   1.399190   1.074606   3.230287   4.248172   3.532629
#2 TraesCS1A02G058400.1  11.172633  14.749862  44.213661  31.611015  39.869216
#3 TraesCS1A02G058400.2   2.415274   2.763716  11.629040   9.894266  12.424415
#4 TraesCS1A02G058400.3   0.000000   0.000000   1.686642   3.219718   1.845818
#5 TraesCS1A02G076900.1   4.402836   3.110870  40.815573  47.283302  45.454829
#6 TraesCS1A02G076900.2   0.000000   0.000000   5.442433   7.497449   5.826613
#pval         qval    log2FC
#1 1.001831e-05 3.923308e-04 -1.532479
#2 1.883208e-07 1.194096e-05 -1.551239
#3 1.435158e-10 1.788731e-08 -2.099192
#4 2.987347e-19 1.045177e-16 -7.308665
#5 7.573474e-71 2.174975e-67 -3.559156
#6 7.691461e-55 1.656645e-51 -8.635240

#preprare "whole" output (just repeat the steps)
fdr <- results_table_wt[,c(1:4)]
fdr$log2FC <- -log(exp(fdr$b),2)   #<<- Add Log2FC column
fdr <- subset(fdr, select=-c(b))   #<<- Remove b column
all <- merge(d, fdr, by="target_id", sort=T)
head(all)

#write-down
write.csv(tmp, file="sleuth-DEG.csv", row.names=FALSE)
write.csv(all, file="sleuth-whole.csv", row.names=FALSE)

##############[Optional] Volcano-plot
#collect the data & roundup numbers
vc <- all[,c("target_id","qval","log2FC")]
vc$qval <- -log(vc$qval, 10)
vc$qval <- round(vc$qval, digits=1)
vc$log2FC <- round(vc$log2FC, digits=1)
head(vc)
#             target_id  qval log2FC
#1 TRAESCS3B02G135400.1  29.0    8.6
#2 TRAESCS5D02G068800.1 101.4    8.0
#3 TRAESCS5B02G059200.1   9.7    7.5
#4 TRAESCS5A02G057500.1   9.4    6.6
#5  TRAESCSU02G050500.1   4.9    6.6
#6  TRAESCSU02G196100.1  49.2    6.2

#Generate & cound DEGs
v1 <- subset(vc, vc$log2FC >= 1 & vc$qval > 1.3) ## UP-DEG
v2 <- subset(vc, vc$log2FC <= -1 & vc$qval > 1.3) ## DN-DEG
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
plot(v3$log2FC,v3$qval, pch=20, cex=0.4, col="grey70", ylim=c(0,40), xlim=c(-10,10), xlab="Log2FC", ylab="-log10FDR")
points(v1$log2FC,v1$qval, pch=20, cex=0.5, col="#ED0422")
points(v2$log2FC,v2$qval, pch=20, cex=0.5, col="#3498DB")
abline(v=0, lty=2)

# Add numbers of DEG on the plot
text(-10,39, up.deg, pos=4, cex=1.4, col="#ED0422")
text(-10,36, dn.deg, pos=4, cex=1.4, col="#3498DB")
