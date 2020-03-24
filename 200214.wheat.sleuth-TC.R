##2020.02.14
## ALRC wheat mRNA-seq
## KALLISTO-SlEUTH, Timecourse DEG & hclust-heatmap

library(sleuth)
s2c <- read.table("meta(se).tmp.txt", header=TRUE, stringsAsFactors=FALSE)
smpl <- s2c[c(3:43),] #exclude the 0-hr data

############# Likelihood ratio test (LRT)
head(smpl)
#    sample cond Time                     path        run
#3    ta04 cool    3 2_kallisto.pe/matsu-ta04 matsu-ta04
#4    ta06 cool    3 2_kallisto.pe/matsu-ta06 matsu-ta06
#5    ta07 cool   24 2_kallisto.pe/matsu-ta07 matsu-ta07
#6    ta08 cool   24 2_kallisto.pe/matsu-ta08 matsu-ta08
#7    ta09 cool   24 2_kallisto.pe/matsu-ta09 matsu-ta09
#8    ta10 cool   27 2_kallisto.pe/matsu-ta10 matsu-ta10
so <- sleuth_prep(smpl)
so <- sleuth_fit(so, ~cond + Time + cond:Time, 'full')
so <- sleuth_fit(so, ~cond + Time, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
results_table_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
table(results_table_lrt[,"qval"] < 0.05)
#FALSE  TRUE 
#53658   292 from PE data
#63174   410 from SR data

#[Optional]PCA plot
plot_pca(so, color_by = "Time")
plot_pca(so, color_by = "cond")

##################generate output table
#DEG data
sig_ids <- results_table_lrt$target_id[which(results_table_lrt$qval < 0.05)]
sig_table <- results_table_lrt[results_table_lrt$target_id %in% sig_ids,][1:3]
head(sig_table)
#target_id         pval         qval
#1 TRAESCS2B02G362400.1 3.450217e-12 1.861392e-07
#2 TRAESCS6B02G302400.1 1.960301e-10 5.287912e-06
#3 TRAESCS7A02G529900.1 1.255258e-09 1.693029e-05
#4 TRAESCS7B02G446900.1 1.238979e-09 1.693029e-05
#5  TRAESCSU02G029400.1 2.523922e-09 2.723312e-05
#6 TRAESCS3B02G563100.1 7.893389e-09 7.097472e-05

#Extract TPM data
tpm_matrix <- sleuth_to_matrix(so, 'obs_raw', 'tpm')
d <- tpm_matrix #for sure
d <- cbind(rownames(d), data.frame(d, row.names=NULL))
colnames(d)[1] <- "target_id"
tpm <- d[d$target_id %in% sig_ids,]
dim(tpm)
#[1] 292  42

#Merging & outputs
out <- merge(tpm, sig_table, by="target_id", sort=T)
dim(out)
#[1] 292  44
write.csv(out, file="sleuth-timecourse.csv", row.names=FALSE)

###### Hclust & Heatmap
#Need to express the data to fold-change against 0-hr point
#re-call the whole data & convert TPM to log2FC
so <- sleuth_prep(s2c)
so_tpm <- sleuth_to_matrix(so, 'obs_raw', 'tpm')
so_tpm <- cbind(rownames(so_tpm), data.frame(so_tpm, row.names=NULL))
colnames(so_tpm)[1] <- "target_id"
deg_tpm <- so_tpm[so_tpm$target_id %in% sig_ids,]

means <- as.vector(rowMeans(deg_tpm[,c("ta02","ta03")])) #generae a vector of 0-h mean
deg_fc <- subset(deg_tpm, select=-c(target_id, ta02, ta03)) # remove not-number columns
deg_fc <- log((deg_fc+1e-12)/(means+1e-12),2)   #Log2FC 1e-12 to prevent Na,Inf

##### gplot drawing
library(gplots)
rownames(deg_fc) <- deg_tpm$target_id #Give rowname and gene_id
d <- data.matrix(deg_fc)

# Set color range
colors = c(seq(-5,-0.5,length=101),seq(-0.5,0.5,length=101),seq(0.5,5,length=101))
colors = unique(colors) ## break vector should composite unique numbers
length(colors)
my_palette <- colorRampPalette(c("green", "black", "red"))(n = 300)

# draw plot
heatmap.2(d, col=my_palette, breaks=colors, density.info="none", trace="none", dendrogram=c("row"), symm=F,symkey=F,symbreaks=T, scale="none", Colv=FALSE)

########### DEG Clustering by cclust (screening N)
library(cclust)
d_matrix <- as.matrix(deg_fc)
d_df <- as.data.frame(d_matrix)
for(i in 2:12){
  kmer <- cclust(d_matrix, i, 100, method="kmeans")
  d_df[[paste("FC_k",i,sep="")]] <- kmer$cluster
}

#Writedown the output
write.csv(d_df, file="sleuth-tc.cclust.csv", row.names=TRUE)


###Seperate clusters & draw line-plots for each
#in case of k=6; seperate tables from df1 to df6

for(i in 1:6){namc <- paste("dfc",i,sep="")
  namh <- paste("dfh",i,sep="")
  assign(namc, subset(d_df, d_df$FC_k6==i)[,c(1:20)])
  assign(namh, subset(d_df, d_df$FC_k6==i)[,c(21:41)])
}
dim(dfh5)
#[1]  13 21

require(data.table)
dfListc <- list(dfc1, dfc2, dfc3, dfc4, dfc5, dfc6) 
dfListh <- list(dfh1, dfh2, dfh3, dfh4, dfh5, dfh6)
varc <- grep("dfc", ls(), value=TRUE)
varh <- grep("dfh", ls(), value=TRUE)
for(i in c(1:length(dfListc))){
  tm1 <- rowMeans(dfListc[[i]][, c(1:2)])
  tm2 <- rowMeans(dfListc[[i]][, c(3:5)])
  tm3 <- rowMeans(dfListc[[i]][, c(6:7)])
  tm4 <- rowMeans(dfListc[[i]][, c(8:10)])
  tm5 <- rowMeans(dfListc[[i]][, c(11:12)])
  tm6 <- rowMeans(dfListc[[i]][, c(13:15)])
  tm7 <- rowMeans(dfListc[[i]][, c(16:17)])
  tm8 <- rowMeans(dfListc[[i]][, c(18:20)])
  assign(varc[i], data.frame(tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8))
}
for(i in c(1:length(dfListh))){
  tm1 <- rowMeans(dfListh[[i]][, c(1:3)])
  tm2 <- rowMeans(dfListh[[i]][, c(4:5)])
  tm3 <- rowMeans(dfListh[[i]][, c(6:7)])
  tm4 <- rowMeans(dfListh[[i]][, c(8:10)])
  tm5 <- rowMeans(dfListh[[i]][, c(11:13)])
  tm6 <- rowMeans(dfListh[[i]][, c(14:15)])
  tm7 <- rowMeans(dfListh[[i]][, c(16:18)])
  tm8 <- rowMeans(dfListh[[i]][, c(19:21)])
  assign(varh[i], data.frame(tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8))
}

#Plot with matplot
hlimit <- c(min(d_df)-5, max(d_df)+5) 

matplot(t(dfc1), main="Cluster1",type="l", col=c("grey40","grey60"), lty=1, lwd=1, ylim = hlimit)
matplot(t(dfh1), type="l", col=c("red","magenta"), lty=1, lwd=1, ylim = hlimit, add=TRUE)







#Wald test
models(so)
#[  full  ]
#formula:  ~cond + Time + cond:Time 
#data modeled:  obs_counts 
#transform sync'ed:  TRUE 
#coefficients:
#(Intercept)
#condhot
#Time
#condhot:Time
#[  reduced  ]
#formula:  ~cond + Time 
#data modeled:  obs_counts 
#transform sync'ed:  TRUE 
#coefficients:
#  (Intercept)
#condhot
#Time






### Generate the spline model
library(splines)
Time <- s2c$Time
Time <- as.numeric(Time) ##願掛け
full_design <- model.matrix(formula(~ ns(Time, df = 3))) ## Adjust df as (No.-2)
head(full_design)
#(Intercept) ns(Time, df = 3)1 ns(Time, df = 3)2 ns(Time, df = 3)3
#1            1        0.00000000        0.00000000        0.00000000
#2            1        0.00000000        0.00000000        0.00000000
#3            1       -0.01929288        0.05725819       -0.03790331
#4            1       -0.01929288        0.05725819       -0.03790331
#5            1       -0.10580188        0.40691256       -0.26936465
#6            1       -0.10580188        0.40691256       -0.26936465

### Run sleuth based on the model
so <- sleuth_prep(s2c, full_model = full_design)
#52953 targets passed the filter
so <- sleuth_fit(so)
so <- sleuth_fit(so, formula = ~ 1, fit_name = "reduced")
so <- sleuth_lrt(so, "reduced", "full")

