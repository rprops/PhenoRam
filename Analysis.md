---
title: "PhenoRam analysis"
author: "Cristina Garcia Timmermans & Ruben Props"
date: "Today"
output:
  html_document:
    code_folding: show
    highlight: haddock
    keep_md: yes
    theme: united
    toc: yes
    toc_float:
      collapsed: no
      smooth_scroll: yes
      toc_depth: 2
    css: report_styles.css
editor_options: 
  chunk_output_type: console
---



# Load libraries


```r
library("Phenoflow")
library("mclust")
library("plyr")
library("dplyr")
library("gridExtra")
library("tidyr")
library("phyloseq")
library("ggplot2")
library("gridExtra")
library("fpc")
library("RColorBrewer")
library("vegan")
library("tsne")
library("sandwich")
library("cluster")
library("grid")
library("egg")
library("MALDIquant")
library("reshape2")
library("caret") # for cross validation
library("foreach") # for parallelization
source("functions.R")
library("flowAI") # for denoising of FCM data
load("workspace.RData")

my.settings <- list(
  strip.background=list(col="transparent"),
  strip.border=list(col="transparent", cex=5),
  gate=list(col="black", fill="lightblue", alpha=0.2,border=NA,lwd=2),
  panel.background=list(col="lightgray"),
  background=list(col="white"))
```

# Cluster analysis
<!-- ## A. Hyperspec-normalized spectra -->
<!-- We will start with the hyperspec processed spectr:   -->

<!-- * Select number of stable clusters in data set -->
<!-- * Map back to original samples -->
<!-- * Calculate phenotypic diversity in each sample using Hill numbers -->

<!-- ```{r determine-clusters, fig.width = 7, fig.height= 6, dpi = 500, warning=FALSE} -->
<!-- # Choose if you want to run PCA prior to clustering -->
<!-- PCA <- FALSE -->

<!-- if(PCA == TRUE){ -->
<!--   # Perform PCA to reduce number of features in fingerprint -->
<!--   pca_bacteria <- prcomp(hs.norm) -->

<!--   # Only retain PC which explain 90% of the variance -->
<!--   thresh <- 0.9 -->
<!--   nr_pc_bacteria <- min(which((cumsum(vegan::eigenvals(pca_bacteria)/sum(vegan::eigenvals(pca_bacteria)))>thresh) == TRUE)) -->
<!--   pc_cluster_bacteria <- pca_bacteria$x[, 1:nr_pc_bacteria] -->
<!-- } else { -->
<!--   pc_cluster_bacteria <- hs.norm -->
<!-- } -->

<!-- # Evaluate number of robust clusters by means of silhouette index -->
<!-- # We limit the search to 50 clusters -->
<!-- tmp.si <- c() -->
<!-- for(i in 2:50){ -->
<!--   if(i%%10 == 0) cat(date(), paste0("---- at k =  ", i, "/",  nrow(pc_cluster_bacteria), "\n")) -->
<!--   tmp.si[i] <- pam(pc_cluster_bacteria, k = i)$silinfo$avg.width -->
<!-- } -->
<!-- nr_clusters_bacteria <- which(tmp.si == max(tmp.si, na.rm = TRUE)) -->

<!-- # Plot Silhouette index distribution -->
<!-- plot(tmp.si, type = "l", ylab = "Silhouette index",  -->
<!--      xlab = "Number of clusters") -->

<!-- # Cluster samples and export cluster labels -->
<!-- clusters_bacteria <- pam(pc_cluster_bacteria, k = nr_clusters_bacteria) -->

<!-- # Extract cluster labels -->
<!-- cluster_labels_pam <- data.frame(Sample = names(clusters_bacteria$clustering), -->
<!--                                       cluster_label = clusters_bacteria$clustering) -->

<!-- # Method 2: the Mclust( ) function in the mclust package selects the optimal model according to BIC for EM initialized by hierarchical clustering for parameterized Gaussian mixture models. -->
<!-- BIC = mclustBIC(pc_cluster_bacteria$spc, G = c(1:20)) -->
<!-- plot(BIC);lines(x = c(13, 13), y = c(0, 10e6), col = "red", lty = 2) -->
<!-- mc_fit <- Mclust(pc_cluster_bacteria, G = 13) -->

<!-- # plot(fit) # plot results  -->
<!-- summary(mc_fit) # display the best model -->

<!-- cluster_labels_mc <- data.frame(Sample = names(clusters_bacteria$clustering), -->
<!--                                       cluster_label = mc_fit$classification) -->

<!-- # To compare both clustering approaches: -->
<!-- # cluster.stats(dist(hs.norm), mc_fit$classification, clusters_bacteria$clustering) -->

<!-- # Extract count table (i.e. "operational phenotypic unit table") for each sample -->
<!-- OPU_hs_pam <- data.frame(table(cluster_labels_pam)) -->
<!-- # print(OPU_hs_pam) -->

<!-- OPU_hs_mc <- data.frame(table(cluster_labels_mc)) -->
<!-- # print(OPU_hs_mc) -->

<!-- # Merge cluster outputs in long format df -->
<!-- OPU_hs_merged <- rbind(OPU_hs_pam, OPU_hs_mc) -->
<!-- OPU_hs_merged <- data.frame(OPU_hs_merged, method = -->
<!--                               c(rep("PAM", nrow(OPU_hs_pam)), -->
<!--                                 rep("Mclust", nrow(OPU_hs_mc))),  -->
<!--                             replicate = do.call(rbind, strsplit(as.character(OPU_hs_merged$Sample), " "))[, 2], -->
<!--                             growth_phase = do.call(rbind, strsplit(as.character(OPU_hs_merged$Sample), " "))[, 1]) -->
<!-- colnames(OPU_hs_merged)[colnames(OPU_hs_merged) == "cluster_label"] <- "OPU" -->
<!-- ``` -->

<!-- # Plot OPU table -->
<!-- ```{r plot-clusters, fig.width = 7, fig.height= 6, dpi = 500, warning=FALSE} -->
<!-- # Plot according to metadata -->
<!-- p1 <- ggplot(OPU_hs_merged, aes(x = replicate, y = Freq, fill = OPU))+ -->
<!--   geom_bar(stat = "identity")+ -->
<!--   scale_fill_brewer(palette = "Paired")+ -->
<!--   theme_bw()+ -->
<!--   facet_grid(method ~ growth_phase, scales = "free")+ -->
<!--    theme(axis.title=element_text(size=16), strip.text=element_text(size=16), -->
<!--         legend.title=element_text(size=15),legend.text=element_text(size=14), -->
<!--         axis.text = element_text(size=14),title=element_text(size=20), -->
<!--         axis.text.x = element_text(angle = 45, hjust = 1), -->
<!--         strip.background=element_rect(fill=adjustcolor("lightgray",0.2)) -->
<!--         #,panel.grid.major = element_blank(), panel.grid.minor = element_blank() -->
<!--         ) -->

<!-- print(p1) -->
<!-- ``` -->

## B. Maldiquant-normalized spectra  

Same analysis as for hyperspec-normalized spectra.


```r
# Convert massSpectrum object to hyperspec
wv_mq <- mass(mq.norm[[1]])
matrix.spectra <- matrix(nrow=length(mq.norm), ncol = length(wv_mq))
for (i in 1:length(mq.norm)){
  matrix.spectra[i,] <- intensity(mq.norm[[i]])
}
hs.mq <- new("hyperSpec", spc = matrix.spectra, wavelength = wv_mq, labels = cell.name)

# Choose if you want to run PCA prior to clustering
PCA <- FALSE
PEAKS <- FALSE

if(PCA == TRUE){
  # Perform PCA to reduce number of features in fingerprint
  pca_bacteria <- prcomp(hs.mq)

  # Only retain PC which explain 90% of the variance
  thresh <- 0.9
  nr_pc_bacteria <- min(which((cumsum(vegan::eigenvals(pca_bacteria)/sum(vegan::eigenvals(pca_bacteria)))>thresh) == TRUE))
  pc_cluster_bacteria <- pca_bacteria$x[, 1:nr_pc_bacteria]
} else if(PEAKS == TRUE){
  # Run peak detection algorithm
    peaks <- detectPeaks(mq.norm, method="MAD", halfWindowSize=1, SNR=0.001)
    plot(mq.norm[[1]], xlim=c(600, 1800))
    points(peaks[[1]], col="red", pch=4)

    # Tolerance for wave number shift
    peaks <- binPeaks(peaks, tolerance = 0.002)
  
    # Filter out intensities at peak wave numbers
    peaks <- filterPeaks(peaks, minFrequency = 0.25)

    pc_cluster_bacteria <- intensityMatrix(peaks, mq.norm)
} else {
  pc_cluster_bacteria <- hs.mq
}

# Evaluate number of robust clusters by means of silhouette index
# We limit the search to 50 clusters
tmp.si <- c()
for(i in 2:50){
  if(i%%10 == 0) cat(date(), paste0("---- at k =  ", i, "/",  nrow(pc_cluster_bacteria), "\n"))
  tmp.si[i] <- pam(pc_cluster_bacteria, k = i)$silinfo$avg.width
}
```

```
## Mon Nov 27 15:22:52 2017 ---- at k =  10/536
## Mon Nov 27 15:22:54 2017 ---- at k =  20/536
## Mon Nov 27 15:22:58 2017 ---- at k =  30/536
## Mon Nov 27 15:23:04 2017 ---- at k =  40/536
## Mon Nov 27 15:23:16 2017 ---- at k =  50/536
```

```r
nr_clusters_bacteria <- which(tmp.si == max(tmp.si, na.rm = TRUE))

# Plot Silhouette index distribution
plot(tmp.si, type = "l", ylab = "Silhouette index", 
     xlab = "Number of clusters")
```

<img src="Figures/cached/determine-clusters-mq-1.png" style="display: block; margin: auto;" />

```r
# Cluster samples and export cluster labels
clusters_bacteria <- pam(pc_cluster_bacteria, k = nr_clusters_bacteria)

# Extract cluster labels
cluster_labels_pam <- data.frame(Sample = cell.name,
                                      cluster_label = clusters_bacteria$clustering)

# Method 2: the Mclust( ) function in the mclust package selects the optimal model according to BIC for EM initialized by hierarchical clustering for parameterized Gaussian mixture models.
if(PEAKS == TRUE){
  mc_fit <- Mclust(as.matrix(pc_cluster_bacteria))
} else {
  mc_fit <- Mclust(pc_cluster_bacteria, G = 13)
}

# plot(fit) # plot results 
summary(mc_fit) # display the best model
```

```
## ----------------------------------------------------
## Gaussian finite mixture model fitted by EM algorithm 
## ----------------------------------------------------
## 
## Mclust VII (spherical, varying volume) model with 13 components:
## 
##  log.likelihood   n   df     BIC     ICL
##         1285302 536 4354 2543243 2543243
## 
## Clustering table:
##  1  2  3  4  5  6  7  8  9 10 11 12 13 
## 37 49 26 20 36 17 41 75 55 59 43 34 44
```

```r
cluster_labels_mc <- data.frame(Sample = cell.name,
                                      cluster_label = mc_fit$classification)

# To compare both clustering approaches:
# cluster.stats(dist(hs.mq), mc_fit$classification, clusters_bacteria$clustering)

# Extract count table (i.e. "operational phenotypic unit table") for each sample
OPU_mq_pam <- data.frame(table(cluster_labels_pam))
# print(OPU_hs_pam)

OPU_mq_mc <- data.frame(table(cluster_labels_mc))
# print(OPU_hs_mc)

# Merge cluster outputs in long format df
OPU_mq_merged <- rbind(OPU_mq_pam, OPU_mq_mc)
OPU_mq_merged <- data.frame(OPU_mq_merged, method =
                              c(rep("PAM", nrow(OPU_mq_pam)),
                                rep("Mclust", nrow(OPU_mq_mc))), 
                            replicate = do.call(rbind, strsplit(as.character(OPU_mq_merged$Sample), " "))[, 2],
                            growth_phase = do.call(rbind, strsplit(as.character(OPU_mq_merged$Sample), " "))[, 1])
colnames(OPU_mq_merged)[colnames(OPU_mq_merged) == "cluster_label"] <- "OPU"
```

# Plot OPU table

```r
p2 <- ggplot(OPU_mq_merged, aes(x = replicate, y = Freq, fill = OPU))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  facet_grid(method ~ growth_phase, scales = "free")+
   theme(axis.title=element_text(size=16), strip.text=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )

print(p2)
```

<img src="Figures/cached/plot-clusters-mq-1.png" style="display: block; margin: auto;" />

# PhenoRam diversity

```r
# Format PAM clusters into otu tables
OPU_pam_table <- OPU_mq_merged %>% filter(method == "PAM") %>% select(c("Sample","OPU","Freq")) %>% tidyr::spread(OPU, Freq)
rownames(OPU_pam_table) <- OPU_pam_table$Sample
OPU_pam_table <- OPU_pam_table[, -1]
OPU_pam_tax <- as.matrix(data.frame(OPU = colnames(OPU_pam_table)))
rownames(OPU_pam_tax) <- OPU_pam_tax[,1]
OPU_pam_table2 <- phyloseq(otu_table(OPU_pam_table, taxa_are_rows = FALSE),
                           tax_table(OPU_pam_tax))

# Format Mclust clusters into otu tables
OPU_mclust_table <- OPU_mq_merged %>% filter(method == "Mclust") %>% select(c("Sample","OPU","Freq")) %>% tidyr::spread(OPU, Freq)
rownames(OPU_mclust_table) <- OPU_mclust_table$Sample
OPU_mclust_table <- OPU_mclust_table[, -1]
OPU_mclust_tax <- as.matrix(data.frame(OPU = colnames(OPU_mclust_table)))
rownames(OPU_mclust_tax) <- OPU_mclust_tax[,1]
OPU_mclust_table <- phyloseq(otu_table(OPU_mclust_table, taxa_are_rows = FALSE),
                           tax_table(OPU_mclust_tax))

div_ram_pam <- Diversity_16S(OPU_mclust_table, R = 100, brea = FALSE, 
                             parallel = TRUE, ncores = 3)
```

```
## 	**WARNING** this functions assumes that rows are samples and columns
##       	are taxa in your phyloseq object, please verify.
## Mon Nov 27 15:23:47 2017 	Using 3 cores for calculations
## Mon Nov 27 15:23:47 2017	Calculating diversity for sample 1/9 --- lag rep1
## Mon Nov 27 15:24:00 2017	Done with sample lag rep1
## Mon Nov 27 15:24:00 2017	Calculating diversity for sample 2/9 --- lag rep2
## Mon Nov 27 15:24:04 2017	Done with sample lag rep2
## Mon Nov 27 15:24:04 2017	Calculating diversity for sample 3/9 --- lag rep3
## Mon Nov 27 15:24:07 2017	Done with sample lag rep3
## Mon Nov 27 15:24:07 2017	Calculating diversity for sample 4/9 --- log rep1
## Mon Nov 27 15:24:11 2017	Done with sample log rep1
## Mon Nov 27 15:24:11 2017	Calculating diversity for sample 5/9 --- log rep2
## Mon Nov 27 15:24:14 2017	Done with sample log rep2
## Mon Nov 27 15:24:14 2017	Calculating diversity for sample 6/9 --- log rep3
## Mon Nov 27 15:24:18 2017	Done with sample log rep3
## Mon Nov 27 15:24:18 2017	Calculating diversity for sample 7/9 --- stat rep1
## Mon Nov 27 15:24:21 2017	Done with sample stat rep1
## Mon Nov 27 15:24:21 2017	Calculating diversity for sample 8/9 --- stat rep2
## Mon Nov 27 15:24:24 2017	Done with sample stat rep2
## Mon Nov 27 15:24:24 2017	Calculating diversity for sample 9/9 --- stat rep3
## Mon Nov 27 15:24:27 2017	Done with sample stat rep3
## Mon Nov 27 15:24:27 2017 	Closing workers
## Mon Nov 27 15:24:27 2017 	Done with all 9 samples
```

```r
div_ram_mclust <- Diversity_16S(OPU_pam_table2, R = 100, brea = FALSE, 
                             parallel = TRUE, ncores = 3)
```

```
## 	**WARNING** this functions assumes that rows are samples and columns
##       	are taxa in your phyloseq object, please verify.
## Mon Nov 27 15:24:28 2017 	Using 3 cores for calculations
## Mon Nov 27 15:24:28 2017	Calculating diversity for sample 1/9 --- lag rep1
## Mon Nov 27 15:24:39 2017	Done with sample lag rep1
## Mon Nov 27 15:24:39 2017	Calculating diversity for sample 2/9 --- lag rep2
## Mon Nov 27 15:24:42 2017	Done with sample lag rep2
## Mon Nov 27 15:24:42 2017	Calculating diversity for sample 3/9 --- lag rep3
## Mon Nov 27 15:24:45 2017	Done with sample lag rep3
## Mon Nov 27 15:24:45 2017	Calculating diversity for sample 4/9 --- log rep1
## Mon Nov 27 15:24:48 2017	Done with sample log rep1
## Mon Nov 27 15:24:48 2017	Calculating diversity for sample 5/9 --- log rep2
## Mon Nov 27 15:24:51 2017	Done with sample log rep2
## Mon Nov 27 15:24:51 2017	Calculating diversity for sample 6/9 --- log rep3
## Mon Nov 27 15:24:54 2017	Done with sample log rep3
## Mon Nov 27 15:24:54 2017	Calculating diversity for sample 7/9 --- stat rep1
## Mon Nov 27 15:24:57 2017	Done with sample stat rep1
## Mon Nov 27 15:24:57 2017	Calculating diversity for sample 8/9 --- stat rep2
## Mon Nov 27 15:25:00 2017	Done with sample stat rep2
## Mon Nov 27 15:25:00 2017	Calculating diversity for sample 9/9 --- stat rep3
## Mon Nov 27 15:25:03 2017	Done with sample stat rep3
## Mon Nov 27 15:25:03 2017 	Closing workers
## Mon Nov 27 15:25:03 2017 	Done with all 9 samples
```

```r
div_ram_merged <- data.frame(Sample = rep(rownames(div_ram_pam),2),
                             rbind(div_ram_pam, div_ram_mclust),
                             method = c(rep("pam", nrow(div_ram_pam)),
                                        rep("mclust", nrow(div_ram_mclust)))
                             )
div_ram_merged <- div_ram_merged[, -c(4:7)]
div_ram_merged$Sample <- as.character(div_ram_merged$Sample)

# Merge with metadata
div_ram_merged$GrowthPhase <- do.call(rbind, strsplit(div_ram_merged$Sample, " "))[,1]

# Plot results
p_ram_div_pam <- div_ram_merged %>% filter(method == "pam") %>% 
  ggplot(aes(x = GrowthPhase, y = D2, fill = GrowthPhase))+
  geom_point(shape = 21, size = 4)+
  geom_boxplot(alpha = 0.4)+
  ggplot2::theme_bw()+
     theme(axis.title=element_text(size=16), strip.text=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )+
  scale_fill_brewer(palette = "Accent")+
  geom_errorbar(aes(ymin = D2 - sd.D2, ymax = D2 + sd.D2), width = 0.025)+
  guides(fill = FALSE)+
  ylab(expression("Phenotypic diversity - D2 (Raman)"))

print(p_ram_div_pam)
```

<img src="Figures/cached/calculate-pheno-div-ram-1.png" style="display: block; margin: auto;" />

```r
p_ram_div_mclust <- div_ram_merged %>% filter(method == "mclust") %>% 
  ggplot(aes(x = GrowthPhase, y = D2, fill = GrowthPhase))+
  geom_point(shape = 21, size = 4)+
  geom_boxplot(alpha = 0.4)+
  ggplot2::theme_bw()+
     theme(axis.title=element_text(size=16), strip.text=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )+
  scale_fill_brewer(palette = "Accent")+
  geom_errorbar(aes(ymin = D2 - sd.D2, ymax = D2 + sd.D2), width = 0.025)+
  guides(fill = FALSE)+
  ylab(expression("Phenotypic diversity - D2 (Raman)"))

print(p_ram_div_mclust)
```

<img src="Figures/cached/calculate-pheno-div-ram-2.png" style="display: block; margin: auto;" />



# Contrast analysis
## Hyperspec normalized  


```r
# ram_contrast(hyprs = hs.norm, comp1 = c("LB rep1", "LB rep2", "LB rep3"), 
# comp2 = c("NB rep1","NB rep2","NB rep3"))
ram.hs_lag_log <- ram_contrast(hs.norm, comp1 = c("lag rep1", "lag rep2", "lag rep3"),
              comp2 = c("log rep1", "log rep2", "log rep3"), plot = FALSE)
```

```
## -----------------------------------------------------------------------------------------------------
##  
## 	 Your cells are distributed over these samples:
## 
##  Samples
##  lag rep1  lag rep2  lag rep3  log rep1  log rep2  log rep3 stat rep1 
##        60        62        61        58        60        59        59 
## stat rep2 stat rep3 
##        56        61 
## -----------------------------------------------------------------------------------------------------
##  
## 	 Returning contrasts between mean spectra for 183 cells of
##  c("lag rep1", "lag rep2", "lag rep3")
## 	 and 177 cells of
##  c("log rep1", "log rep2", "log rep3")
## -----------------------------------------------------------------------------------------------------
## 
```

```r
ram.hs_lag_stat <- ram_contrast(hs.norm, comp1 = c("lag rep1", "lag rep2", "lag rep3"),
              comp2 = c("stat rep1", "stat rep2", "stat rep3"), plot = FALSE)
```

```
## -----------------------------------------------------------------------------------------------------
##  
## 	 Your cells are distributed over these samples:
## 
##  Samples
##  lag rep1  lag rep2  lag rep3  log rep1  log rep2  log rep3 stat rep1 
##        60        62        61        58        60        59        59 
## stat rep2 stat rep3 
##        56        61 
## -----------------------------------------------------------------------------------------------------
##  
## 	 Returning contrasts between mean spectra for 183 cells of
##  c("lag rep1", "lag rep2", "lag rep3")
## 	 and 176 cells of
##  c("stat rep1", "stat rep2", "stat rep3")
## -----------------------------------------------------------------------------------------------------
## 
```

```r
ram.hs_log_stat <- ram_contrast(hs.norm, comp1 = c("log rep1", "log rep2", "log rep3"),
              comp2 = c("stat rep1", "stat rep2", "stat rep3"), plot = FALSE)
```

```
## -----------------------------------------------------------------------------------------------------
##  
## 	 Your cells are distributed over these samples:
## 
##  Samples
##  lag rep1  lag rep2  lag rep3  log rep1  log rep2  log rep3 stat rep1 
##        60        62        61        58        60        59        59 
## stat rep2 stat rep3 
##        56        61 
## -----------------------------------------------------------------------------------------------------
##  
## 	 Returning contrasts between mean spectra for 177 cells of
##  c("log rep1", "log rep2", "log rep3")
## 	 and 176 cells of
##  c("stat rep1", "stat rep2", "stat rep3")
## -----------------------------------------------------------------------------------------------------
## 
```

```r
ram.hs_merged <- data.frame(rbind(ram.hs_lag_log, ram.hs_lag_stat, ram.hs_log_stat),
                            Comparison = rep(c("lag-log", "lag-stat", "log-stat"), 
                                               each = nrow(ram.hs_lag_stat))
)

v.hs <- ggplot2::ggplot(ram.hs_merged, ggplot2::aes(x = Wavenumber, y = Density, fill = Density))+
  ggplot2::geom_point(shape = 21, colour="black", alpha = 1.0,
                          size = 3)+
  geom_line(color = "black", alpha = 0.6)+
  facet_grid(.~Comparison)+
  ggplot2::scale_fill_distiller(palette="RdBu", na.value="white", limits = c(-0.22,0.22)) +
  scale_x_continuous(breaks = seq(600,1800,200), labels = seq(600,1800,200))+
  ggplot2::theme_bw()+
     theme(axis.title=element_text(size=16), strip.text=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )

print(v.hs)
```

<img src="Figures/cached/plot-contrasts-hs-1.png" style="display: block; margin: auto;" />

## Maldiquant normalized  


```r
# ram_contrast(hyprs = hs.norm, comp1 = c("LB rep1", "LB rep2", "LB rep3"), 
# comp2 = c("NB rep1","NB rep2","NB rep3"))
ram.mq_lag_log <- ram_contrast(hs.mq, comp1 = c("lag rep1", "lag rep2", "lag rep3"),
              comp2 = c("log rep1", "log rep2", "log rep3"), plot = FALSE)
```

```
## -----------------------------------------------------------------------------------------------------
##  
## 	 Your cells are distributed over these samples:
## 
##  Samples
##  lag rep1  lag rep2  lag rep3  log rep1  log rep2  log rep3 stat rep1 
##        60        62        61        58        60        59        59 
## stat rep2 stat rep3 
##        56        61 
## -----------------------------------------------------------------------------------------------------
##  
## 	 Returning contrasts between mean spectra for 183 cells of
##  c("lag rep1", "lag rep2", "lag rep3")
## 	 and 177 cells of
##  c("log rep1", "log rep2", "log rep3")
## -----------------------------------------------------------------------------------------------------
## 
```

```r
ram.mq_lag_stat <- ram_contrast(hs.mq, comp1 = c("lag rep1", "lag rep2", "lag rep3"),
              comp2 = c("stat rep1", "stat rep2", "stat rep3"), plot = FALSE)
```

```
## -----------------------------------------------------------------------------------------------------
##  
## 	 Your cells are distributed over these samples:
## 
##  Samples
##  lag rep1  lag rep2  lag rep3  log rep1  log rep2  log rep3 stat rep1 
##        60        62        61        58        60        59        59 
## stat rep2 stat rep3 
##        56        61 
## -----------------------------------------------------------------------------------------------------
##  
## 	 Returning contrasts between mean spectra for 183 cells of
##  c("lag rep1", "lag rep2", "lag rep3")
## 	 and 176 cells of
##  c("stat rep1", "stat rep2", "stat rep3")
## -----------------------------------------------------------------------------------------------------
## 
```

```r
ram.mq_log_stat <- ram_contrast(hs.mq, comp1 = c("log rep1", "log rep2", "log rep3"),
              comp2 = c("stat rep1", "stat rep2", "stat rep3"), plot = FALSE)
```

```
## -----------------------------------------------------------------------------------------------------
##  
## 	 Your cells are distributed over these samples:
## 
##  Samples
##  lag rep1  lag rep2  lag rep3  log rep1  log rep2  log rep3 stat rep1 
##        60        62        61        58        60        59        59 
## stat rep2 stat rep3 
##        56        61 
## -----------------------------------------------------------------------------------------------------
##  
## 	 Returning contrasts between mean spectra for 177 cells of
##  c("log rep1", "log rep2", "log rep3")
## 	 and 176 cells of
##  c("stat rep1", "stat rep2", "stat rep3")
## -----------------------------------------------------------------------------------------------------
## 
```

```r
ram.mq_merged <- data.frame(rbind(ram.mq_lag_log, ram.mq_lag_stat, ram.mq_log_stat),
                            Comparison = rep(c("lag-log", "lag-stat", "log-stat"), 
                                               each = nrow(ram.mq_lag_stat))
)

v.mq <- ggplot2::ggplot(ram.mq_merged, ggplot2::aes(x = Wavenumber, y = Density, fill = Density))+
  ggplot2::geom_point(shape = 21, colour="black", alpha = 1.0,
                          size = 3)+
  geom_line(color = "black", alpha = 0.6)+
  facet_grid(.~Comparison)+
  ggplot2::scale_fill_distiller(palette="RdBu", na.value="white", limits = c(-0.22,0.22)) +
  scale_x_continuous(breaks = seq(600,1800,200), labels = seq(600,1800,200))+
  ggplot2::theme_bw()+
     theme(axis.title=element_text(size=16), strip.text=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )

print(v.mq)
```

<img src="Figures/cached/plot-contrasts-mq-1.png" style="display: block; margin: auto;" />

# Flow cytometry data

### Alpha diversity analysis  



```r
# Import data in FCS format
fs <- flowCore::read.flowSet(path = "./FCSfiles", pattern = ".fcs")

# Extract metadata from sample names
meta_fs <- do.call(rbind, strsplit(flowCore::sampleNames(fs), " "))[,2]
meta_fs <- data.frame(Sample = flowCore::sampleNames(fs),
                      Organism = do.call(rbind, strsplit(meta_fs, "_"))[,1],
                      GrowthPhase = do.call(rbind, strsplit(meta_fs, "_"))[,2],
                      Replicate= do.call(rbind, strsplit(meta_fs, "_"))[,3],
                      Dilution= do.call(rbind, strsplit(meta_fs, "_"))[,4])
meta_fs$Dilution <- as.numeric(gsub(".fcs", "", meta_fs$Dilution))

# Transform data with asinh
# Select phenotypic features of interest and transform parameters
fs <- flowCore::transform(fs,`FL1-H`=asinh(`FL1-H`), 
                                   `SSC-H`=asinh(`SSC-H`), 
                                   `FL3-H`=asinh(`FL3-H`), 
                                   `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")

# Denoise data
### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(7.75,7.75,14,14,3,7.75,14,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

###  Gating quality check
xyplot(`FL3-H` ~ `FL1-H`, data=fs[1], filter=polyGate1,
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(6,16))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)
```

<img src="Figures/cached/FCM-analysis-1-1.png" style="display: block; margin: auto;" />

```r
### Isolate only the cellular information based on the polyGate1
fs <- Subset(fs, polyGate1)

#Normalize
summary <- fsApply(x = fs, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
maxval <- max(summary[,9])
mytrans <- function(x) x/maxval
fs <- transform(fs,`FL1-H`=mytrans(`FL1-H`),
                                  `FL3-H`=mytrans(`FL3-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))

# Calculate phenotypic diversity
fs_div <- Diversity_rf(fs, param = param, R.b = 100, R = 100, cleanFCS = FALSE,
                       parallel = TRUE, ncores = 3)
```

```
## --- parameters are already normalized at: 1
## Mon Nov 27 15:25:19 2017 --- Using 3 cores for calculations
## Mon Nov 27 15:29:08 2017 --- Closing workers
## Mon Nov 27 15:29:08 2017 --- Alpha diversity metrics (D0,D1,D2) have been computed after 100 bootstraps
## -----------------------------------------------------------------------------------------------------
## 
```

```r
## Plot alpha diversity vs growth phase
fs_div <- dplyr::left_join(fs_div, meta_fs, by = c("Sample_names"="Sample"))

p_fs_div <- ggplot(fs_div, aes(x = GrowthPhase, y = D2, fill = GrowthPhase))+
  geom_point(shape = 21, size = 4)+
  geom_boxplot(alpha = 0.4)+
  ggplot2::theme_bw()+
     theme(axis.title=element_text(size=16), strip.text=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )+
  scale_fill_brewer(palette = "Accent")+
  geom_errorbar(aes(ymin = D2 - sd.D2, ymax = D2 + sd.D2), width = 0.025)+
  guides(fill = FALSE)+
  ylab(expression("Phenotypic diversity - D2 (FCM)"))

# Compare FCM diversity (left hand side) with Raman diversity (right hand side)
grid.arrange(p_fs_div, p_ram_div_mclust, ncol = 2)
```

<img src="Figures/cached/FCM-analysis-1-2.png" style="display: block; margin: auto;" />

### Contrast analysis  


```r
# Calculate fingerprint
fbasis <- flowBasis(fs, param, nbin=128, 
                   bw=0.01,normalize=function(x) x)

# Calculate contrasts
fp_c_lag_log <- fp_contrasts(fbasis, comp1 = meta_fs$GrowthPhase=="lag", comp2 = 
               meta_fs$GrowthPhase=="log", thresh = 0.05)
```

```
## 	Region used for contrasts 1 16384
## 	Returning contrasts for A01 Ecoli_lag_rep1_1000.fcs A01 Ecoli_log_rep1_10000.fcs
##  	Returning contrasts for A02 Ecoli_lag_rep2_1000.fcs A02 Ecoli_log_rep2_10000.fcs
##  	Returning contrasts for A03 Ecoli_lag_rep3_1000.fcs A03 Ecoli_log_rep3_100000.fcs
```

```r
fp_c_lag_stat <- fp_contrasts(fbasis, comp1 = meta_fs$GrowthPhase=="lag", comp2 = 
               meta_fs$GrowthPhase=="stat", thresh = 0.05)
```

```
## 	Region used for contrasts 1 16384
## 	Returning contrasts for A01 Ecoli_lag_rep1_1000.fcs A01 Ecoli_stat_rep1_10000.fcs
##  	Returning contrasts for A02 Ecoli_lag_rep2_1000.fcs A02 Ecoli_stat_rep2_10000.fcs
##  	Returning contrasts for A03 Ecoli_lag_rep3_1000.fcs A03 Ecoli_stat_rep3_10000.fcs
```

```r
fp_c_log_stat <- fp_contrasts(fbasis, comp1 = meta_fs$GrowthPhase=="log", comp2 = 
               meta_fs$GrowthPhase=="stat", thresh = 0.05)
```

```
## 	Region used for contrasts 1 16384
## 	Returning contrasts for A01 Ecoli_log_rep1_10000.fcs A01 Ecoli_stat_rep1_10000.fcs
##  	Returning contrasts for A02 Ecoli_log_rep2_10000.fcs A02 Ecoli_stat_rep2_10000.fcs
##  	Returning contrasts for A03 Ecoli_log_rep3_100000.fcs A03 Ecoli_stat_rep3_10000.fcs
```

```r
fp_c_merge <- data.frame(rbind(fp_c_lag_log, fp_c_lag_stat, fp_c_log_stat),
                         Comparison = c(rep("lag-log", nrow(fp_c_lag_log)),
                                        rep("lag-stat", nrow(fp_c_lag_stat)),
                                        rep("log-stat", nrow(fp_c_log_stat))))
# Plot contrasts
v.fp_c <- ggplot2::ggplot(fp_c_merge, ggplot2::aes(x = FL1.H, y = FL3.H, fill = Density, 
                                                   z = Density))+
  geom_tile(aes(fill=Density)) + 
  geom_point(colour="gray", alpha=0.7)+
  scale_fill_distiller(palette="RdYlBu", na.value="white") + 
  stat_contour(aes(fill=..level..), geom="polygon", binwidth=0.1)+
  theme_bw()+
  # geom_line()+
  facet_grid(.~Comparison)+
  # scale_x_continuous(breaks = seq(600,1800,200), labels = seq(600,1800,200))+
  theme(axis.title=element_text(size=16), strip.text=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )+
  labs(x="Green fluorescence intensity (a.u.)", y="Red fluorescence intensity (a.u.)")

print(v.fp_c)
```

<img src="Figures/cached/FCM-analysis-2-1.png" style="display: block; margin: auto;" />

### Beta diversity analysis


```r
# Calculate beta diversity
beta.div <- beta_div_fcm(fbasis, ord.type="PCoA")

# Plot ordination
plot_beta_fcm(beta.div, color = meta_fs$GrowthPhase, labels="Growth phase") + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.5)+
  ggtitle("")
```

<img src="Figures/cached/FCM-analysis-3-1.png" style="display: block; margin: auto;" />

```r
# Juxtaposition with beta-diversity based on Raman data
```

# tSNE

## Flow cytometry


```r
# Calculate the dissimilarity matrix
dist_fbasis_sample <- dist(fbasis@basis)
tsne_fbasis <- tsne::tsne(dist_fbasis_sample, 
                          perplexity = 5)
```

```
## sigma summary: Min. : 12.5709245723508 |1st Qu. : 19.2135377420559 |Median : 20.9792132965026 |Mean : 20.9191497359823 |3rd Qu. : 22.082868101704 |Max. : 28.1623460431449 |
```

```
## Epoch: Iteration #100 error is: 16.2166585199963
```

```
## Epoch: Iteration #200 error is: 0.800916250981897
```

```
## Epoch: Iteration #300 error is: 0.521594437296171
```

```
## Epoch: Iteration #400 error is: 0.494656031369735
```

```
## Epoch: Iteration #500 error is: 0.480375937539916
```

```
## Epoch: Iteration #600 error is: 0.343160688367551
```

```
## Epoch: Iteration #700 error is: 0.227939522084319
```

```
## Epoch: Iteration #800 error is: 0.106205994334686
```

```
## Epoch: Iteration #900 error is: 0.0898148877257923
```

```
## Epoch: Iteration #1000 error is: 0.0660148814083378
```

```r
tsne_fbasis <- data.frame(tsne_fbasis, Sample = rownames(fbasis@basis))
colnames(tsne_fbasis)[1:2] <- c("Axis1", "Axis2")
tsne_fbasis <- dplyr::left_join(tsne_fbasis, meta_fs, by = "Sample")

p_tsne_fbasis <- ggplot2::ggplot(tsne_fbasis, ggplot2::aes(x = Axis1, y = Axis2, 
        fill = GrowthPhase)) + 
  ggplot2::geom_point(alpha = 0.7, 
        size = 8, shape = 21, colour = "black") + 
  ggplot2::scale_fill_manual(values = c("#a65628", "red", 
            "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta")) +
  ggplot2::labs(x = paste0("Axis1"), y = paste0("Axis2"))+
  theme_bw()+
    theme(axis.title=element_text(size=16), strip.text=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2)))

print(p_tsne_fbasis)
```

<img src="Figures/cached/tsne-1-1.png" style="display: block; margin: auto;" />

## Raman data


```r
# At single cell level
dist_ram_cells <- dist(hs.mq@data)
tsne_ram.c <- tsne::tsne(dist_ram_cells)
```

```
## sigma summary: Min. : 0.0131363067079677 |1st Qu. : 0.0200191227308468 |Median : 0.021848832655061 |Mean : 0.0239405053832849 |3rd Qu. : 0.025763593149787 |Max. : 0.0505763424606676 |
```

```
## Epoch: Iteration #100 error is: 10.3473032362903
```

```
## Epoch: Iteration #200 error is: 0.276972705489419
```

```
## Epoch: Iteration #300 error is: 0.234388809023767
```

```
## Epoch: Iteration #400 error is: 0.218658452883787
```

```
## Epoch: Iteration #500 error is: 0.213552120191546
```

```
## Epoch: Iteration #600 error is: 0.210866926308676
```

```
## Epoch: Iteration #700 error is: 0.209103898440317
```

```
## Epoch: Iteration #800 error is: 0.207774749990583
```

```
## Epoch: Iteration #900 error is: 0.206720669680325
```

```
## Epoch: Iteration #1000 error is: 0.205859746088686
```

```r
tsne_ram.c <- data.frame(tsne_ram.c, GrowthPhase = do.call(rbind, strsplit(cell.name, " "))[,1], replicate = do.call(rbind, strsplit(cell.name, " "))[,2])
colnames(tsne_ram.c)[1:2] <- c("Axis1", "Axis2")

p_tsne_ram.c <- ggplot2::ggplot(tsne_ram.c, ggplot2::aes(x = Axis1, y = Axis2, 
        fill = GrowthPhase, shape = Replicate)) + 
    ggplot2::scale_fill_manual(values = c("#a65628", "red", 
            "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta")) +
  ggplot2::geom_point(alpha = 0.7, 
        size = 4, colour = "black") + 
  scale_shape_manual(values = c(21,22,24))+
  ggplot2::labs(x = paste0("Axis1"), y = paste0("Axis2"))+
  theme_bw()+
  theme(axis.title=element_text(size=16), strip.text=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2)))

print(p_tsne_ram.c)
```

<img src="Figures/cached/tsne-2-1.png" style="display: block; margin: auto;" />
