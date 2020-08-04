#RaceID4 
library(tidyverse)
library(viridis)
library(RaceID)
library(Seurat)
library(Matrix)
library(clustree)

date = Sys.Date()

#load andnormalize dataset 
all <- read.delim("data/GSE146926_sample_TPM_mouse_microglia_region_May2018.txt", sep = "\t")
all <- all[!duplicated(all$gene_name),]
rownames(all) <- all$gene_name

all <- all[,-1] %>%
  CreateSeuratObject(project="wt", min.cells = 10, min.features = 500) %>%
  SCTransform(variable.features.n = 10000) %>%
  RunPCA(features=grep("^(Gm|Rpl|Rps|Fos|Jun|mt-|Dusp|Zfp36|Hsp)|(Rik)$", VariableFeatures(.), value=T, invert = T)) 

#
ElbowPlot(all)
all<- all %>% 
  RunUMAP(dims=1:15) %>%
  FindNeighbors(dims=1:15) %>%
  FindClusters(resolution=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1))

#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(all)

pdf("plots/others/overview-cluster-resolutions.pdf")
clustree(all)
dev.off()

all<-FindClusters(all,resolution=.3)


DimPlot(all, label = TRUE) + NoLegend()

save(all, file = "data/all.RData")

#find cluster markers
all.markers<-FindAllMarkers(all,only.pos=T,min.pct=.25,logfc.threshold=.25,return.thresh = 0.05)

a <- all.markers %>% 
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_logFC)

save(all.markers, file = "data/diffgenes.RData")
write_csv(all.markers, "data/diffgenes.csv")

#compare miroglia and cams
micr.markers <- FindMarkers(all, "1", "0")

save(micr.markers, file = "data/comparison_micr_cams.RData")
write.csv(micr.markers, "data/comparison_micr_cams.csv")
#write_csv(data.frame(ID=names(all@active.ident)[all@active.ident==14]), "data/nkt-cells.csv")
