library(tidyverse)
library(viridis)
library(readxl)
library(RaceID)
library(RColorBrewer)
library(pheatmap)
library(fishualize)
library(assertthat)
library(Seurat)
library(ggpubr)
library(ggrepel)

#export sessions info
#writeLines(capture.output(sessionInfo()), "data/sessionInfo.txt")

date = Sys.Date()
#load('data/harmony-integration-SCTransform.RData')
load("data/all.RData")
source("bin/functions.R")

#build data frame
if (!file.exists("data/metadata-myeloid-cells.RData")) {
  metadata <- all@meta.data
  retain_cl <- levels(all)
  save(retain_cl, file = "data/retain_cl-myeloid-cells.RData")
  ord_clust <- levels(all)
  save(ord_clust, file = "data/ord_clust-myeloid-cells.RData")
  
  save(metadata, file = "data/metadata-myeloid-cells.RData")
} else {
  load("data/metadata-myeloid-cells.RData")
  load("data/retain_cl-myeloid-cells.RData")
  load('data/ord_clust-myeloid-cells.RData')
}


#Plot cluster tsne map
       tsne <- DimPlot(all, label = T, cols = c(colors_pat, colors_many, colors_fig)[-2], pt.size = 6) +
          theme_void()

        tsne
        
        ggsave(paste0('plots/umap/', date, '-clusters-umap.pdf'), useDingbats=F, width = 8.57, height = 5.79)  
        
        #no legend
        #tsne <- umap_plot_seurat(data = metadata, FILL = metadata$Cluster, fill_colors =   c(colors_many), point_size = 1.5) +
         # theme(legend.position = "None")
        tsne <- DimPlot(all, label = T, cols = c(colors_pat, colors_many, colors_fig)[-2], pt.size = 6) +
          theme_void() +
          NoLegend() 
        
        tsne
        
        ggsave(paste0('plots/umap/', date, '-clusters-umap-no-legend.pdf'), useDingbats=F, width = 8.57, height = 5.79)  
        
  
#plot heatmap with markers
        #diffgenes
        if (!file.exists("data/diffgenes-clusters-myeloid-cells.csv")) {
        
          all.markers<-FindAllMarkers(all,only.pos=F,min.pct=.2,logfc.threshold=.2,return.thresh = 0.05) #,only.pos=T,min.pct=.25,logfc.threshold=.25,return.thresh = 0.05
         
          save(all.markers, file = "data/diffgenes-clusters-myeloid-cells.RData")
          write_csv(all.markers, "data/diffgenes-clusters-myeloid-cells.csv")
     
          
          #all.markers.roc <- FindAllMarkers(all,test.use = "roc", only.pos = T)
          #all.markers.bimod <- FindAllMarkers(all,test.use = "bimod", only.pos = T)
          #all.markers.t <- FindAllMarkers(all,test.use = "t", only.pos = T)
          
        
          
        } else {
        all.markers <- read_csv( "data/diffgenes-clusters-myeloid-cells.csv")
        }
        #levels(all) <- as.character(ord_clust)
        #levels(all@active.ident) <- paste0("C", 1:length(ord_clust))
        
        all.markers <- all.markers[!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", all.markers$gene),]
        #all.markers$cluster <- factor(all.markers$cluster, levels = ord_clust)
        
        top20 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) #filter(!cluster %in% c("5","8")) %>% 
        
        heat <- DoHeatmap(all,features = top20$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2]) #[,names(all@active.ident)[!all@active.ident %in% c("5", "8")]]
        heat + #+ NoLegend()
        #scale_fill_fish(option = "Lepomis_megalotis") #option = "Hypsypops_rubicundus"
        scale_fill_viridis(option = "A")
          
        ggsave(paste0("plots/heatmaps/",date,"top10-gene-heatmap-viridis-A-myeloid-cells.pdf"), width = 12, height = 8)
        
        postscript(paste0("plots/heatmaps/",date,"top15-gene-heatmap-viridis-A-myeloid-cells.ps"), width = 30, height = 20)
        heat + #+ NoLegend()
          #scale_fill_fish(option = "Lepomis_megalotis") #option = "Hypsypops_rubicundus"
          scale_fill_viridis(option = "A")
        dev.off()
        
        heat <- DoHeatmap(all,features = top20$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2]) #[,names(all@active.ident)[!all@active.ident %in% c("5", "8")]]
        DoHeatmap(all,features = top20$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2], group.bar = F)+
          theme_void() +
          scale_fill_viridis(option = "A")  + 
          NoLegend() 
        
        ggsave(paste0("plots/heatmaps/",date,"top10-gene-heatmap-viridis-A-myeloid-cells-heatmap-only.png"), width = 12, height = 8)
        
        
        
        heat
        ggsave(paste0("plots/heatmaps/",date,"top20-gene-heatmap-original-color-myeloid-cells.pdf"), width = 30, height = 20)
       
        #text heatmap color scales
        my_cols <- toupper(rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')))
        heat + scale_fill_gradientn(colors=my_cols)
        
        ggsave(paste0("plots/heatmaps/",date,"top20-gene-heatmap-blue-to-red-color.pdf"), width = 10, height = 8)
        ggsave(paste0("plots/heatmaps/",date,"top20-gene-heatmap-oblue-to-red-color.png"), width = 12, height = 8)
        
        #plot cell signatures
        signature_genes <-  read_excel('/Users/romansankowski/Documents/single_cell_analysis/EAE_Final/Cluster-information-dimRed.xlsx', 'Core signature',skip = 2)
        
        #take out calm1
        signature_genes$T_cells[4] <- NA
        #signature_genes <- data.frame("ahr"=c("Ahr", "Ahrr", "Cyp1a1", "Cyp1b1", "Entpd1", "Hif1a", "Il1b", "Il27", "Klf4", "Pparg", "Stat1", "Stat3", "Tiparp", "Vegfa" ), stringsAsFactors = F)
        
        
        for (i in colnames(signature_genes)) {
          tryCatch({
            pdf(paste0('plots/umap/', date, as.character(i), '-gene_signature.pdf'), width = 8.57, height = 5.79, useDingbats = F)
            pl <- plot_expmap_seurat(features=c(na.omit(signature_genes[[i]])), point_size = 6, .retain_cl = unique(all$seurat_clusters))
            print(pl)
            dev.off()
            
            pdf(paste0('plots/umap/', date, as.character(i), '-gene_signature-logsc.pdf'), width = 8.57, height = 5.79, useDingbats = F)
            pl <- plot_expmap_seurat(features=c(na.omit(signature_genes[[i]])), point_size = 6, logsc = T,.retain_cl = unique(all$seurat_clusters))
            print(pl)
            dev.off()
            
          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
          on.exit(dev.off())
        }     
        

for (i in colnames(signature_genes)) {
  tryCatch({
    pdf(paste0('plots/umap/', date,"-" ,as.character(i), '-gene_signature.pdf'), width = 8.57, height = 5.79, useDingbats = F)
    pl <- plot_expmap_seurat(na.omit(signature_genes[[i]]), object = all,point_size = 1, reduction = "NetUMAP", .retain_cl = levels(all)) #+
    #scale_color_viridis_c(option = "A")
    #scale_color_gradientn(colors=my_cols)
    print(pl)
    dev.off()
    
    pdf(paste0('plots/umap/', date, "-", as.character(i), '-gene_signature-logsc.pdf'), width = 8.57, height = 5.79, useDingbats = F)
    pl <- plot_expmap_seurat(na.omit(signature_genes[[i]]), object = all,point_size = 1, logsc = T, reduction = "NetUMAP", .retain_cl = levels(all)) #+
    #scale_color_viridis_c(option = "A")
    #scale_color_gradientn(colors=my_cols)
    print(pl)
    dev.off()
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  on.exit(dev.off())
}     


#plot single cell gene expression
#my_cols <- toupper(c("#2166ac", "#f7f7f7", "#ca0020"))
#my_cols <- toupper(rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')))

my_cols <- c("darkblue","lightblue2","yellow","red2")
     
genes_to_plot <- c("CCR2", "S100A6", "CD1C")
#genes_to_plot <- unique(top20$gene)

for (i in genes_to_plot) {
              tryCatch({
                pdf(paste0('plots/umap/', i, '.pdf'), width = 8.57, height = 5.79, useDingbats = F)
                pl <- FeaturePlot(all, i, pt.size = 1) +
                  #scale_color_viridis_c(option = "A")
                  scale_color_gradientn(colors=my_cols) +
                  theme_void()
                print(pl)
                dev.off()
                
                
              }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
              on.exit(dev.off())
            }    

#test colors
my_cols <- toupper(c("#2166ac", "#f7f7f7", "#67001f"))
plot_expmap_seurat(i, object = all,point_size = 1.5) +
  scale_color_gradientn(colors=my_cols)
my_cols <- toupper(c("#2d004b", "#f7f7f7", "#7f3b08"))
plot_expmap_seurat(i, object = all,point_size = 1.5) +
  scale_color_gradientn(colors=my_cols)
my_cols <- toupper(c("#542788", "#f7f7f7", "#7f3b08"))
plot_expmap_seurat(i, object = all,point_size = 1.5) +
  scale_color_gradientn(colors=my_cols)
my_cols <- toupper(c("#2166ac", "#f7f7f7", "#ca0020"))
plot_expmap_seurat(i, object = all,point_size = 1.5) +
  scale_color_gradientn(colors=my_cols)

#log scaled genes
          for (i in unique(all.markers$gene)) {
            tryCatch({
              svg(paste0('plots/umap/logsc-', i, '.svg'), width = 8.57, height = 5.79)
              pl <- plot_expmap_seurat(i, object = all,point_size = 1.5, logsc = T)+
                #scale_color_viridis_c(option = "A")
                scale_color_gradientn(colors=my_cols)
              print(pl)
              dev.off()
              
            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            on.exit(dev.off())
          }     

#volcano plots
mhc_hilo <- FindMarkers(all, ident.1 = "C3",
                        ident.2 = "C2b",
                        logfc.threshold = 0.01,
                        min.pct = 0.01) %>%
  rownames_to_column(var="gene")  %>%
  mutate(p_val_adj = ifelse(p_val_adj == 0, 2.225074e-308, p_val_adj))

save(mhc_hilo, file = "data/mhc_hilo_early_diffgenes.RData")
load("data/mhc_hilo_early_diffgenes.RData")

mhc_hilo <- mhc_hilo[!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", mhc_hilo$gene),]
top5_wt <- mhc_hilo %>%
  filter(abs(.$avg_logFC) > .25 & p_val_adj < .05) %>%
  top_n(10, avg_logFC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes)

top_5_both <- mhc_hilo %>%
  filter(abs(.$avg_logFC) > .25 & p_val_adj < .05) %>%
  top_n(-10, avg_logFC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes) %>%
  bind_rows(top5_wt) 

mhc_hilo <- mhc_hilo %>%
  left_join(top_5_both) %>%
  mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_logFC) > .25, "sig.", "not sig."))


mhc_hilo_volcano <- ggplot(mhc_hilo, aes(x=avg_logFC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  geom_point(size=5) + #aes(size=avg_logFC)
  geom_text_repel(size=7, box.padding=1) +
  expand_limits(x=c(-2.25, 2.25)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. logFC", y="-log10 transf. adj. p-value")

mhc_hilo_volcano 

ggsave("plots/others/volcano_plot_early_mhc_hilo.pdf", width = 7, height = 7, useDingbats=F)

#mhc_hi_ccr2
        mhc_hi_ccr2 <- FindMarkers(all, ident.1 = "C2b",
                                   ident.2 = "C2a",
                                   logfc.threshold = 0.01,
                                   min.pct = 0.01) %>%
          rownames_to_column(var="gene")  %>%
          mutate(p_val_adj = ifelse(p_val_adj == 0, 2.225074e-308, p_val_adj))
        
        save(mhc_hi_ccr2, file = "data/mhc_hi_ccr2_diffgenes.RData")
        load("data/mhc_hi_ccr2_diffgenes.RData")
        
        mhc_hi_ccr2 <- mhc_hi_ccr2[!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", mhc_hi_ccr2$gene),]
        top5_wt <- mhc_hi_ccr2 %>%
          filter(abs(.$avg_logFC) > .25 & p_val_adj < .05) %>%
          top_n(15, avg_logFC) %>%
          mutate(show_genes = gene) %>%
          select(gene, show_genes)
        
        top_5_both <- mhc_hi_ccr2 %>%
          filter(abs(.$avg_logFC) > .25 & p_val_adj < .05) %>%
          top_n(-15, avg_logFC) %>%
          mutate(show_genes = gene) %>%
          select(gene, show_genes) %>%
          bind_rows(top5_wt) 
        
        mhc_hi_ccr2 <- mhc_hi_ccr2 %>%
          left_join(top_5_both) %>%
          mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_logFC) > .25, "sig.", "not sig."))
        
        
        mhc_hi_ccr2_volcano <- ggplot(mhc_hi_ccr2, aes(x=avg_logFC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
          geom_point(size=5) + #aes(size=avg_logFC)
          geom_text_repel(size=7, box.padding=1) +
          expand_limits(x=c(-2.25, 2.25)) +
          theme_bw() +
          theme(panel.grid = element_blank(),
                text = element_text(size=25)) +
          scale_color_manual(values = c("light grey", "black")) +
          NoLegend() +
          labs(x="avg. logFC", y="-log10 transf. adj. p-value")
        
        mhc_hi_ccr2_volcano 
        
        ggsave("plots/others/volcano_plot_early_mhc_hi_ccr2.pdf", width = 7, height = 7, useDingbats=F)

#mhc_hi_monos
        mhc_hi_monos <- FindMarkers(all, ident.1 = "C2b",
                                    ident.2 = "C1",
                                    logfc.threshold = 0.01,
                                    min.pct = 0.01) %>%
          rownames_to_column(var="gene")  %>%
          mutate(p_val_adj = ifelse(p_val_adj == 0, 2.225074e-308, p_val_adj))
        
        save(mhc_hi_monos, file = "data/mhc_hi_monos_diffgenes.RData")
        load("data/mhc_hi_monos_diffgenes.RData")
        
        mhc_hi_monos <- mhc_hi_monos[!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", mhc_hi_monos$gene),]
        top5_wt <- mhc_hi_monos %>%
          filter(abs(.$avg_logFC) > .25 & p_val_adj < .05) %>%
          top_n(15, avg_logFC) %>%
          mutate(show_genes = gene) %>%
          select(gene, show_genes)
        
        top_5_both <- mhc_hi_monos %>%
          filter(abs(.$avg_logFC) > .25 & p_val_adj < .05) %>%
          top_n(-15, avg_logFC) %>%
          mutate(show_genes = gene) %>%
          select(gene, show_genes) %>%
          bind_rows(top5_wt) 
        
        mhc_hi_monos <- mhc_hi_monos %>%
          left_join(top_5_both) %>%
          mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_logFC) > .25, "sig.", "not sig."))
        
        
        mhc_hi_monos_volcano <- ggplot(mhc_hi_monos, aes(x=avg_logFC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
          geom_point(size=5) + #aes(size=avg_logFC)
          geom_text_repel(size=7, box.padding=1) +
          expand_limits(x=c(-2.25, 2.25)) +
          theme_bw() +
          theme(panel.grid = element_blank(),
                text = element_text(size=25)) +
          scale_color_manual(values = c("light grey", "black")) +
          NoLegend() +
          labs(x="avg. logFC", y="-log10 transf. adj. p-value")
        
        mhc_hi_monos_volcano 
        
        ggsave("plots/others/volcano_plot_early_mhc_hi_monos.pdf", width = 7, height = 7, useDingbats=F)
        
#cross species comparison
        load("data/mhc_hilo_early_diffgenes.RData")
        mhc_hilo_hum <- mhc_hilo 
        load("/home/roman/Documents/Single cell analysis/eae-raceid4/data/mhc_hilo_early_diffgenes.RData")

        mhc_hilo_hum <- mhc_hilo_hum %>%
          mutate(sign_log_pval_hum = (avg_logFC / abs(avg_logFC)) * (-log(p_val_adj)),
                 species_2 = "Human") %>%
          dplyr::select(gene, species_2, sign_log_pval_hum)

        mhc_hilo <- mhc_hilo %>%
          mutate(sign_log_pval_ms = (avg_logFC / abs(avg_logFC)) * (-log(p_val_adj)),
                 gene = toupper(gene),
                 species_1 = "Mouse") %>%
          dplyr::select(gene, species_1, sign_log_pval_ms)
        
         mhc_hilo_both <- mhc_hilo %>%
           full_join(mhc_hilo_hum) %>% 
           mutate(gene_sig = ifelse(sign_log_pval_ms > 1 | sign_log_pval_hum > 1, gene, NA))
           
         ggplot(mhc_hilo_both, aes(sign_log_pval_ms, sign_log_pval_hum, label=gene_sig)) +
           geom_point(pch = 21) +
           scale_x_log10() +
           scale_y_log10() +
           ggrepel::geom_text_repel()
        
#Cumulative cell types per cluster
data.frame(table(metadata$Cluster, metadata$Compartment)) %>%
  group_by(Var2) %>%
  mutate(cum_sum = cumsum(Freq)/sum(Freq)) %>%
  ggplot(aes(factor(Var1),cum_sum, color = Var2, group=Var2)) +
  theme_minimal() +
  geom_step(size=2)  +
  theme(legend.text=element_text(size = 25),
        axis.text = element_text(size=20),
        axis.title = element_text(size = 20)) +
  scale_color_brewer("",palette = "Set1") +
  labs(x="Cluster", y = "Cum. Density")#+
  #scale_y_log10() +
  #facet_wrap(~Var2)

ggsave("plots/others/ecdf-compartments.pdf")
ggsave("plots/others/ecdf-compartments.png")


cl_dist <- data.frame(table(metadata$Cluster, metadata$Compartment)) %>%
  group_by(Var2) %>%
  mutate(cum_sum = cumsum(Freq)/sum(Freq)) %>%
  select(-Freq) %>%
  pivot_wider(names_from = Var2, values_from = cum_sum) %>%
  mutate(Var1 =paste0("C", Var1)) %>%
  as.data.frame()

rownames(cl_dist) <- cl_dist$Var1
cl_dist <- cl_dist[,-1]

plot(hclust(dist(as.matrix(t(cl_dist))), method = "single"), xlab = " ")

pdf("plots/others/dendrogram-compartments.pdf")
plot(hclust(dist(as.matrix(t(cl_dist))), method = "single"), xlab = " ")
dev.off()


png("plots/others/dendrogram-compartments.png")
plot(hclust(dist(as.matrix(t(cl_dist))), method = "single"), xlab = " ")
dev.off()

#stat testing
cl_dist <- data.frame(table(metadata$Cluster, metadata$Compartment)) %>%
  group_by(Var2) %>%
  mutate(cum_sum = cumsum(Freq)/sum(Freq))
ks.test(cl_dist$cum_sum[cl_dist$Var2=="CP"], cl_dist$cum_sum[cl_dist$Var2=="Micr"])
ks.test(cl_dist$cum_sum[cl_dist$Var2=="PVM"], cl_dist$cum_sum[cl_dist$Var2=="Micr"])
ks.test(cl_dist$cum_sum[cl_dist$Var2=="Dura"], cl_dist$cum_sum[cl_dist$Var2=="Micr"])
ks.test(cl_dist$cum_sum[cl_dist$Var2=="MenM"], cl_dist$cum_sum[cl_dist$Var2=="Micr"])

#k-sample Anderson-Darling test
library(kSamples) #url https://stats.stackexchange.com/questions/35461/is-there-a-multiple-sample-version-or-alternative-to-the-kolmogorov-smirnov-test
ad.test(cl_dist$cum_sum[cl_dist$Var2=="Micr"], cl_dist$cum_sum[cl_dist$Var2=="PVM"], cl_dist$cum_sum[cl_dist$Var2=="CP"], cl_dist$cum_sum[cl_dist$Var2=="MenM"], cl_dist$cum_sum[cl_dist$Var2=="Dura"])
ad.test(cl_dist$cum_sum[cl_dist$Var2=="PVM"], cl_dist$cum_sum[cl_dist$Var2=="CP"], cl_dist$cum_sum[cl_dist$Var2=="MenM"], cl_dist$cum_sum[cl_dist$Var2=="Dura"])
ad.test(cl_dist$cum_sum[cl_dist$Var2=="CP"], cl_dist$cum_sum[cl_dist$Var2=="MenM"], cl_dist$cum_sum[cl_dist$Var2=="Dura"])
ad.test(cl_dist$cum_sum[cl_dist$Var2=="Micr"], cl_dist$cum_sum[cl_dist$Var2=="PVM"])

#gene violin plot

df <- all.markers %>%
  filter(cluster %in% c(6,4,11,1,0,9,3,5)) %>%
  group_by(cluster) %>%
  top_n(10, wt=avg_logFC) 

genes <- c("MRC1", "MS4A7", "CD163", "LYVE1", "APOE", 'P2RY12', 'CX3CR1', 'CSF1R', 'TMEM119', 'SLC2A5', "S100A6", "CD74", "HLA-DRA", "SIGLEC1", "CD1C", "SPP1", "CCL2")
#df2 <- all@assays$RNA@scale.data[rownames(all@assays$RNA@scale.data) %in% c(unique(df$gene),"S100A6"), names(all@active.ident)[all@active.ident %in% c(6,4,11,1,0,9,3,5)]] 
df2 <- all@assays$RNA@scale.data[rownames(all@assays$RNA@scale.data) %in% genes, names(all@active.ident)[all@active.ident %in% c(6,4,11,1,0,9,3,5)]] 
df2 <- as.data.frame(as.matrix(t(df2)))
gene_order <- hclust(dist(t(scale(df2))))
gene_order <- colnames(df2)[gene_order$order]

df2$Cluster <- factor(all@active.ident[all@active.ident %in% c(6,4,11,1,0,9,3,5)], levels = c(6,4,11,1,0,9,3,5)) 

df2 <- df2 %>%
  #pivot_longer(A2M:VIM) %>% filter(name != "RND3")
  pivot_longer(APOE:TMEM119) 

df2$name <- factor(df2$name, levels = gene_order)


df2 %>%
  ggplot(aes(name, value, fill=Cluster)) +
  geom_violin(scale = "width", lwd=0.5) +
  facet_grid(facets = ~Cluster,
             drop = TRUE, 
             #space = "free", 
             #scales = "free_x", 
             #switch = "x",
             space = "free_x") +
  coord_flip() +
  theme_pubr() +
  labs(y="Expression (A.U.)",x=element_blank()) +
  #scale_fill_brewer(palette = "Set1", guide=F) + 
  scale_fill_manual(values = colors_pat, guide=F) +
  theme(
    strip.text.x = element_text(
      size = 12, color = "black", face="bold"
    ))

ggsave("plots/others/gene-violin-plots.pdf", height = 5, width = 10)
ggsave("plots/others/gene-violin-plots.png", height = 5, width = 10)


#plot quadprog
load("data/2020-02-27qp_micr_vs_cams_weights_based_on_top150cells_expressing_cell_type_signatures.RData")

pvM <- as.data.frame(pvM) %>%
  rownames_to_column(var="ID")
df2 <- left_join(metadata, pvM)    

plot_continuous("mhclo_cam", na.omit(df2), point_size = 1)
plot_continuous("mhchi_cam", na.omit(df2), point_size = 1)
plot_continuous("micr", df2, point_size = 1)
plot_continuous("cams", df2, point_size = 1)

plot_continuous("mhc_lo_micr_c3", na.omit(df2), point_size = 1)
plot_continuous("mhc_hi_micr_c2", na.omit(df2), point_size = 1)
plot_continuous("proinfl_micr_c1589", na.omit(df2), point_size = 1)
plot_continuous("aging_micr_c6", na.omit(df2), point_size = 1)

#mhclo
        plot_continuous("cams", df2, point_size = 1)
        ggsave("plots/umap/quadprog_alignment_cams.pdf", width = 8.48, height = 5.76)
        
        svg("plots/umap/quadprog_alignment_cams.svg", width = 8.48, height = 5.76)
        plot_continuous("micr", df2, point_size = 1)
        dev.off()

#mhchi
        plot_continuous("micr", df2, point_size = 1)
        ggsave("plots/umap/quadprog_alignment_micr.pdf", width = 8.48, height = 5.76)
        
        svg("plots/umap/quadprog_alignment_micr.svg", width = 8.48, height = 5.76)
        plot_continuous("micr", df2, point_size = 1)
        dev.off()

#mhclo
        plot_continuous("micr", df2, point_size = 1)
        ggsave("plots/umap/quadprog_alignment_micr.pdf", width = 8.48, height = 5.76)
        
        svg("plots/umap/quadprog_alignment_micr.svg", width = 8.48, height = 5.76)
        plot_continuous("micr", df2, point_size = 1)
        dev.off()

load("data/quad_programming_microglia_weights_based_on_nat-neurosci_dataset.RData")
        
        pvM <- as.data.frame(pvM) %>%
          rownames_to_column(var="ID")
        df2 <- left_join(metadata, pvM) %>%
          na.omit()
        
        
        plot_continuous("mhc_lo_micr_c3", df2, point_size = 1)
        plot_continuous("mhc_hi_micr_c2", df2, point_size = 1)
        plot_continuous("proinfl_micr_c1589", df2, point_size = 1)
        plot_continuous("aging_micr_c6", df2, point_size = 1)
        
        #mhclo
        plot_continuous("mhc_lo_micr_c3", df2, point_size = 1)
        ggsave("plots/umap/quadprog_alignment_microglia_mhc_lo_micr_c3.pdf", width = 8.48, height = 5.76)
        
        svg("plots/umap/quadprog_alignment_microglia_mhc_lo_micr_c3.svg", width = 8.48, height = 5.76)
        plot_continuous("mhc_lo_micr_c3", df2, point_size = 1)
        dev.off()
        
        #mhchi
        plot_continuous("mhc_hi_micr_c2", df2, point_size = 1)
        ggsave("plots/umap/quadprog_alignment_microglia_mhc_hi_micr_c2.pdf", width = 8.48, height = 5.76)
        
        svg("plots/umap/quadprog_alignment_microglia_mhc_hi_micr_c2.svg", width = 8.48, height = 5.76)
        plot_continuous("mhc_hi_micr_c2", df2, point_size = 1)
        dev.off()
        
        #inflamm
        plot_continuous("proinfl_micr_c1589", df2, point_size = 1)
        ggsave("plots/umap/quadprog_alignment_microglia_proinfl_micr_c1589.pdf", width = 8.48, height = 5.76)
        
        svg("plots/umap/quadprog_alignment_microglia_proinfl_micr_c1589.svg", width = 8.48, height = 5.76)
        plot_continuous("proinfl_micr_c1589", df2, point_size = 1)
        dev.off()
        
        #aging
        plot_continuous("aging_micr_c6", df2, point_size = 1)
        ggsave("plots/umap/quadprog_alignment_microglia_aging_micr_c6.pdf", width = 8.48, height = 5.76)
        
        svg("plots/umap/quadprog_alignment_microglia_aging_micr_c6.svg", width = 8.48, height = 5.76)
        plot_continuous("aging_micr_c6", df2, point_size = 1)
        dev.off()
        
#line plot

for (i in c("WT", "T138", "T143")) {
  
  df3 <- df2 %>%
    arrange(Cluster,desc(df2[[i]]))
  
  df3$rowID <- factor(1:nrow(df3))
  
  df3$Cluster <- reorder(df3$Cluster, desc(df3[[i]]), FUN = median)
  
  line_plot <- ggplot(df3, aes(x=rowID, y=df3[[i]], color = Cluster, fill = Cluster, group=Cluster)) +
    geom_bar(stat = 'identity') + #width = 0.1, 
    facet_grid(facets = ~Cluster, 
               drop = TRUE, 
               #space = "free", 
               scales = "free_x", 
               switch = "x",
               space = "free_x") +
    labs(title = i, y = 'Gene Expression', x = 'Cluster') +
    theme_minimal() +
    theme(axis.line = element_blank(), 
          #axis.title.y = element_blank(),
          #axis.title.x = element_blank(),
          axis.ticks.y = element_blank(), 
          strip.text.x = element_text(),
          #axis.text.y = element_text(size = 10), 
          axis.text.x = element_blank(),#element_text(size = cex.col), 
          strip.background = element_blank(), 
          panel.grid = element_blank(),
          panel.spacing.x = unit(0.2, units = 'line'),
          legend.position = 'None') +
    tidytext::scale_x_reordered() +
    geom_boxplot(fill="white", color="black", outlier.shape = NA) +
    scale_fill_manual(values =c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3')) +
    scale_color_manual(values = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3')) +
    expand_limits(y=c(0,1))
  print(line_plot)
  ggsave(paste0("plots/others/",i,"-line-boxplot.pdf"), width = 8.48, height = 5.76)
  svg(paste0("plots/others/",i,"-line-boxplot.svg"), width = 8.48, height = 5.76)
  print(line_plot)
  dev.off()
}


#donut plot
          metadata %>%
            group_by(Diagnosis, Cluster) %>%
            summarise(freq = length(Cluster)) %>%
            write_csv("cluster-diagnosis-quant.csv") %>%
            ggplot(aes(x=2, y=freq,fill=Cluster)) +
            geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
            coord_polar(theta='y', start=0) +
            theme_void() +
            scale_fill_manual(values=colors_many) +
            facet_wrap(~Diagnosis) +
            xlim(0.2, 2.5)
          ggsave("plots/others/donutplots-diagnoses-clusters.pdf")
          
          metadata %>%
            group_by(Diagnosis, Cluster) %>%
            summarise(freq = length(Cluster)) %>%
            ungroup() %>%
            mutate(Diagnosis = factor(Diagnosis, levels = c("Healthy", "IDHmutGBM", "IDHwtGBM")),
                   Cluster = tidytext::reorder_within(Cluster, freq, Diagnosis)) %>%
            ggplot(aes(x=2, y=freq,fill=Cluster)) +
            geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
            coord_polar(theta='y', start=0) +
            theme_void() +
            scale_fill_manual(values=c(colors_many, colors_fig)) +
            facet_wrap(~Diagnosis) +
            xlim(0.2, 2.5) 
          
          ggsave("plots/others/ordered-donutplots-diagnoses-clusters.pdf")

