library(vegan)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(gplots)

vf<-read.table('data/90vf_count.txt',head=T, sep='\t', 
               stringsAsFactors = FALSE,check.names = F,row.names = 1)
head(vf)
dim(vf)
head(meta)
vf107 <- vf %>%
  select(any_of(meta$label))

sample<-read.table('/Users/huiluocao/Desktop/HKU-microbiology/Lactococcus_garvieve/data/sample_species.csv',head=T, sep='\t', 
                   stringsAsFactors = FALSE,check.names = F)
head(sample)

vf14_list<-read.table('/Users/huiluocao/Desktop/HKU-microbiology/Lactococcus_garvieve/data/vf15.txt',head=F, sep=',', 
                      stringsAsFactors = FALSE,check.names = F)
dim(vf14_list)
vf14 <- vf[rownames(vf) %in% vf14_list$V1, ]  
vf14<-vf14[, which(colSums(vf14) != 0)]
dim(vf14)

dist.mat <- vegdist(t(vf107),na.rm = T)
cmds <- cmdscale(dist.mat, k=3, eig=TRUE)
eigen <- cmds$eig / sum(cmds$eig) * 100

fit= envfit(cmds, t(vf107))
goodness_of_fit <- fit$vectors$r
top_variables <- names(sort(goodness_of_fit, decreasing = TRUE))[1:16]

spp.scrs <- as.data.frame(scores(fit, display = "vectors"))
spp.scrs <- cbind(spp.scrs, VFs = rownames(spp.scrs))

spp.scrs_top10<-spp.scrs[spp.scrs$VFs %in% top_variables,]

dat.merged <- (merge(cmds$points, meta, by.x=0, by.y="tip", all.x=TRUE))
head(dat.merged)
dat.merged$sample <- dat.merged$Row.names

#dat.merged$region <- relevel(factor(dat.merged$region), 'HKHuman')
adonis_result_species<-adonis2(t(vf14) ~ species+source, data = sample[,-(4)],permutations = 999,method="bray")
adonis_result_source<-adonis2(t(vf14) ~ source, data = sample[,-c(2,4)],permutations = 999,method="bray")
species_source_adonis <- paste0("adonis R2 by species: ",round(adonis_result_species$R2,4), "; P-value: ", adonis_result_species$`Pr(>F)`,";","adonis R2 by source: ",round(adonis_result_source$R2,4), "; P-value: ", adonis_result_source$`Pr(>F)`)

dat.merged<-column_to_rownames(dat.merged,var = 'Row.names')
dat.merged <- dat.merged %>% mutate(coordinates = paste(V1, V2, sep = ","))
vf_counts<-table(dat.merged$coordinates)%>%
  as.data.frame()
colnames(vf_counts)<-c('coordinates','freq')
dat.merged_new <- (merge(dat.merged, vf_counts, by.x='coordinates', by.y="coordinates", all.x=TRUE))
dat.merged_new<-column_to_rownames(dat.merged_new,var = 'Row.names')
dat.merged_new <- dat.merged_new %>% mutate(log2freq = (log2(freq)+1))
head(dat.merged_new)

ggplot(data = dat.merged_new, aes(V1, V2)) + 
  geom_density_2d(aes(x=V1, y=V2), inherit.aes = FALSE, col='grey', lwd=1) +
  geom_point(aes(col=Source, shape=Contenient, size=log2freq)) +
  geom_text_repel(
    data = subset(dat.merged_new, grepl("^Pse", sample)),
    aes(label = sample),
    size = 4,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.3,
    max.overlaps = 20,
    force = 3
  ) +
  labs(x="MDS1", y="MDS2") +
  scale_size(range = c(4, 8), name='log2(count)+1') +
  theme_bw() +
  theme(text = element_text(size = 14, family = "Times"),
        legend.position = "right")

ggplot(data = dat.merged_new, aes(V1, V2)) + 
  geom_density_2d(aes(x=V1, y=V2), inherit.aes = FALSE, col='grey', lwd=1) +
  geom_point(aes(col=Source,shape=Continent,size=log2freq))+
  labs(x=paste0('MDS1 (',round(eigen[1], 1),'%)'),
       y=paste0('MDS2 (',round(eigen[2], 1),'%)')) +
  #scale_shape_discrete(name = "Species",labels = c("L.formosensis", "L.garvieae","L.petauri"))+
  scale_size(range = c(4, 8),name='log2(count)+1')+
  #scale_color_manual(name="Source",values = c("Human infection"="#FF0000","Human fecal"="#FCC5C0","Others"="#BDBDBD","Bovines"="#84E7F7","Rainbow trout"="#74C476"))+
  #geom_segment(data = spp.scrs_top10,
  #             aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
  #             arrow = arrow(length = unit(0.1, "cm")), colour = "grey")+
  #geom_text(data=spp.scrs_top10,aes(x=Dim1+0.02,y=Dim2-0.02,label=VFs),size=5,family = "Times")+
  #geom_text(aes(label=samples, y=V2+0.01,x=V1+0.01, vjust=0),size=3)+
  geom_text(
    data = subset(dat.merged_new, grepl("^Pse", samples)),
    aes(label = samples),
    size = 4,
    fontface = "bold",
    hjust = -0.1,
    vjust = -0.5,
    colour = "red"
  )+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold",family = "Times"),
        text = element_text(size = 14,family = "Times"),
        axis.title = element_text(face="bold",family = "Times"),
        axis.text.x=element_text(size = 14,family = "Times",angle = 45,hjust = 1.0),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        #axis.title.x=element_blank(),
        legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5)))

#geom_text(data=vec.sp.df,aes(x=MDS1,y=MDS2,label=species),size=5)

#### PCOA with permutation test
otu.distance <- vegdist(t(vf), method = 'bray')
otu.distance

pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
#p <- ggplot(pc12,aes(x=V1, y=V2))+
#  geom_point(size=3)+theme_bw()
#p

df <- merge(pc12,group,by="samples")

dat.merged <- (merge(pc12, sample, by.x='samples', by.y="strain", all.x=TRUE))
head(dat.merged)

write.csv(dat.merged,file = "df.csv",row.names = TRUE)

set.seed(2)
adonis_result<-adonis2(t(vf) ~ species, data = sample[,-3],permutations = 999,method="bray")
adonis_result
#color=c('#367EB7','#E3181A')
color = c("#771155","#114477", "#DDAA77")

dune_adonis <- paste0("adonis R2: ",round(adonis_result$R2,4), "; P-value: ", adonis_result$`Pr(>F)`)
pcoa<-ggplot(data=dat.merged,aes(x=V1,y=V2,
                                 color=species))+
  theme_bw()+
  geom_point(size=8)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_text(aes(label=samples, y=V2+0.01,x=V1+0.01, vjust=0),size=3)+
  #geom_text(aes(label=samples),size=3)+
  #guides(color=guide_legend(title=NULL))+
  labs(x=paste0("PCoA1 ","(",pc[1],"%)"),
       y=paste0("PCoA2 ","(",pc[2],"%)"),
       title=dune_adonis)+
  scale_color_manual(values = color) +
  scale_fill_manual(values = color)+
  theme(axis.title.x=element_text(size=16,family = "Times"),
        axis.title.y=element_text(size=16,angle=90,family = "Times"),
        axis.text.y=element_text(size=16,family = "Times"),
        axis.text.x=element_text(size=16,family = "Times"),
        plot.title = element_text(size = 14, face = "bold",family = "Times"),
        legend.position = "bottom",
        text=element_text(family="Times",size=16),
        panel.grid=element_blank()) + 
  stat_ellipse(data=dat.merged,geom = "polygon",level=0.95,linetype = 2, linewidth=0.2,aes(fill=species),alpha=0.2,show.legend = FALSE)
pcoa

library(cowplot)
plot_grid(pnmds,pcoa,label_size = 20,rel_widths=c(10,10),labels = c('(A)', '(B)'),label_fontfamily = 'Times')

ggsave(pcoa,filename = "PCoA.pdf",width = 10, height = 10)


### plot upset
head(vf)
t(vf)
vf_all <- (merge(t(vf), sample, by.x=0, by.y="strain", all.x=TRUE))
head(vf_all)
result <- vf_all %>%
  group_by(species) %>%
  summarise_all(sum)

pheatmap::pheatmap(vf)

