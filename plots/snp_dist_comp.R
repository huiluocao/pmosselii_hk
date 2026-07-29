library(ggplot2)
library(ggsci)
library(scales)
library(data.table)
colors=pal_npg('nrc')(10) ###library(scales)
show_col(colors)
pal_simpsons("springfield")(16)
show_col(pal_simpsons("springfield")(16))

### https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html#npg
### https://rdrr.io/cran/ape/man/dist.dna.html

### snp_dist
d <- read.table("./data/phylogeny/hku_pm107_parsnp_dist.tsv", head = TRUE, sep = "\t",
              row.names = 1,check.names = F)

ind <- which(upper.tri(d, diag = TRUE), arr.ind = TRUE)
nn <- dimnames(d)
snp_diff<-data.frame(row = nn[[1]][ind[, 1]],
                     col = nn[[2]][ind[, 2]],
                     val = d[ind])
colnames(snp_diff) <- c("G1", "G2","snp")
snp_diff$G1 <- gsub("\\.fna\\.ref$", "", snp_diff$G1)
snp_diff$G2 <- gsub("\\.fna\\.ref$", "", snp_diff$G2)
head(snp_diff)
write.table(snp_diff,'./data/phylogeny/hku_pm107_parsnp_snp_diff1.txt', quote=F, sep='\t', row.names = F, col.names = T)

clusters_genome <- read.table("../../../Bacteroides/data/Bf2/snp_diff/cluster1.txt", head = TRUE, sep = "\t", 
                              check.names = F)

head(meta)
snp_diff_G1<-merge(snp_diff,meta[c('label','BAPS1')],by.x='G1',by.y='label')
dim(snp_diff)
dim(snp_diff_G1)
names(snp_diff_G1)[4]<-"G1_C"
snp_diff_G1<-merge(snp_diff_G1,meta[c('label','BAPS1')],by.x='G2',by.y='label')
names(snp_diff_G1)[5]<-"G2_C"
head(snp_diff_G1)
snp_diff_G1 <- snp_diff_G1 %>%
  mutate(Group = ifelse(G1_C == G2_C, "within", "between"))

snp_group <- read.table("../../../Bacteroides/data/Bf2/snp_diff/snp_group1.txt", head = TRUE, sep = "\t", 
                        check.names = F)
snp_diff_G1$Group <- as.factor(snp_diff_G1$Group)
# Basic violin plot
p <- ggplot(snp_diff_G1, aes(x=Group, y=snp,fill=Group)) + 
  geom_violin(alpha=0.3)+
  theme_bw()+
  coord_flip() + 
  theme(panel.border = element_rect(colour="black", size=1),
        plot.title = element_text(hjust = 0.5,face='bold',family="Times",size = 24),
        axis.text.x = element_text(face='bold',family="Times",size = 22),
        axis.text.y=element_text(face='bold',family="Times",size = 22),
        axis.title = element_text(family="Times",face="bold",size = 24),
        legend.title = element_text(colour="black", family="Times",size=24, face="bold"),
        legend.text = element_text(colour="black", family="Times",size = 22))+ #axis.ticks.x = element_blank()
  xlab("BAPS cluster level1")+
  ggtitle("Pseudomonas mosselii")+
  scale_y_continuous(name = "SNP distance")+
  labs(fill = "Comparison") ###p+scale_x_discrete(labels=c("B" = "Between", "W" = "Within"))
p+geom_boxplot(width=0.1,alpha=0.1)
###p+scale_fill_manual(values=c("#E69F00", "#56B4E9"))
p
median(subset(snp_diff_G1,Group=="between")$snp)
median(subset(snp_diff_G1,Group=="within")$snp)

##### ani distance
d2 <- read.table("./data/ani/109fna_ANIm/ANIm_percentage_identity.tab", head = TRUE, sep = "\t",
                row.names = 1,check.names = F)

ind <- which(upper.tri(d2, diag = TRUE), arr.ind = TRUE)
nn <- dimnames(d2)
ani_diff<-data.frame(row = nn[[1]][ind[, 1]],
                     col = nn[[2]][ind[, 2]],
                     val = d2[ind])
colnames(ani_diff) <- c("G1", "G2","ani")
#snp_diff$G1 <- gsub("\\.fna\\.ref$", "", snp_diff$G1)
#snp_diff$G2 <- gsub("\\.fna\\.ref$", "", snp_diff$G2)
head(ani_diff)
write.table(snp_diff,'./data/ani/hku_pm107_ani_diff1.txt', quote=F, sep='\t', row.names = F, col.names = T)

clusters_genome <- read.table("../../../Bacteroides/data/Bf2/snp_diff/cluster1.txt", head = TRUE, sep = "\t", 
                              check.names = F)

head(meta)
ani_diff_G1<-merge(ani_diff,meta[c('tip','BAPS_level1')],by.x='G1',by.y='tip')
dim(ani_diff)
dim(ani_diff_G1)
head(ani_diff_G1)
names(ani_diff_G1)[4]<-"G1_C"
ani_diff_G1<-merge(ani_diff_G1,meta[c('tip','BAPS_level1')],by.x='G2',by.y='tip')
names(ani_diff_G1)[5]<-"G2_C"
head(ani_diff_G1)
ani_diff_G1 <- ani_diff_G1 %>%
  mutate(Group = ifelse(G1_C == G2_C, "within", "between"))

snp_group <- read.table("../../../Bacteroides/data/Bf2/snp_diff/snp_group1.txt", head = TRUE, sep = "\t", 
                        check.names = F)
ani_diff_G1$Group <- as.factor(ani_diff_G1$Group)
# Basic violin plot
p <- ggplot(ani_diff_G1, aes(x=Group, y=ani,fill=Group)) + 
  geom_violin(alpha=0.3)+
  theme_bw()+
  theme(panel.border = element_rect(colour="black", size=1),
        plot.title = element_text(hjust = 0.5,face='bold',family="Times",size = 24),
        axis.text.x = element_text(face='bold',family="Times",size = 22),
        axis.text.y=element_text(face='bold',family="Times",size = 22),
        axis.title = element_text(family="Times",face="bold",size = 24),
        legend.title = element_text(colour="black", family="Times",size=24, face="bold"),
        legend.text = element_text(colour="black", family="Times",size = 22))+ #axis.ticks.x = element_blank()
  xlab("BAPS cluster level1")+
  ggtitle("Pseudomonas mosselii")+
  scale_y_continuous(name = "ani distance")+
  labs(fill = "Comparison") ###p+scale_x_discrete(labels=c("B" = "Between", "W" = "Within"))
p+geom_boxplot(width=0.1,alpha=0.1)
###p+scale_fill_manual(values=c("#E69F00", "#56B4E9"))
p
median(subset(snp_diff_G1,Group=="between")$snp)
median(subset(snp_diff_G1,Group=="within")$snp)



ggplot(ani_diff_G1, aes(x = ani, fill = Group, color = Group)) +
  geom_density(alpha = 0.35, adjust = 1.2, size = 1.2) +
  
  # Specific colors
  scale_fill_manual(values = c("within" = "#E41A1C", "between" = "#377EB8")) +
  scale_color_manual(values = c("within" = "#E41A1C", "between" = "#377EB8")) +
  
  # To add grid lines (like the image)
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", size = 0.5)
  ) +
  
  # Set y-axis limit similar to image
  scale_y_continuous(limits = c(0, 350), expand = c(0, 0))
