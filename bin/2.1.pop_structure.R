# Populations genomics of Quercus macdougallii
## MNS, PCA and DAPC analysis for population structure
## Graphics of fastStructure and ADMIXTURE results
## Nelly J. Pacheco Cruz
## August 2024

# Cargar bibliotecas
library(vcfR)

#Carga el archivo VCF en R
VCF_ref_gen_qrob <- read.vcfR("../data/1.3.assembly_variant_calling/denovo_trim01_1_sorted.vcf")

CGUM_79_1.VCF = VCF_ref_gen_qrob

#Cargar el archivo de metadatos
pop.metadata <- read.csv("../metadata/Qmacdougalli_79ind_.csv")
head(pop.metadata)

#Verifica que los nombres de los individuos en el VCF coincidan con los metadatos
all(colnames(VCF_ref_gen_qrob@gt)[-1] == pop.metadata$ID)

#Converting the dataset to a genlight object
ref_gen_qrob <- vcfR2genlight(VCF_ref_gen_qrob)
ref_gen_qrob_pop <- vcfR2genlight(VCF_ref_gen_qrob)
#Ploidy
ploidy(ref_gen_qrob) <- 2
ref_gen_qrob@ploidy


#Agregar el nivel de pop con los metadatos de SITE y ZONE
pop(cgum_79_1_vcf) <- pop.data$SITE
pop(cgum_79_1_vcf_pop) <- pop.data$POP
class(cgum_79_1_vcf_pop)
###Convertir de genlight a genind y a genclone (DartR)
cgum_79_1_genind <- gl2gi(cgum_79_1_vcf, v=1)
cgum_79_1_genclone <- as.genclone(cgum_79_1_genind)
#Para POP
cgum_79_1_vcf_pop_genind <- gl2gi(cgum_79_1_vcf_pop, v=1)
cgum_79_1_vcf_pop_genclone <- as.genclone(cgum_79_1_vcf_pop_genind)
#
#glfst_genind <- gl2gi(glfst, v=1)
#glfst_genclone <- as.genclone(glfst_genind)


#pop(glfst) <- pop.data$POP
#ploidy(glfst) <- 2

############################# Distance matrices
cgum_79_1.dist <- dist(cgum_79_1_vcf)
cgum_79_1.dist.x <- poppr::bitwise.dist(cgum_79_1_vcf  , percent = T, euclidean = F, mat = T)





cols <- c("#7570B3", "#075277","#00B1E8","#1FC944",
          "#E6AB02", "#E7298A","#E07E34", "#F15858", "red")


cgum_79_1.pca <- glPca(cgum_79_1_vcf, nf=80)

glPca(cgum_79_1_vcf)
sum(100*cgum_79_1.pca$eig/sum(cgum_79_1.pca$eig))
(cgum_79_1.pca$eig[1]/sum(cgum_79_1.pca$eig))*100
(cgum_79_1.pca$eig[2]/sum(cgum_79_1.pca$eig))*100
(cgum_79_1.pca$eig[3]/sum(cgum_79_1.pca$eig))*100
(cgum_79_1.pca$eig[4]/sum(cgum_79_1.pca$eig))*100
#Barplot de los eigenvalores
barplot(100*cgum_79_1.pca$eig/sum(cgum_79_1.pca$eig), main="PCA Eigenvalores")
title(ylab="Porcentaje de la varianza/explicada", line = 2)
title(xlab="Eigenvalores", line = 1)
#Reajustar los márgenes
par("mar")
par(mar=c(4,4,4,4))
#Scores
cgum_79_1.pca.scores <- as.data.frame(cgum_79_1.pca$scores)
cgum_79_1.pca.scores$pop <- pop(cgum_79_1_vcf)
cgum_79_1.pca.scores$zone <- pop(cgum_79_1_vcf_pop)

# PLOT PCA
# Definir colores
cols <- c("#00B1E8","#075277","#E7298A","#E07E34","#F15858","#E6AB02","#7570B3","#1FC944", "red")

# Asignar colores basados en las poblaciones
colores <- cols[as.numeric(cgum_79_1.pca.scores$pop)]

# Crear el plot
set.seed(12345)
p1 <- ggplot(cgum_79_1.pca.scores, aes(x=PC1, y=PC2, colour=as.factor(pop), shape=as.factor(zone), fill=as.factor(pop))) + 
  geom_point(size=4, alpha=0.7) + 
  scale_color_manual(values = cols, name="Poblaciones") +
  scale_fill_manual(values = cols, name="Poblaciones") +
  scale_shape_manual(name = "Sitios", values = c(21, 24, 25)) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 18),
        axis.title.x = element_text(size=20), 
        axis.text.x = element_text(size=16),
        axis.title.y = element_text(size=20), 
        axis.text.y = element_text(size=16)) + 
  xlab("PC1 %4.62") + 
  ylab("PC2 %2.76") + 
  geom_text(aes(label=cgum_79_1_vcf$ind.names), hjust=-0.4, vjust=0.8, size=5)

# Mostrar el plot
print(pca_ref_gen_qrob)

