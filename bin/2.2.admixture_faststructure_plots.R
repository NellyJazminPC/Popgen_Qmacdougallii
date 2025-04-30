
#Load libraries
library(ggplot2)
library(tidyr)
setwd("/Users/nelly/bioinfo/Popgen_Qmacdougallii/bin")
#Load the results from ADMIXTURE
tbl_cv <- read.table("../data/1.5.structure/admixture_output/chooseK.txt", 
                     header = F)
tbl_cv

# Cargar el archivo de texto
tbl_cv <- read.table("../data/1.5.structure/admixture_output/chooseK.txt", header = FALSE, sep = "\n", stringsAsFactors = FALSE)

# Extraer los valores de K y CV error usando expresiones regulares
tbl_cv_parsed <- data.frame(
  V1 = "K",
  V2 = as.numeric(gsub(".*\\(K=(\\d+)\\).*", "\\1", tbl_cv$V1)),
  V3 = as.numeric(gsub(".*:\\s*([0-9.]+)$", "\\1", tbl_cv$V1))
)

# Verificar el resultado
print(tbl_cv_parsed)

# Cross validation error plot
CV_error <- ggplot(data=tbl_cv, aes(x=K, y=V3)) +  # Usar la columna K en lugar de V2
  geom_line(size=1) + geom_point(size=3, alpha=1) +
  theme_bw() + ylab("CV error") + xlab("K") + scale_x_continuous(breaks=seq(0, 10, 1)) +
  theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), 
        axis.text.x  = element_text(vjust=0.5, size=14, face = "bold")) +
  theme(text = element_text(size=15)) + expand_limits(y=c(0.3,0.5)) + 
  scale_y_continuous(breaks = seq(0.25, 0.5, 0.05)) +
  theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),"cm"), text = element_text(size=20)) 

# Plot
CV_error

# Guardar el gráfico CV_error
ggsave(filename = "../results/CV_error_plot.png", plot = CV_error, width = 8, height = 6, dpi = 300)



# Next

# Load the Q files from ADMIXTURE

admix2.1=read.table(paste0("../data/1.5.structure/admixture_output/qmacd_ref_gen_rob.1.Q"))
admix2.2=read.table(paste0("../data/1.5.structure/admixture_output/qmacd_ref_gen_rob.2.Q"))
admix2.3=read.table(paste0("../data/1.5.structure/admixture_output/qmacd_ref_gen_rob.3.Q"))

# Load the metadata
qmacd=read.csv("../metadata/Qmacdougalli_79ind_.csv")
head(qmacd)
tail(qmacd)

# Rename columns
admix2.1<-`colnames<-`(admix2.1, c("K1"))
head(admix2.1)
admix2.2<-`colnames<-`(admix2.2, c("K2","K1"))
head(admix2.2)
admix2.3<-`colnames<-`(admix2.3, c("K2", "K3","K1"))
head(admix2.3)

# Using cbind to combine by rows and columns
admix2.1<- cbind(admix2.1,qmacd)
head(admix2.1)
admix2.2<- cbind(admix2.2,qmacd)
head(admix2.2)
admix2.3<- cbind(admix2.3,qmacd)
head(admix2.3)

# Using gather
admix2.1_gather<- gather (admix2.1, key= K, value=admixture, K1)
tail(admix2.1_gather)
admix2.2_gather<- gather (admix2.2, key= K, value=admixture, K1:K2)
tail(admix2.2_gather)
admix2.3_gather<- gather (admix2.3, key= K, value=admixture, K2:K1)
tail(admix2.3_gather)

### Plots

# K = 1

plot_admix2.1 <- ggplot(data=admix2.1_gather, aes(x=NUM_SAMPLE, y=admixture, fill=K)) + 
  geom_bar(stat="identity") + scale_fill_manual("K", values = c("K1" = "#E07E34")) + 
  ylab("")+ xlab("SITE")+ theme_bw() +
  theme(axis.title.x = element_text(size=16), 
        axis.title.y =element_text(size=16),
        axis.text.x  = element_text(size=10, angle = 90)) +
  theme(legend.title= element_blank(),text = element_text(size=20))+
  geom_vline(aes(xintercept=10.5)) + geom_vline(aes(xintercept=20.5)) + geom_vline(aes(xintercept=30.5)) + 
  geom_vline(aes(xintercept=40.5)) + geom_vline(aes(xintercept=50.5)) + geom_vline(aes(xintercept=54.5)) + 
  geom_vline(aes(xintercept=66.5)) + geom_vline(aes(xintercept=76.5))

plot_admix2.1 


## K = 2

plot_admix2.2 <- ggplot(data=admix2.2_gather, aes(x=NUM_SAMPLE, y=admixture, fill=K)) + 
  geom_bar(stat="identity") + scale_fill_manual("K", values = c("K2" = "#40B95B", "K1" = "#E07E34")) + 
  ylab("")+ xlab("SITE")+ theme_bw() +
  theme(axis.title.x = element_text(size=16), 
        axis.title.y =element_text(size=16),
        axis.text.x  = element_text(size=10, angle = 90)) +
  theme(legend.title= element_blank(),text = element_text(size=20))+
  geom_vline(aes(xintercept=10.5)) + geom_vline(aes(xintercept=20.5)) + geom_vline(aes(xintercept=30.5)) + 
  geom_vline(aes(xintercept=40.5)) + geom_vline(aes(xintercept=50.5)) + geom_vline(aes(xintercept=54.5)) + 
  geom_vline(aes(xintercept=66.5)) + geom_vline(aes(xintercept=76.5))

plot_admix2.2 

## K = 3

plot_admix2.3 <- ggplot(data=admix2.3_gather, aes(x=NUM_SAMPLE, y=admixture, fill=K)) + 
  geom_bar(stat="identity") + scale_fill_manual("K", values = c("K1" = "#00B1E8", "K3" = "#E07E34", "K2" = "#40B95B")) + 
  ylab("")+ xlab("SITE")+ theme_bw() +
  theme(axis.title.x = element_text(size=16), 
        axis.title.y =element_text(size=16),
        axis.text.x  = element_text(size=10, angle = 90)) +
  theme(legend.title= element_blank(),text = element_text(size=20))+
  geom_vline(aes(xintercept=10.5)) + geom_vline(aes(xintercept=20.5)) + geom_vline(aes(xintercept=30.5)) + 
  geom_vline(aes(xintercept=40.5)) + geom_vline(aes(xintercept=50.5)) + geom_vline(aes(xintercept=54.5)) + 
  geom_vline(aes(xintercept=66.5)) + geom_vline(aes(xintercept=76.5))

plot_admix2.3

# Exportar plot_admix2.1
ggsave(filename = "../results/plot_admix2.1.png", plot = plot_admix2.1, width = 8, height = 6, dpi = 300)

# Exportar plot_admix2.2
ggsave(filename = "../results/plot_admix2.2.png", plot = plot_admix2.2, width = 8, height = 6, dpi = 300)

# Exportar plot_admix2.3
ggsave(filename = "../results/plot_admix2.3.png", plot = plot_admix2.3, width = 8, height = 6, dpi = 300)

# Instalar patchwork si no está instalado
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}

# Cargar la librería patchwork
library(patchwork)

# Combinar los gráficos en una sola imagen
combined_plot <- plot_admix2.1 / plot_admix2.2 / plot_admix2.3 + 
  plot_layout(ncol = 1) + 
  plot_annotation(title = "Admixture Plots for K = 1, 2, and 3")

# Mostrar el gráfico combinado
print(combined_plot)

# Exportar el gráfico combinado
ggsave(filename = "../results/combined_admixture_plots.png", plot = combined_plot, width = 8, height = 18, dpi = 300)

# -------------------------------
# PLOT ALTITUDE
# -------------------------------

#"#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"
#VERDE     NARANJA   MORADO    ROSA       VERDE V   AMARILLO   MARRON   GRIS



plot_pop_alt_3<- ggplot(data=admix2.3_gather, aes(x=NUM_SAMPLE, y=ALT)) + 
  geom_line(stat="identity", color="#00cb5f") + geom_point(aes(colour=SITE_NAME, size=ALT, alpha=10/20)) +
  scale_colour_manual("SITE_NUM_NOM", values = c( "1CZ" = "#7570B3", "2MT" = "#075277", "3MC" = "#00B1E8","4MB" = "#1FC944","5CY" = "#E6AB02", "6LS" = "#E7298A","7PZ" = "#E07E34", "8CR" = "#F15858", "9IT" = "blue")) +
  ylab("Altitud msnm")+ xlab("Individuos")+ theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=12),
        axis.title.y =element_text(size=18),
        axis.text.x  = element_blank()) +
  geom_vline(aes(xintercept=50.5))
plot_pop_alt_3

# Exportar plot_pop_alt_3
ggsave(filename = "../results/plot_pop_alt_3.png", plot = plot_pop_alt_3, width = 8, height = 6, dpi = 300)



# -------------------------------
# FAST STRUCTURE PLOTS
# -------------------------------

#Load libraries
#library(ggplot2)
#library(tidyr)

# Load databases 
# Simple mode
fast_1.1_simple=read.table(paste0("../data/1.5.structure/faststructure_output/qmacd_ref_gen_rob.simple.1.meanQ"))
fast_1.2_simple=read.table(paste0("../data/1.5.structure/faststructure_output/qmacd_ref_gen_rob.simple.2.meanQ"))
fast_1.3_simple=read.table(paste0("../data/1.5.structure/faststructure_output/qmacd_ref_gen_rob.simple.3.meanQ"))

# Logistic mode
fast_1.1_log=read.table(paste0("../data/1.5.structure/faststructure_output/qmacd_ref_gen_rob.logistic.1.meanQ"))
fast_1.2_log=read.table(paste0("../data/1.5.structure/faststructure_output/qmacd_ref_gen_rob.logistic.2.meanQ"))
fast_1.3_log=read.table(paste0("../data/1.5.structure/faststructure_output/qmacd_ref_gen_rob.logistic.3.meanQ"))


# Load the metadata
#qmacd=read.csv("../metadata/metadata_qmacdo_79ind_gen.csv")
head(qmacd)
tail(qmacd)

# Rename columns
fast_1.1_simple<-`colnames<-`(fast_1.1_simple, c("K1"))
head(fast_1.1_simple)
fast_1.2_simple<-`colnames<-`(fast_1.2_simple, c("K2","K1"))
head(fast_1.2_simple)
fast_1.3_simple<-`colnames<-`(fast_1.3_simple, c("K3", "K2","K1"))
head(fast_1.3_simple)

fast_1.1_log<-`colnames<-`(fast_1.1_log, c("K1"))
head(fast_1.1_log)
fast_1.2_log<-`colnames<-`(fast_1.2_log, c("K2","K1"))
head(fast_1.2_log)
fast_1.3_log<-`colnames<-`(fast_1.3_log, c("K3", "K2","K1"))
head(fast_1.3_log)


# Using cbind to combine by rows and columns

# Simple mode
fast_1.1_simple<- cbind(fast_1.1_simple,qmacd)
head(fast_1.1_simple)
fast_1.2_simple<- cbind(fast_1.2_simple,qmacd)
head(fast_1.2_simple)
fast_1.3_simple<- cbind(fast_1.3_simple,qmacd)
head(fast_1.3_simple)

# Logistic mode

fast_1.1_log<- cbind(fast_1.1_log,qmacd)
head(fast_1.1_log)
fast_1.2_log<- cbind(fast_1.2_log,qmacd)
head(fast_1.2_log)
fast_1.3_log<- cbind(fast_1.3_log,qmacd)
head(fast_1.3_log)


# Using gather

# Simple mode

fast_1.1_simple_gather<- gather (fast_1.1_simple, key= K, value=admixture, K1)
tail(fast_1.1_simple_gather)
fast_1.2_simple_gather<- gather (fast_1.2_simple, key= K, value=admixture, K2:K1)
tail(fast_1.2_simple_gather)
fast_1.3_simple_gather<- gather (fast_1.3_simple, key= K, value=admixture, K3:K1)
tail(fast_1.3_simple_gather)

# Logistic mode

fast_1.1_log_gather<- gather (fast_1.1_log, key= K, value=admixture, K1)
tail(fast_1.1_log_gather)
fast_1.2_log_gather<- gather (fast_1.2_log, key= K, value=admixture, K2:K1)
tail(fast_1.2_log_gather)
fast_1.3_log_gather<- gather (fast_1.3_log, key= K, value=admixture, K3:K1)
tail(fast_1.3_log_gather)

### Plots

# K = 1

# We can change from "_log_" to "_simple_" to obtain the plot from the simple mode with K=1

plot_fast_1.1 <- ggplot(data=fast_1.1_log_gather, aes(x=NUM_SAMPLE, y=admixture, fill=K)) + 
  geom_bar(stat="identity") + scale_fill_manual("K", values = c("K1" ="#E07E34" )) + 
  ylab(" ")+ xlab("SITES")+ theme_bw() +
  theme(axis.title.x = element_text(size=18), 
        axis.title.y =element_blank(),
        axis.text.x  = element_text(size=10, angle = 90)) +
  theme(legend.title= element_blank(),text = element_text(size=20))+
  geom_vline(aes(xintercept=10.5)) + geom_vline(aes(xintercept=20.5)) + geom_vline(aes(xintercept=30.5)) + 
  geom_vline(aes(xintercept=40.5)) + geom_vline(aes(xintercept=50.5)) + geom_vline(aes(xintercept=54.5)) + 
  geom_vline(aes(xintercept=66.5)) + geom_vline(aes(xintercept=76.5)) + theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),"cm"))

plot_fast_1.1 

# K = 2
# We can change from "_log_" to "_simple_" to obtain the plot from the simple mode with K=2

plot_fast_1.2 <- ggplot(data=fast_1.2_log_gather, aes(x=NUM_SAMPLE, y=admixture, fill=K)) + 
  geom_bar(stat="identity") + scale_fill_manual("K", values = c("K1" ="#E07E34" , "K2" = "#40B95B")) + 
  ylab(" ")+ xlab("SITES")+ theme_bw() +
  theme(axis.title.x = element_text(size=18), 
        axis.title.y =element_blank(),
        axis.text.x  = element_text(size=10, angle = 90)) +
  theme(legend.title= element_blank(),text = element_text(size=20))+
  geom_vline(aes(xintercept=10.5)) + geom_vline(aes(xintercept=20.5)) + geom_vline(aes(xintercept=30.5)) + 
  geom_vline(aes(xintercept=40.5)) + geom_vline(aes(xintercept=50.5)) + geom_vline(aes(xintercept=54.5)) + 
  geom_vline(aes(xintercept=66.5)) + geom_vline(aes(xintercept=76.5)) + theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),"cm"))

plot_fast_1.2


# K = 3
# We can change from "_log_" to "_simple_" to obtain the plot from the simple mode with K=3

plot_fast_1.3 <- ggplot(data=fast_1.3_log_gather, aes(x=NUM_SAMPLE, y=admixture, fill=K)) + 
  geom_bar(stat="identity") + scale_fill_manual("K", values = c("K1" = "#00B1E8", "K3" = "#E07E34", "K2" = "#40B95B")) + 
  ylab(" ")+ xlab("SITES")+ theme_bw() +
  theme(axis.title.x = element_text(size=18), 
        axis.title.y =element_blank(),
        axis.text.x  = element_text(size=10, angle = 90)) +
  theme(legend.title= element_blank(),text = element_text(size=20))+
  geom_vline(aes(xintercept=10.5)) + geom_vline(aes(xintercept=20.5)) + geom_vline(aes(xintercept=30.5)) + 
  geom_vline(aes(xintercept=40.5)) + geom_vline(aes(xintercept=50.5)) + geom_vline(aes(xintercept=54.5)) + 
  geom_vline(aes(xintercept=66.5)) + geom_vline(aes(xintercept=76.5)) + theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),"cm"))

plot_fast_1.3

# Exportar plot_fast_1.1
ggsave(filename = "../results/plot_fast_1.1.png", plot = plot_fast_1.1, width = 8, height = 6, dpi = 300)

# Exportar plot_fast_1.2
ggsave(filename = "../results/plot_fast_1.2.png", plot = plot_fast_1.2, width = 8, height = 6, dpi = 300)

# Exportar plot_fast_1.3
ggsave(filename = "../results/plot_fast_1.3.png", plot = plot_fast_1.3, width = 8, height = 6, dpi = 300)

# Exportar los tres gráficos en una sola imagen
# Cargar la librería patchwork
library(patchwork)

# Combinar los gráficos en una sola imagen
combined_fast_plot <- plot_fast_1.1 / plot_fast_1.2 / plot_fast_1.3 + 
  plot_layout(ncol = 1) + 
  plot_annotation(title = "Fast Structure Plots for K = 1, 2, and 3")

# Mostrar el gráfico combinado
print(combined_fast_plot)

# Exportar el gráfico combinado
ggsave(filename = "../results/combined_fast_structure_plots.png", plot = combined_fast_plot, width = 8, height = 18, dpi = 300)