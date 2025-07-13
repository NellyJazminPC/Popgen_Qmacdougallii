
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
CV_error <- ggplot(data=tbl_cv_parsed, aes(x=V2, y=V3)) +  # Usar la columna K en lugar de V2
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
  geom_bar(stat="identity") + scale_fill_manual("K", values = c("K1" = "#0072B2")) + 
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
  geom_bar(stat="identity") + scale_fill_manual("K", values = c("K2" = "#E69F00", "K1" = "#0072B2")) + 
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
  geom_bar(stat="identity") + scale_fill_manual("K", values = c("K1" = "#0072B2", "K3" = "#E69F00", "K2" = "#40B95B")) + 
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
  geom_bar(stat="identity") + scale_fill_manual("K", values = c("K1" ="#0072B2" )) + 
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
  geom_bar(stat="identity") + scale_fill_manual("K", values = c("K1" ="#0072B2" , "K2" = "#E69F00")) + 
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
  geom_bar(stat="identity") + scale_fill_manual("K", values = c("K1" = "#0072B2", "K3" = "#E69F00", "K2" = "#40B95B")) + 
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




# -------------------------------
# Plot fastStructure K=2 with scatterpie on a map
# -------------------------------

library(dplyr)
library(tidyr)
library(scatterpie)
library(ggplot2)

# Crear los gráficos de pastel para cada sitio.


# Calcular el promedio por sitio de K1 y K2
site_props_avg <- fast_1.2_log_gather %>%
  group_by(SITE_NAME) %>%
  summarise(
    Q1 = mean(admixture[K == "K1"], na.rm = TRUE),
    Q2 = mean(admixture[K == "K2"], na.rm = TRUE)
  )

# Convertir a formato largo para ggplot2
site_props_long <- site_props_avg %>%
  pivot_longer(cols = c(Q1, Q2), names_to = "Cluster", values_to = "Prop")

# Guardar cada gráfico de pastel como PNG transparente
unique_sites <- unique(site_props_long$SITE_NAME)

for (site in unique_sites) {
  pie_data <- filter(site_props_long, SITE_NAME == site)
  pie_plot <- ggplot(pie_data, aes(x = "", y = Prop, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c("Q1" = "#0072B2", "Q2" = "#E69F00")) +
    theme_void() +
    theme(legend.position = "none", plot.background = element_rect(fill = "transparent", color = NA)) +
    labs(title = site)
  
  ggsave(
    filename = paste0("../results/piechart_", site, ".png"),
    plot = pie_plot,
    width = 3, height = 3, dpi = 300, bg = "transparent"
  )
}

library(ggplot2)
library(ggspatial)

ggplot(qmacd, aes(x = long, y = lat)) +
  annotation_map_tile(type = "esri", zoom = 5) +  # Fondo satelital ESRI
  geom_point(aes(color = SITE_NAME), size = 3, alpha = 0.8) +
  geom_text(aes(label = SITE_NAME), vjust = -1, size = 3, color = "white") +
  scale_color_manual(values = c(
    "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#E69F00", "#000000", "#999999"
  )) +
  theme_minimal() +
  labs(title = "Localización de individuos", x = "Longitud", y = "Latitud", color = "Sitio")



  library(ggplot2)
  library(ggspatial)
  
  ggplot(qmacd, aes(x = long, y = lat)) +
    annotation_map_tile(type = "stamen", zoom = 10, source = "terrain") +  # Fondo de relieve
    geom_point(aes(color = SITE_NAME), size = 3, alpha = 0.8) +
    geom_text(aes(label = SITE_NAME), vjust = -1, size = 3, color = "black") +
    scale_color_manual(values = c(
      "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#E69F00", "#000000", "#999999"
    )) +
    theme_minimal() +
    labs(title = "Localización de individuos", x = "Longitud", y = "Latitud", color = "Sitio")



    expand <- 0.05
    xlim <- c(min(qmacd$long) - expand, max(qmacd$long) + expand)
    ylim <- c(min(qmacd$lat) - expand, max(qmacd$lat) + expand)
    
    ggplot(qmacd, aes(x = long, y = lat)) +
      annotation_map_tile(type = "esri", zoom = 7) +
      geom_point(aes(color = SITE_NAME), size = 3, alpha = 0.8) +
      geom_text(aes(label = SITE_NAME), vjust = -1, size = 3, color = "white") +
      scale_color_manual(values = c(
        "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#E69F00", "#000000", "#999999"
      )) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme_minimal() +
      labs(title = "Localización de individuos", x = "Longitud", y = "Latitud", color = "Sitio")


library(ggplot2)
library(ggspatial)

# Ejemplo de coordenadas en CDMX
cdmx_df <- data.frame(
  SITE_NAME = c("Sitio1", "Sitio2", "Sitio3"),
  long = c(-99.1332, -99.1450, -99.1200),
  lat = c(19.4326, 19.4400, 19.4200)
)

# Fondo satelital ESRI
ggplot(cdmx_df, aes(x = long, y = lat)) +
  annotation_map_tile(type = "esri", zoom = 12) +
  geom_point(aes(color = SITE_NAME), size = 4) +
  geom_text(aes(label = SITE_NAME), vjust = -1, size = 4, color = "white") +
  scale_color_manual(values = c("#E69F00", "#0072B2", "#56B4E9")) +
  theme_minimal() +
  labs(title = "Ejemplo: Sitios en Ciudad de México", x = "Longitud", y = "Latitud", color = "Sitio")


    library(ggplot2)
  library(ggspatial)
  
  cdmx_df <- data.frame(
    SITE_NAME = c("Sitio1", "Sitio2", "Sitio3"),
    long = c(-99.1332, -99.1450, -99.1200),
    lat = c(19.4326, 19.4400, 19.4200)
  )
  
  ggplot(cdmx_df, aes(x = long, y = lat)) +
    annotation_map_tile(type = "esri", zoom = 12) +
    geom_point(aes(color = SITE_NAME), size = 4) +
    geom_text(aes(label = SITE_NAME), vjust = -1, size = 4, color = "white") +
    scale_color_manual(values = c("#E69F00", "#0072B2", "#56B4E9")) +
    coord_sf(crs = 4326, datum = NA) +
    theme_minimal() +
    labs(title = "Ejemplo: Sitios en Ciudad de México", x = "Longitud", y = "Latitud", color = "Sitio")


# -------------------------------
library(ggplot2)
library(ggspatial)

cdmx_df <- data.frame(
  SITE_NAME = c("Sitio1", "Sitio2", "Sitio3"),
  long = c(-99.1332, -99.1450, -99.1200),
  lat = c(19.4326, 19.4400, 19.4200)
)

ggplot(cdmx_df, aes(x = long, y = lat)) +
  annotation_map_tile(type = "esri", zoom = 12) +
  geom_point(aes(color = SITE_NAME), size = 4) +
  geom_text(aes(label = SITE_NAME), vjust = -1, size = 4, color = "white") +
  scale_color_manual(values = c("#E69F00", "#0072B2", "#56B4E9")) +
  coord_sf(crs = 4326, datum = NA) +
  theme_minimal() +
  labs(title = "Ejemplo: Sitios en Ciudad de México", x = "Longitud", y = "Latitud", color = "Sitio")

# -------------------------------
ggplot(cdmx_df, aes(x = long, y = lat)) +
  annotation_map_tile(type = "stamen", zoom = 12, source = "terrain") +
  geom_point(aes(color = SITE_NAME), size = 4) +
  geom_text(aes(label = SITE_NAME), vjust = -1, size = 4, color = "black") +
  scale_color_manual(values = c("#E69F00", "#0072B2", "#56B4E9")) +
  coord_sf(crs = 4326, datum = NA) +
  theme_minimal() +
  labs(title = "Ejemplo: Sitios en Ciudad de México", x = "Longitud", y = "Latitud", color = "Sitio")


# -------------------------------


# Librerías necesarias
libs <- c("sf", "terra", "elevatr", "rayshader", "httr", "jsonlite", "png")
invisible(lapply(libs, require, character.only = TRUE))

# 1. Bounding box como polígono sf
bbox <- st_as_sfc(st_bbox(c(xmin = -96.6, xmax = -96.2,
                            ymin = 17.4, ymax = 17.8),
                          crs = 4326))
bbox_sf <- st_sf(geometry = bbox)

# 2. Descarga el DEM
dem <- get_elev_raster(bbox_sf, z = 10, clip = "locations")

# 3. Convierte el raster a matriz para rayshader
elmat <- raster_to_matrix(dem)

# Reemplaza NA o valores fuera de rango por 0
elmat[is.na(elmat)] <- 0
elmat[elmat < 0] <- 0

# 4. Sombrea y renderiza en 3D
hill <- sphere_shade(elmat, texture = "desert")
plot_3d(elmat, hill, zscale = 12, windowsize = c(1200, 800))
render_snapshot("sierra_juarez_3d.png", title_text = "Study area")
# Si quieres cerrar la ventana 3D después:
# rgl::rgl.close()


str(hill)
range(hill)