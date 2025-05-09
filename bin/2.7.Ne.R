# Ne
install.packages("strataG")  # Si no está instalado

install.packages("remotes")
remotes::install_github("ericarcher/strataG")


# make sure you have devtools installed
if (!require('devtools')) install.packages('devtools')
# install strataG latest version
devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)

library(strataG)

library(adegenet)

# Convertir un objeto genind a gtypes
gtypes_data <- genind2gtypes(qmacd_genind, strata = pop(qmacd_genind))


# Calcular Ne basado en el desequilibrio de ligamiento
ne_results <- linkNe(gtypes_data, ci = 0.95, nrep = 1000)

# Ver los resultados
print(ne_results)


library(PopGenome)

# Cargar datos desde un archivo VCF
genome_data <- readVCF("../data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf", numcols = 1000)

# Calcular diversidad nucleotídica (pi) como proxy para Ne
diversity <- genome_data@nucleotide.diversity


library(hierfstat)

# Calcular F-statistics
f_stats <- basic.stats(qmacd_genind)

# Usar Fis para inferir Ne (requiere fórmulas adicionales)
fis_values <- f_stats$Fis

