# Cargar las librerías
library(adegenet)  # Manejo de objetos genéticos y estadísticas de diversidad
library(vcfR)      # Lectura y conversión de archivos VCF
library(dartR)     # Conversión de genlight a genind y genclone

# Load the VCF file into R
qmacd_vcf <- read.vcfR("../data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf")

# Convert the VCF to a genlight object
qmacd_genlight <- vcfR2genlight(qmacd_vcf)

# Load the metadata file
pop.metadata <- read.csv("../metadata/Qmacdougalli_79ind_.csv")
head(pop.metadata)

# Verify that the individual names in the VCF match the metadata
all(colnames(qmacd_vcf@gt)[-1] == pop.metadata$ID)

# Set ploidy
ploidy(qmacd_genlight) <- 2
qmacd_genlight@ploidy

# Add population levels using SITE and ZONE metadata
pop(qmacd_genlight) <- pop.metadata$SITE  # Assign SITE as population
qmacd_genlight_pop <- qmacd_genlight  # Create a copy to use POP
pop(qmacd_genlight_pop) <- pop.metadata$POP  # Assign POP as population

# Convert from genlight to genind and genclone (DartR)
qmacd_genind <- gl2gi(qmacd_genlight, v=1)
qmacd_genclone <- as.genclone(qmacd_genind)

# For POP
qmacd_genind_pop <- gl2gi(qmacd_genlight_pop, v=1)
qmacd_genclone_pop <- as.genclone(qmacd_genind_pop)


# Calcular estadísticas de diversidad genética

# ajustar

snps_divers <- summary(qmacd_genind)
names(snps_divers)

plot(snps_divers$n.by.pop, snps_divers$pop.n.all)

barplot(snps_divers$loc.n.all)

barplot(snps_divers$Hexp-snps_divers$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")
barplot(snps_divers$pop.eff, main="Sample sizes per population",
        ylab="Number of genotypes",las=3)

bartlett.test(list(snps_divers$Hexp,snps_divers$Hobs))

t.test(snps_divers$Hexp,snps_divers$Hobs,pair=T,var.equal=TRUE,alter="greater")



