# Ensamble e identificación de las variantes (SNPs) con ipyRAD

---

![alt text](Diagrama_flujo_ensambles_ipyRAD.png)

- Se utilizaron las modalidades _de novo_ y con _genoma de referencia_. 

- Se utilizaron los genomas de _Quercus lobata_ y _Q. robur_.

- Para los ensambles en ipyRAD se hicieron dos ramas (branches).

Los parámetros que se modificaron son:

```python
# 1
("clust_threshold", "0.70")
('pop_assign_file', 'popmap_9sites_2zones_trim01.txt')

# 2
("min_samples_locus", 79)
("clust_threshold", "0.90")
```

- El `clust_threshold` es el nivel de similitud de secuencia; en el que dos secuencias se identifican como homólogas y, por lo tanto, se agrupan. El valor se ingresa como un decimal (por ejemplo, 0.90). Por otro lado, el programa no recomienda utilizar valores superiores a 0,95 ya que es posible que las secuencias homólogas no se agrupen en un umbral tan alto debido a la presencia de N, indeles, errores de secuenciación o polimorfismos.
- El `min_samples_locus` es el número mínimo de muestras que deben tener datos en un locus determinado para que se conserven en el conjunto de datos final. El programa indica que si se ingresa un número igual al número total de muestras en su conjunto de datos, solo devolverá loci que tengan datos compartidos en todas las muestras, mientras que si ingresa un valor más bajo, devolverá una matriz más dispersa, incluidos los loci para los que al menos cuatro muestras contengan datos por ejemplo. Este parámetro se anula si se ingresan valores `min_samples en el archivo de popmap, como fue en el primer caso.

---

- A los archivos VCF generados por los ensambles se les aplicó un filtro con TASSEL v5.2.93:

> Filter > Filter Genotype Table Sites 

- Site Min Count: 79
- Site Min Allele Freq: 0,00
- Site Max Allele Freq: 1,00
- Remove Minor SNP States
- Remove Sites With Indels

---

El número de SNPs que se identificaron para los 18 archivos VCF varió entre las modalidades de novo y con genoma de referencia, por lo que seguimos explorando los VCF obtenidos con los genomas de referencia.

[**Archivo .CSV**](SNPs_18_VCF_ipyrad.csv)

| VCF_SNPs | pipeline | mode          | input | ipyrad_branches    | wo_filter | 1rst_filter | 2nd_filter |
|----------|----------|---------------|-------|--------------------|-----------|-------------|------------|
| 1        | ipyRAD   | denovo        | trim01| denovo_trim01_1    | 4087      | 3471        | 3429       |
| 2        | ipyRAD   | denovo        | trim02| denovo_trim02_1    | 986       | 638         | 628        |
| 3        | ipyRAD   | denovo        | trim03| denovo_trim03_1    | 3582      | 3002        | 2965       |
| 4        | ipyRAD   | denovo        | trim01| denovo_trim01_2    | 4087      | 3471        | 3429       |
| 5        | ipyRAD   | denovo        | trim02| denovo_trim02_2    | 986       | 638         | 628        |
| 6        | ipyRAD   | denovo        | trim03| denovo_trim03_2    | 3582      | 3002        | 2965       |
| 7        | ipyRAD   | ref_gen_qlobata| trim01| ref_gen_qlob_trim01_1| 7559   | 5488        | 5415       |
| 8        | ipyRAD   | ref_gen_qlobata| trim02| ref_gen_qlob_trim02_1| 1946   | 1094        | 1076       |
| 9        | ipyRAD   | ref_gen_qlobata| trim03| ref_gen_qlob_trim03_1| 5487   | 4435        | 4374       |
| 10       | ipyRAD   | ref_gen_qlobata| trim01| ref_gen_qlob_trim01_2| 7559   | 5488        | 5415       |
| 11       | ipyRAD   | ref_gen_qlobata| trim02| ref_gen_qlob_trim02_2| 1946   | 1094        | 1076       |
| 12       | ipyRAD   | ref_gen_qlobata| trim03| ref_gen_qlob_trim03_2| 5487   | 4435        | 4374       |
| 13       | ipyRAD   | ref_gen_qrobur | trim01| ref_gen_qrob_trim01_1| 7611   | 5501        | 5426       |
| 14       | ipyRAD   | ref_gen_qrobur | trim02| ref_gen_qrob_trim02_1| 1815   | 1048        | 1027       |
| 15       | ipyRAD   | ref_gen_qrobur | trim03| ref_gen_qrob_trim03_1| 5382   | 4291        | 4230       |
| 16       | ipyRAD   | ref_gen_qrobur | trim01| ref_gen_qrob_trim01_2| 7611   | 5501        | 5426       |
| 17       | ipyRAD   | ref_gen_qrobur | trim02| ref_gen_qrob_trim02_2| 1815   | 1048        | 1027       |
| 18       | ipyRAD   | ref_gen_qrobur | trim03| ref_gen_qrob_trim03_2| 5382   | 4291        | 4230       |
