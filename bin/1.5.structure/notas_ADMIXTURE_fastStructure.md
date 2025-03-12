# Ejecución del programa ADMIXTURE 1.3.0

## Requisitos previos

1. Instalación de [ADMIXTURE 1.3.0.](https://dalexander.github.io/admixture/download.html) 
2. Conversión de los datos en formato PLINK (.bed, .bim, .fam) 

    2.1 Se convirtip a `.ped`y `.map` con TASSEL

    2.2 Con el programa PLINK se convirtio a `.bed`, `.bim` y `.fam`


    ## Conversión de archivos `.ped` y `.map` a `.bed`, `.bim` y `.fam`

    Para convertir archivos `.ped` y `.map` a `.bed`, `.bim` y `.fam` utilizando PLINK, puedes usar el siguiente comando:

    ```bash
    plink --file <input_file> --make-bed --out <output_file>
    ```

    - `<input_file>`: Nombre del archivo de entrada sin la extensión.
    - `<output_file>`: Nombre del archivo de salida sin la extensión.

    ### Ejemplo de uso

    ```bash
    plink --file mydata --make-bed --out mydata_converted
    ```

    Este comando convertirá los archivos `mydata.ped` y `mydata.map` en `mydata_converted.bed`, `mydata_converted.bim` y `mydata_converted.fam`.

## Comando básico

```bash
admixture <input_file>.bed <K>
```

- `<input_file>`: Nombre del archivo de entrada sin la extensión.
- `<K>`: Número de grupos ancestrales.

## Ejemplo de uso

```bash
admixture mydata.bed 3
```

Este comando ejecutará ADMIXTURE en el archivo `mydata.bed` asumiendo 3 grupos ancestrales.

## Opciones adicionales

- `-B`: Realiza bootstrapping.
- `--cv`: Realiza validación cruzada.

### Ejemplo con validación cruzada

```bash
admixture --cv mydata.bed 3
```

Este comando ejecutará ADMIXTURE con validación cruzada para 3 grupos ancestrales.

## Referencias

Para más información, consulta la [documentación oficial de ADMIXTURE](https://dalexander.github.io/admixture/).

