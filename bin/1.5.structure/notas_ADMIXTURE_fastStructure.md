# Ejecución del programa ADMIXTURE 1.3.0

## Requisitos previos

1. Instalación de [ADMIXTURE 1.3.0.](https://dalexander.github.io/admixture/download.html) 
2. Conversión de los datos en formato PLINK (.bed, .bim, .fam) con TASSEL


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
