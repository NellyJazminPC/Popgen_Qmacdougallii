# Compute minimum year and minimum Ne_median for each demography summary file
import pandas as pd

files = [
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_4.2e8_100y_demography_analysis.final.summary',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_1.01e8_100y_demography_analysis.final.summary',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_5.2e8_100y_demography_analysis.final.summary',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_4.2e8_50y_demography_analysis.final.summary',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_5.2e8_50y_demography_analysis.final.summary',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_1.01e8_50y_demography_analysis.final.summary'
]

results = []
for f in files:
    df = pd.read_csv(f, sep='\t', comment='#')
    min_year = df['year'].min()
    min_ne = df['Ne_median'].min()
    label = f.split('/')[-1].replace('pop1_','').replace('_demography_analysis.final.summary','')
    results.append({'file': label, 'min_year': min_year, 'min_Ne_median': min_ne})

summary_df = pd.DataFrame(results)
print(summary_df)

# Exportar el DataFrame a CSV en el directorio especificado
summary_df.to_csv('/home/n311pc/bioinfo/Popgen_Qmacdougallii/results/demography/min_Ne_summary.csv', index=False)
print('Resumen exportado como min_Ne_summary.csv en results/demography')

#######
# Plot six lines with additional paleo events (Bølling–Allerød, Younger Dryas, Clovis, Megafauna extinction)

import pandas as pd, matplotlib.pyplot as plt

files = [
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_1.01e8_50y_demography_analysis.final.summary',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_1.01e8_100y_demography_analysis.final.summary',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_4.2e8_50y_demography_analysis.final.summary',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_4.2e8_100y_demography_analysis.final.summary',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_5.2e8_50y_demography_analysis.final.summary',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_5.2e8_100y_demography_analysis.final.summary'
]

# color mapping consistent with previous
color_map = {
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_1.01e8_50y_demography_analysis.final.summary': 'firebrick',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_1.01e8_100y_demography_analysis.final.summary': 'darkred',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_4.2e8_50y_demography_analysis.final.summary': 'gold',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_4.2e8_100y_demography_analysis.final.summary': 'goldenrod',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_5.2e8_50y_demography_analysis.final.summary': 'royalblue',
    '/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/pop1_5.2e8_100y_demography_analysis.final.summary': 'navy'
}

plt.figure(figsize=(13,8))
for f in files:
    df = pd.read_csv(f, sep='\t', comment='#')
    x = df['year']/1000
    ne = df['Ne_median']/1000
    ci_low = df['Ne_2.5%']/1000
    ci_up = df['Ne_97.5%']/1000
    # Extraer solo el nombre del archivo y construir la etiqueta corta
    label_core = f.split('/')[-1].replace('pop1_','').replace('_demography_analysis.final.summary','')
    parts = label_core.split('_')
    mu, gen = parts[0], parts[1]
    mu_fmt = mu.replace('1.01e8', r'1.01$\times 10^{-8}$') \
               .replace('4.2e8', r'4.2$\times 10^{-8}$') \
               .replace('5.2e8', r'5.2$\times 10^{-8}$')    
    label = f"{mu_fmt}_{gen}"
    c = color_map[f]
    plt.fill_between(x, ci_low, ci_up, color=c, alpha=0.15)
    plt.plot(x, ne, lw=2, color=c, alpha=0.8, label=label)

# Existing markers: LGM and Holocene start
plt.axvspan(19,26, color='grey', alpha=0.2, label='LGM')
plt.axvline(11.7, color='green', ls='--', lw=2, alpha=0.7, label='Holocene start (11.7 ka)')


plt.xscale('log')
plt.yscale('log')
plt.xlim(0.05, 20000)
plt.ylim(0.1, 300000)
plt.xlabel('ka Before Present', fontsize=15)
plt.ylabel('$N_e$ ($\ x10^{3}$)', fontsize=15)

# Ejes con números completos
plt.xticks([0.1, 1, 10, 100, 1000, 10000], ['0.1', '1', '10', '100', '1000', '10000'])
plt.yticks([0.1, 1, 10, 100, 1000, 10000, 100000], ['0.1', '1', '10', '100', '1000', '10000', '100000'])

plt.legend(fontsize=12, loc='upper left', ncol=2, framealpha=0.95)
plt.tight_layout()

# Exportar el gráfico al directorio de resultados
plt.savefig('/home/n311pc/bioinfo/Popgen_Qmacdougallii/results/demography/demography_paleo_events.png', dpi=300)
plt.show()
