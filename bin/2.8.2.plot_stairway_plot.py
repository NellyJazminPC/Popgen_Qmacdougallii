# Load the stairway plot results file and inspect the first few rows
import pandas as pd

# Read the file, skipping any comment lines if present
with open('/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/one_pop/pop1_demography_analysis.final.summary', 'r') as f:
    lines = f.readlines()

# Find the first non-comment line (header)
header_idx = 0
for i, line in enumerate(lines):
    if not line.startswith('#') and line.strip() != '':
        header_idx = i
        break

# Load the dataframe from the correct header
summary_df = pd.read_csv('/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.8.demography/stairway/one_pop/pop1_demography_analysis.final.summary', sep='\t', header=header_idx)
# Show the head of the dataframe to confirm correct loading
print(summary_df.head())


# Prepare and plot the stairway plot results in the style of the reference image
import matplotlib.pyplot as plt
import numpy as np

# Extract relevant columns
x = summary_df['year'] / 1000  # Convert years to ka (thousands of years)
y = summary_df['Ne_median'] / 1000  # Convert Ne to thousands
ci_lower = summary_df['Ne_2.5%'] / 1000
ci_upper = summary_df['Ne_97.5%'] / 1000

# Set up the plot
plt.figure(figsize=(7, 6))

# Plot confidence interval as a filled area
plt.fill_between(x, ci_lower, ci_upper, color='mediumpurple', alpha=0.3, label='95% CI')

# Plot the median Ne line
plt.plot(x, y, color='purple', lw=2, label='Median $N_e$')

# Add vertical reference line at 20 ka BP
#plt.axvline(20, color='blue', linestyle='--', lw=2)

# Add annotation for SNPs (using the value from your image)
plt.text(0.08, 0.4, '5426 SNPs', fontsize=12, color='black', transform=plt.gca().transAxes)

# Add panel label and title
#plt.text(-0.15, 1.05, '(c)', fontsize=16, fontweight='bold', transform=plt.gca().transAxes)
plt.title('One single pop', fontsize=16, pad=20)

# Franja del LGM (26-19 ka BP)
plt.axvspan(19, 26, color='gray', alpha=0.3, label='LGM (26-19 ka BP)')
plt.text(22.5, plt.ylim()[1]/2, 'LGM', color='black', fontsize=13, ha='center', va='center', rotation=90, alpha=0.7)


# Set log-log axes
plt.xscale('log')
plt.yscale('log')

# Set axis labels
plt.xlabel('Years before present', fontsize=14)
plt.ylabel('$N_e$ ($\	x 10^3$)', fontsize=14)

# Set axis limits and ticks to match the reference
plt.xlim(0.05, 1000000)
plt.ylim(0.25, 20000)
plt.xticks([0.5, 2, 10, 100, 500, 1000000], ['0.5k', '2k', '10k','100k','500k','1M'])
plt.yticks([0.25, 1, 4, 16, 30, 40, 64, 256, 4000], ['0.25', '1', '4', '16', '30', '40', '64', '256', '4000'])

# Remove top and right spines for a cleaner look
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

# Add grid for better readability
plt.tight_layout()
# Guardar la figura en el directorio especificado
plt.savefig('/home/n311pc/bioinfo/Popgen_Qmacdougallii/results/demography/one_single_pop.png', dpi=300)
plt.show()

print('Stairway plot styled, displayed y exportado como one_single_pop.png en results/demography')
# ...existing code...
print('Stairway plot styled and displayed as in the reference image from paper Ortego et al., 2023')