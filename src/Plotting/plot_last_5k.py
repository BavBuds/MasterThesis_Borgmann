import matplotlib.pyplot as plt
import pandas as pd

# Input and output paths
plot_txt = "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/PSMC/Plotting/Bolinus_brandaris.plot.txt"
output_plot = "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/PSMC/Plotting/Bolinus_brandaris_last_5k.pdf"

# Read the file, keeping only rows starting with 'h'
data = pd.read_csv(
    plot_txt,
    sep="\t",
    header=None,
    names=["Raw"],
    comment="T",  # Skip lines starting with "T", "R", or "N"
    skiprows=lambda x: not str(x).startswith("h"),  # Skip non-'h' rows
    engine="python",
)

# Split 'h' rows into Time and Ne
data[["Time", "Ne"]] = data["Raw"].str.split(expand=True)[[1, 2]].astype(float)

# Crop to the last 5,000 years
data_cropped = data[data["Time"] <= 5000]

# Create the plot
plt.figure(figsize=(8, 6))
plt.plot(data_cropped["Time"], data_cropped["Ne"], label="Effective Population Size", color='blue')

# Customize ticks and labels
plt.xticks(range(0, 5001, 500))
plt.xlabel("Time (years before present)")
plt.ylabel("Effective Population Size (Ne)")
plt.title("PSMC - Last 5,000 Years")
plt.grid(True)
plt.legend()

# Save the plot
plt.savefig(output_plot, format="pdf")
print(f"Plot saved to {output_plot}")
