# 🧬 RNA Bootstrap Analysis for miRNA Expression
📌 Overview

This repository contains an R script for performing bootstrap analysis on miRNA expression data to assess statistical significance using permutation tests. This approach helps mitigate false positives and ensures robust results when analyzing small sample sizes.

🔬 Methodology

Simulated Data Generation:

Creates randomly generated biological data using normal distributions.

Three groups simulated: ED, Ward, and ICU.

Ensures reproducibility by setting a random seed.

Bootstrap Analysis:

Performs 1000 permutations by randomly selecting samples.

Conducts t-tests for each miRNA across groups.

Calculates p-values for statistical assessment.

Ranking and Visualization:

Computes mean and median p-values.

Determines frequency of statistically significant comparisons.

Boxplots illustrate the distribution of p-values.

📦 Dependencies

Ensure you have the following R libraries installed:

install.packages("ggplot2")

🛠️ Usage

1️⃣ Clone this repository:

git clone https://github.com/claudiacastrillon/RNA_Bootstrap_Analysis.git

2️⃣ Navigate to the project folder:

cd RNA_Bootstrap_Analysis

3️⃣ Run the R script:

source("bootstrapping_boxplots.r")

📊 Output

Boxplots of p-values for each miRNA across 1000 permutations.

Ranked list of miRNAs based on statistical significance.

PNG file (boxplot_p_values_14micros.png) generated with the visualization.

🤝 Contributions

Feel free to contribute by submitting pull requests or reporting issues!

📜 License

This project is open-source. See LICENSE for details.

📩 Contact

For inquiries, contact claudiacastrillon via GitHub. 💡


