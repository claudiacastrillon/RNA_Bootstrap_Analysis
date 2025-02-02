# Install and load necessary packages
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

# List of microRNAs
micrornas <- c(
  "hsa-miR-142-5p", "hsa-miR-148a-3p", "hsa-miR-4508", "hsa-miR-199a-3p",
  "hsa-miR-197-3p", "hsa-miR-193a-5p", "hsa-miR-223-5p", "hsa-miR-1307-3p",
  "hsa-miR-1246", "hsa-let-7f-5p", "hsa-miR-221-3p", "hsa-miR-484", 
  "hsa-miR-143-3p", "hsa-miR-223-3p"
)

# Number of permutations
iterations <- 1000
p_value_threshold <- 0.05

# Simulate data
set.seed(42)  # For reproducibility

# Generate random expression data for ED, ward (10 samples each), and ICU (10 samples) patients
ed_data <- matrix(rnorm(10 * length(micrornas), mean = 10, sd = 2), nrow = 10, ncol = length(micrornas))
ward_data <- matrix(rnorm(10 * length(micrornas), mean = 10, sd = 2), nrow = 10, ncol = length(micrornas))
icu_data <- matrix(rnorm(10 * length(micrornas), mean = 12, sd = 2), nrow = 10, ncol = length(micrornas))
colnames(ed_data) <- micrornas
colnames(ward_data) <- micrornas
colnames(icu_data) <- micrornas

# Combine ED and ward data
ed_ward_data <- rbind(ed_data, ward_data)

# Create a matrix to store p-values
p_values <- matrix(NA, nrow = iterations, ncol = length(micrornas))
colnames(p_values) <- micrornas

# Perform permutations and t-tests
for (i in 1:iterations) {
  # Randomly select 10 ED/ward samples and 5 ICU samples without replacement
  selected_ed_ward <- ed_ward_data[sample(1:20, 10), ]
  selected_icu <- icu_data[sample(1:10, 5), ]
  
  # Perform t-tests for each miRNA
  for (j in 1:length(micrornas)) {
    t_test <- t.test(selected_ed_ward[, j], selected_icu[, j])
    p_values[i, j] <- t_test$p.value
  }
}

# Convert to data frame for easier manipulation
p_values_df <- as.data.frame(p_values)

# Calculate mean p-value and frequency of significance for each microRNA
mean_p_values <- colMeans(p_values_df)
significant_counts <- colSums(p_values_df < p_value_threshold)
summary_df <- data.frame(microRNA = names(mean_p_values), mean_p_value = mean_p_values, significant_count = significant_counts)

# Sort microRNAs based on mean p-value and frequency of significance
summary_df <- summary_df %>%
  arrange(mean_p_value, desc(significant_count))

# Create a long format data frame for boxplot
p_values_long <- stack(p_values_df)

# Generate the boxplot for all microRNAs sorted by mean p-value
sorted_micrornas <- summary_df$microRNA

p_boxplot <- ggplot(p_values_long, aes(x = values, y = factor(ind, levels = sorted_micrornas))) +
  geom_boxplot(outlier.size = 0.5, width = 0.3) + # Thinner boxes
  geom_vline(xintercept = p_value_threshold, linetype = "dashed", color = "black") + # Add threshold line
  labs(title = "Box Plot of P-Values for Each microRNA Across 1000 Permutations", x = "P-Value", y = "microRNAs") +
  theme_bw()

# Save the box plot
ggsave("boxplot_p_values_14micros.png", plot = p_boxplot, width = 12, height = 8)

print(p_boxplot)