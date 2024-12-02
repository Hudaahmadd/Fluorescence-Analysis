# ============================
# Fluorescence Intensity Analysis
# Author: Huda Ahmad
# Date: [2024-12-03]
# Purpose: Analyze fluorescence intensities for glucose and arabinose conditions in a 384-well plate.
# ============================

# ============================
# 1. Required Libraries
# ============================
library(imager)      # For image processing
library(magick)      # For advanced image handling
library(reshape2)    # For reshaping data
library(ggplot2)     # For data visualization
library(ggbeeswarm)  # For enhanced swarm plots

# ============================
# 2. Data Input
# ============================
# Set input and output paths
input_path <- "data/"      # Replace with your input folder
output_path <- "results/"  # Replace with your output folder

# Specify file names
glucose_file <- paste0(input_path, "Glucose.png.grid.jpg")
arabinose_file <- paste0(input_path, "Arabinose.png.grid.jpg")

# Check if files exist
if (!file.exists(glucose_file) || !file.exists(arabinose_file)) {
  stop("Image files not found. Please check the input path.")
}

# Load images
glucose_img <- load.image(glucose_file)
arabinose_img <- load.image(arabinose_file)

# ============================
# 3. Data Processing
# ============================
# Resize images to consistent dimensions
glucose_img_resized <- resize(glucose_img, 1115, 740)
arabinose_img_resized <- resize(arabinose_img, 1115, 740)

# Define plate dimensions
n_rows <- 16
n_cols <- 24
row_height <- dim(glucose_img_resized)[1] / n_rows
col_width <- dim(glucose_img_resized)[2] / n_cols

# Initialize a data frame for fluorescence intensity data
intensity_data <- data.frame(Row = integer(), Column = integer(),
                             Glucose_Intensity = numeric(), Arabinose_Intensity = numeric())

# Loop through wells to extract mean intensities
for (i in 1:n_rows) {
  for (j in 1:n_cols) {
    x_min <- round((j - 1) * col_width) + 1
    x_max <- round(j * col_width)
    y_min <- round((i - 1) * row_height) + 1
    y_max <- round(i * row_height)
    
    glucose_intensity <- mean(as.numeric(glucose_img_resized[y_min:y_max, x_min:x_max]))
    arabinose_intensity <- mean(as.numeric(arabinose_img_resized[y_min:y_max, x_min:x_max]))
    
    intensity_data <- rbind(intensity_data, data.frame(Row = i, Column = j,
                                                       Glucose_Intensity = glucose_intensity,
                                                       Arabinose_Intensity = arabinose_intensity))
  }
}

# Save extracted data to a CSV file
write.csv(intensity_data, paste0(output_path, "intensity_results.csv"), row.names = FALSE)

# ============================
# 4. Visualization
# ============================
# Convert data to long format for plotting
intensity_long <- melt(intensity_data, id.vars = c("Row", "Column"),
                       measure.vars = c("Glucose_Intensity", "Arabinose_Intensity"),
                       variable.name = "Condition", value.name = "Intensity")
intensity_long$Condition <- gsub("Glucose_Intensity", "Glucose", intensity_long$Condition)
intensity_long$Condition <- gsub("Arabinose_Intensity", "Arabinose", intensity_long$Condition)

# Create and save plots
# Swarm Plot
swarm_plot <- ggplot(intensity_long, aes(x = Condition, y = Intensity, color = Condition)) +
  geom_quasirandom() +
  theme_minimal() +
  labs(x = "Condition", y = "Fluorescence Intensity") +
  scale_color_manual(values = c("Glucose" = "blue", "Arabinose" = "red"))
ggsave(paste0(output_path, "swarm_plot.png"), swarm_plot, dpi = 300)

# Violin Plot
violin_plot <- ggplot(intensity_long, aes(x = Condition, y = Intensity, fill = Condition)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(x = "Condition", y = "Fluorescence Intensity") +
  scale_fill_manual(values = c("Glucose" = "lightblue", "Arabinose" = "lightcoral"))
ggsave(paste0(output_path, "violin_plot.png"), violin_plot, dpi = 300)

# Box Plot
box_plot <- ggplot(intensity_long, aes(x = Condition, y = Intensity, fill = Condition)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  theme_minimal() +
  labs(x = "Condition", y = "Fluorescence Intensity") +
  scale_fill_manual(values = c("Glucose" = "lightblue", "Arabinose" = "lightcoral"))
ggsave(paste0(output_path, "box_plot.png"), box_plot, dpi = 300)

# Density Plot
density_plot <- ggplot(intensity_long, aes(x = Intensity, fill = Condition)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(x = "Fluorescence Intensity", y = "Density") +
  scale_fill_manual(values = c("Glucose" = "lightblue", "Arabinose" = "lightcoral"))
ggsave(paste0(output_path, "density_plot.png"), density_plot, dpi = 300)

# ============================
# 5. Statistical Analysis
# ============================
# Perform a t-test
t_test_result <- t.test(Intensity ~ Condition, data = intensity_long)

# Print results to console
cat("T-test results:\n")
print(t_test_result)

# Save the results to a text file
stat_results_path <- paste0(output_path, "statistical_analysis_results.txt")
capture.output({
  cat("T-test results:\n")
  print(t_test_result)
}, file = stat_results_path)

cat("Statistical analysis results saved to:", stat_results_path, "\n")

# ============================
# End of Script
# ============================

