###############################################
# FLUORESCENCE PIPELINE (TWO OUTPUTS)
#
# Output folders (per run_name):
#  1) fluorescence_plots
#  2) fluorescence_plots_with_stats
#     - Same plot + Tukey HSD brackets (only p.adj < 0.05 shown)
#
# Inputs:
#  - Plate_info.txt: row, column, strain (and optionally Condition)
#  - fluorescence.csv: Row/row, Column/column, fluorescence OR intensity (and optionally Condition)
#
# Outputs:
#  - fluorescence_plots/Flu_<Condition>.tiff
#  - fluorescence_plots_with_stats/Flu_withStats_<Condition>.tiff
#  - CSV_outputs/fluorescence_data.csv
#  - CSV_outputs/TukeyHSD_Flu_<Condition>.csv
###############################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggbeeswarm)
  library(rstatix)
  library(ggpubr)
})

# ============================================================
# USER INPUT / OUTPUT (EDIT ONLY THIS PART)
# ============================================================

flu_path     <- "data/plate_reader/fluorescence.csv"
plate_path   <- "data/Plate_info/Plate_info.txt"
base_out_dir <- "outputs"

run_name <- "RUN_fluorescence"

# If fluorescence.csv does NOT contain Condition, set a default (e.g. "TSA")
default_condition_if_missing <- NA_character_

# Optional: labels for nicer plot titles (add more if needed)
condition_labels <- c("TSA" = "TSA")

# ============================================================
# DO NOT EDIT BELOW THIS LINE
# ============================================================

# ----------------------------
# Step 0 — Helper functions
# ----------------------------
save_tiff_plot <- function(plot_obj, filename, width_mm = 180, height_mm = 140, res = 300) {
  grDevices::tiff(
    filename = filename,
    width = width_mm, height = height_mm, units = "mm",
    res = res, compression = "lzw"
  )
  print(plot_obj)
  grDevices::dev.off()
}

stars_from_p <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

stack_brackets <- function(stat_df, y_vals, step_frac = 0.10) {
  ymax <- max(y_vals, na.rm = TRUE)
  yrng <- diff(range(y_vals, na.rm = TRUE))
  if (!is.finite(yrng) || yrng == 0) yrng <- 1
  
  stat_df %>%
    dplyr::arrange(p.adj) %>%
    dplyr::mutate(y.position = ymax + dplyr::row_number() * step_frac * yrng)
}

theme_pub <- function() {
  theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold")
    )
}

pretty_condition <- function(cond) condition_labels[[cond]] %||% cond

# ----------------------------
# Step 1 — Create output folders
# ----------------------------
out_dir_final <- file.path(base_out_dir, "plate_reader_data")
run_dir <- file.path(out_dir_final, run_name)

out_dir_plots <- file.path(run_dir, "fluorescence_plots")
out_dir_plots_stats <- file.path(run_dir, "fluorescence_plots_with_stats")
out_dir_csv <- file.path(run_dir, "CSV_outputs")

dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plots_stats, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_csv, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Step 2 — Read Plate_info.txt and detect ALL strain/OD levels
# ----------------------------
if (!file.exists(plate_path)) stop("plate_path not found: ", plate_path)

plate_info <- read.table(plate_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  dplyr::rename(Row = row, Column = column) %>%   # explicit dplyr avoids rename conflicts
  dplyr::mutate(
    Row = as.character(Row),
    Column = as.integer(Column),
    strain = format(as.numeric(strain), trim = TRUE, scientific = FALSE)
  )

all_strains_chr <- plate_info %>%
  dplyr::distinct(strain) %>%
  dplyr::mutate(strain_num = suppressWarnings(as.numeric(strain))) %>%
  dplyr::arrange(strain_num) %>%
  dplyr::pull(strain)

message("Detected OD/strain levels in plate_info: ", paste(all_strains_chr, collapse = ", "))

order_strain <- function(x) factor(as.character(x), levels = all_strains_chr)

# ----------------------------
# Step 3 — Read fluorescence and ensure required columns exist
# ----------------------------
if (!file.exists(flu_path)) stop("flu_path not found: ", flu_path)

flu_df <- read.csv(flu_path) %>%
  # accept both Row/row and Column/column
  dplyr::rename_with(~"Row", matches("^row$", ignore.case = TRUE)) %>%
  dplyr::rename_with(~"Column", matches("^col(umn)?$", ignore.case = TRUE)) %>%
  # accept either 'fluorescence' or 'intensity' and standardise to 'Intensity'
  dplyr::rename_with(~"Intensity", matches("^(fluorescence|intensity)$", ignore.case = TRUE)) %>%
  dplyr::mutate(
    Row = as.character(Row),
    Column = as.integer(Column)
  )

required_cols <- c("Row", "Column", "Intensity")
missing_cols <- setdiff(required_cols, names(flu_df))
if (length(missing_cols) > 0) {
  stop(
    "fluorescence.csv is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    "\nExpected: Row/row, Column/column, and fluorescence OR intensity."
  )
}

# If Condition is missing, use a default
if (!("Condition" %in% names(flu_df))) {
  if (is.na(default_condition_if_missing)) {
    stop(
      "fluorescence.csv has no 'Condition' column.\n",
      "Either add it OR set default_condition_if_missing (e.g. 'TSA')."
    )
  }
  flu_df <- flu_df %>% dplyr::mutate(Condition = default_condition_if_missing)
}

# ----------------------------
# Step 4 — Join strain/OD onto fluorescence using Row+Column
# ----------------------------
flu_for_plot <- flu_df %>%
  dplyr::left_join(plate_info, by = c("Row", "Column")) %>%
  dplyr::select(Condition, Row, Column, strain, Intensity) %>%
  dplyr::mutate(
    Condition = as.character(Condition),
    strain = order_strain(strain),
    Intensity = as.numeric(Intensity)
  ) %>%
  dplyr::filter(!is.na(strain), is.finite(Intensity), !is.na(Condition))

if (nrow(flu_for_plot) == 0) {
  stop("No valid rows after join. Check Row/Column matching and strain mapping in Plate_info.")
}

# ----------------------------
# Step 5 — Save the exact fluorescence data used for plotting
# ----------------------------
write.csv(
  flu_for_plot,
  file.path(out_dir_csv, "fluorescence_data.csv"),
  row.names = FALSE
)

# ----------------------------
# Step 6 — Plot builder (shared by both plot types)
# ----------------------------
build_plot <- function(dfc) {
  ggplot(dfc, aes(x = strain, y = Intensity, color = strain)) +
    geom_beeswarm(size = 1.6, alpha = 0.85, cex = 1.1) +
    stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 0.6) +
    stat_summary(aes(group = 1), fun = median, geom = "point", size = 2) +
    scale_x_discrete(drop = FALSE) +
    labs(x = "OD", y = "Fluorescence (Intensity)") +
    theme_pub()
}

# ----------------------------
# Step 7 — Fluorescence plots (NO stats)
# ----------------------------
plot_flu <- function(cond) {
  dfc <- flu_for_plot %>%
    dplyr::filter(Condition == cond) %>%
    dplyr::mutate(strain = order_strain(strain))
  
  p <- build_plot(dfc) +
    labs(title = paste0(pretty_condition(cond), " (Typhoon Fluorescence)"))
  
  save_tiff_plot(p, file.path(out_dir_plots, paste0("Flu_", cond, ".tiff")))
}

# ----------------------------
# Step 8 — Fluorescence plots (WITH Tukey HSD stats)
#         - Only p.adj < 0.05 shown on the plot
#         - Full Tukey table saved to CSV
# ----------------------------
plot_flu_with_stats <- function(cond) {
  dfc <- flu_for_plot %>%
    dplyr::filter(Condition == cond) %>%
    dplyr::mutate(strain = order_strain(strain))
  
  if (nlevels(droplevels(dfc$strain)) < 2) return(invisible(NULL))
  
  tuk <- rstatix::tukey_hsd(dfc, Intensity ~ strain) %>%
    dplyr::mutate(p.signif = stars_from_p(p.adj)) %>%
    stack_brackets(y_vals = dfc$Intensity, step_frac = 0.10)
  
  write.csv(
    tuk,
    file.path(out_dir_csv, paste0("TukeyHSD_Flu_", cond, ".csv")),
    row.names = FALSE
  )
  
  tuk_used <- tuk %>%
    dplyr::filter(p.adj < 0.05) %>%
    dplyr::select(group1, group2, p.adj, p.signif, y.position)
  
  p <- build_plot(dfc) +
    labs(title = paste0(pretty_condition(cond), " (Typhoon Fluorescence)"))
  
  if (nrow(tuk_used) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      tuk_used,
      label = "p.signif",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      size = 4
    )
  }
  
  save_tiff_plot(p, file.path(out_dir_plots_stats, paste0("Flu_withStats_", cond, ".tiff")))
}

# ----------------------------
# Step 9 — Run for all conditions
# ----------------------------
conds <- sort(unique(flu_for_plot$Condition))

purrr::walk(conds, plot_flu)
purrr::walk(conds, plot_flu_with_stats)

# ----------------------------
# Step 10 — Save run notes
# ----------------------------
writeLines(
  paste0(
    "Run: ", run_name, "\n",
    "Strains/ODs detected and included: ", paste(all_strains_chr, collapse = ", "), "\n",
    "Conditions plotted: ", paste(conds, collapse = ", "), "\n",
    "Output dir: ", normalizePath(run_dir), "\n"
  ),
  con = file.path(out_dir_csv, "run_notes.txt")
)

cat("\nALL DONE.\nSaved in:\n", normalizePath(run_dir), "\n", sep = "")
