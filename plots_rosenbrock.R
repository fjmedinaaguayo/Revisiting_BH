# Required libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Define the directory containing the .RData files
data_dir <- "./" # Change to your directory if needed

# Get a list of .RData files containing "mixmodel" in their names
rdata_files <- list.files(data_dir, pattern = "TEMP_rosenbrock_b9.*\\.RData$", full.names = TRUE)

# Initialize an empty list to store the data
all_data <- list()

# Define the methods corresponding to the three settings
methods <- c("RB", "BH", "AIS", "fr-SMC", "ar-SMC")

# Loop through each file and load the data
for (file in rdata_files) {
  # Load the variables from the .RData file
  load(file)
  
  # Extract metadata from the filename (e.g., "TEMP_mixmodel_N500_d2.RData")
  file_info <- strsplit(basename(file), "_|\\.")[[1]]
  N <- as.numeric(gsub("N", "", file_info[4]))
  d <- as.numeric(gsub("d", "", file_info[5]))
  
  # Prepare data frames for each metric
  logZ_df <- as.data.frame(logZ_vec) %>%
    mutate(Method = c("RB", "BH", "AIS", "fr-SMC", "ar-SMC")) %>%
    pivot_longer(cols = starts_with("V"), names_to = "Repetition", values_to = "logZ") %>%
    mutate(Repetition = as.integer(gsub("V", "", Repetition)))  # Convert repetition to integers
  
  time_df <- as.data.frame(time_vec) %>%
    mutate(Method = c("MIS", "AIS", "fr-SMC", "ar-SMC")) %>%
    pivot_longer(cols = starts_with("V"), names_to = "Repetition", values_to = "Time") %>%
    mutate(Repetition = as.integer(gsub("V", "", Repetition)))
  
  n_dists_df <- as.data.frame(n_dists_vec) %>%
    mutate(Method = c("AIS", "fr-SMC", "ar-SMC")) %>%
    pivot_longer(cols = starts_with("V"), names_to = "Repetition", values_to = "n_dists") %>%
    mutate(Repetition = as.integer(gsub("V", "", Repetition)))
  
  # Merge the three metrics into a single data frame
  merged_data <- logZ_df %>%
    left_join(time_df, by = c("Method", "Repetition")) %>%
    left_join(n_dists_df, by = c("Method", "Repetition")) %>%
    mutate(d = d, N = N, K = K)
  
  # Add the merged data to the list
  all_data[[length(all_data) + 1]] <- merged_data
}

# Combine all data frames into a single data frame
final_data <- bind_rows(all_data)

# Preview the combined data
head(final_data)




library(ggplot2)

# Filter data for a specific method (e.g., AIS)
method_data <- final_data %>% filter(d == 3)

method_data %>%
  mutate(Method = forcats::fct_relevel(Method, c("RB", "BH", "AIS", "fr-SMC", "ar-SMC"))) %>%
  ggplot(aes(x = Method, y = logZ)) +
  geom_boxplot(fill = "white", color = "black") +  # Ensures white boxes with black borders
  geom_hline(yintercept = lZ_true, linetype = "dashed") +  # True value as dashed line
  labs(x = "Method", y = "logZ") +
  coord_cartesian(ylim = c(-500, NA)) +
  theme_minimal()

method_data %>%
  mutate(Method = forcats::fct_relevel(Method, c("RB", "BH", "AIS", "fr-SMC", "ar-SMC"))) %>%
  ggplot(aes(x = Method, y = logZ)) +
  geom_boxplot(fill = "white", color = "black") +
  geom_hline(aes(yintercept = lZ_true, linetype = "True value for log(Z)")) +  # <--- Agrega leyenda
  scale_linetype_manual(values = c("True value for log(Z)" = "dashed")) +       # Define cómo se ve
  labs(x = "Method", y = "logZ", linetype = NULL) +                  # Oculta título de leyenda
  coord_cartesian(ylim = c(-500, NA)) +
  theme_minimal() +
  theme(
    legend.position = c(0.8, 0.1),
    legend.background = element_rect(fill = alpha("white", 0.0), color = NA),
    legend.key = element_blank()
  )

ggplot(df, aes(x = Method, y = logZ)) +
  geom_boxplot(fill = "white", color = "black") +
  geom_hline(aes(yintercept = lZ_true, linetype = "True value for log(Z)")) +  # add linetype mapping
  scale_linetype_manual(name = NULL, values = c("True value for log(Z)" = "dashed")) +
  labs(x = "", y = "logZ") +
  theme_minimal() +
  theme(
    legend.position = c(0.5, 0.1),  # x, y coordinates (inside top-right)
    legend.background = element_rect(fill = alpha("white", 0.0), color = NA),
    legend.key = element_blank()
  )

ggsave("boxplot_ex2_bw_all.png", width = 5, height = 4, dpi = 300, units = "in")

method_data %>% filter(Method %in% c("fr-SMC", "ar-SMC")) %>%
  mutate(Method = forcats::fct_relevel(Method, c("fr-SMC", "ar-SMC"))) %>%
  ggplot(aes(x = Method, y = logZ)) +
  geom_boxplot(fill = "white", color = "black") +  # Ensures white boxes with black borders
  geom_hline(yintercept = lZ_true, linetype = "dashed") +  # True value as dashed line
  labs(x = "Method", y = "logZ") +
  coord_cartesian(ylim = c(NA, NA)) +
  theme_minimal()

ggsave("boxplot_ex2_bw_sub.png", width = 2.5, height = 4, dpi = 300, units = "in")

# Create boxplots of logZ for the specified method across dimensions
final_data %>%
  mutate(logZ = case_when(logZ - quantile(logZ)[4] > 1.5*IQR(logZ) ~ NA_real_,
                          quantile(logZ)[2] - Time > 1.5*IQR(logZ) ~ NA_real_,
                          TRUE ~ logZ)) %>%
  ggplot(aes(x = factor(d), y = logZ, colour = Method)) +
  geom_boxplot() +
  labs(x = "dimension (d)", y = "log(Z)") +
  coord_cartesian(ylim=c(-1000, NA)) +
  theme_minimal()

ggplot(final_data, aes(x = factor(d), y = logZ, colour = Method)) +
  geom_boxplot() +
  labs(x = "dimension (d)", y = "log(Z)") +
  theme_minimal()

ggplot(final_data, aes(x = factor(d), y = n_dists, colour = Method)) +
  geom_boxplot() +
  labs(x = "dimension (d)", y = "num. inter. distns.") +
  theme_minimal()

ggplot(final_data, aes(x = factor(d), y = Time, colour = Method)) +
  geom_boxplot() +
  labs(x = "dimension (d)", y = "time") +
  theme_minimal()

final_data %>%
  mutate(Time = case_when(Time - quantile(Time)[4] > 1.5*IQR(Time) ~ NA_real_,
                          quantile(Time)[2] - Time > 1.5*IQR(Time) ~ NA_real_,
                          TRUE ~ Time)) %>%
  ggplot(aes(x = factor(d), y = Time, colour = Method)) +
  geom_boxplot() +
  labs(title = "Time across dimensions",
       x = "dimension (d)", y = "log(Z)") +
  theme_minimal()

