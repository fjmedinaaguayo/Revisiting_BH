# Required libraries
library(dplyr)
library(tidyr)

# Define the directory containing the .RData files
data_dir <- "./" # Change to your directory if needed

# Get a list of .RData files containing "mixmodel" in their names
rdata_files <- list.files(data_dir, pattern = "TEMP_mixmodel_b99.*\\.RData$", full.names = TRUE)

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
    inner_join(time_df, by = c("Method", "Repetition")) %>%
    inner_join(n_dists_df, by = c("Method", "Repetition")) %>%
    mutate(d = d, N = N, K = K)
  
  # Add the merged data to the list
  all_data[[length(all_data) + 1]] <- merged_data
}

# Combine all data frames into a single data frame
final_data <- bind_rows(all_data) %>% mutate(k=d/2)

# Preview the combined data
head(final_data)




library(ggplot2)

# Filter data for a specific method (e.g., AIS)
method_data <- final_data %>% filter(Method == "AIS")

# Create boxplots of logZ for the specified method across dimensions
ggplot(final_data, aes(x = factor(k), y = logZ, colour = Method)) +
  geom_boxplot() +
  labs(x = "k", y = "logZ") +
  theme_minimal()

ggplot(final_data, aes(x = factor(k), y = logZ, fill = Method)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("white", "gray85", "gray50")) +
  labs(x = "k", y = "logZ") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("boxplot_logZ_ex3_bw.png", width = 10, height = 4, dpi = 300, units = "in")

ggplot(final_data, aes(x = factor(k), y = n_dists, colour = Method)) +
  geom_boxplot() +
  labs(x = "k", y = "num. inter. distns.") +
  theme_minimal()

ggplot(final_data, aes(x = factor(k), y = n_dists, fill = Method)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("white", "gray85", "gray50")) +
  labs(x = "k", y = "num. inter. distr.") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("boxplot_nid_ex3_bw.png", width = 10, height = 4, dpi = 300, units = "in")

ggplot(final_data, aes(x = factor(k), y = Time, colour = Method)) +
  geom_boxplot() +
  labs(x = "k", y = "time") +
  theme_minimal()

ggplot(final_data, aes(x = factor(k), y = Time, fill = Method)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = c("white", "gray85", "gray50")) +
  labs(x = "k", y = "time (secs.)") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("boxplot_time_ex3_bw.png", width = 10, height = 4, dpi = 300, units = "in")


final_data %>%
  mutate(Time = case_when(Time - quantile(Time)[4] > 1.5*IQR(Time) ~ NA_real_,
                                  quantile(Time)[2] - Time > 1.5*IQR(Time) ~ NA_real_,
                                  TRUE ~ Time)) %>%
  ggplot(aes(x = factor(k), y = Time, colour = Method)) +
  geom_boxplot() +
  labs(x = "k", y = "time") +
  theme_minimal()
