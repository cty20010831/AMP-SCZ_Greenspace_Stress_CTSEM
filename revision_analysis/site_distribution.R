# Site Distribution Analysis
# Extract site from first two characters of Participant_ID

library(dplyr)
library(tidyr)

# Read the data
data <- read.csv("ndvi_stress_data.csv")

# Extract site (first 2 characters of Participant_ID)
data$Site <- substr(data$Participant_ID, 1, 2)

# Create mapping of site codes to institution names
site_mapping <- c(
  "OR" = "University of Oregon",
  "SF" = "University of California, San Francisco",
  "LA" = "University of California, Los Angeles",
  "IR" = "University of California, Irvine",
  "SD" = "University of California, San Diego",
  "WU" = "Washington University in St. Louis",
  "GA" = "University of Georgia",
  "NN" = "Northwestern University",
  "OH" = "Ohio State University",
  "PI" = "University of Pittsburgh",
  "NC" = "University of North Carolina at Chapel Hill",
  "UR" = "Rochester University",
  "TE" = "Temple University",
  "PA" = "University of Pennsylvania",
  "NL" = "The Zucker Hillside Hospital (Northwell Health)",
  "SI" = "Icahn School of Medicine at Mount Sinai",
  "YA" = "Yale University",
  "OL" = "Olin Neuropsychiatry Research Center / Hartford Hospital",
  "HA" = "Harvard University / Brigham & Women's Hospital",
  "BM" = "Beth Israel Deaconess Medical Center"
)

# Get unique participants per site
site_distribution <- data %>%
  select(Participant_ID, Site) %>%
  distinct() %>%
  group_by(Site) %>%
  summarise(N = n()) %>%
  mutate(Percentage = sprintf("%.2f", N / sum(N) * 100)) %>%
  arrange(desc(N))

# Map site codes to institution names
site_distribution$Institution <- site_mapping[site_distribution$Site]

# Add total row
total_row <- data.frame(
  Site = "Total",
  N = sum(site_distribution$N),
  Percentage = "100.00",
  Institution = "All Sites"
)

site_distribution_table <- rbind(site_distribution, total_row)

# Reorder columns for better readability
site_distribution_table <- site_distribution_table %>%
  select(Site, Institution, N, Percentage)

# Print the table
cat("\nSite Distribution Analysis\n")
cat("==========================\n\n")
print(site_distribution_table, row.names = FALSE)

# Save the table
write.csv(site_distribution_table,
          "revision_analysis/site_distribution_table.csv",
          row.names = FALSE)

cat("\nTable saved to: site_distribution_table.csv\n")
