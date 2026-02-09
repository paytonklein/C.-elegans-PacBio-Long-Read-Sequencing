
library(dplyr)

#input_file <- "/Users/jleighton32/Library/CloudStorage/OneDrive-Personal/C. elegans FCD-2/Long Read Sequencing/Sniffles/Sniffles_merged.tsv"
#output_file <- "/Users/jleighton32/Library/CloudStorage/OneDrive-Personal/C. elegans FCD-2/Long Read Sequencing/Sniffles/FilteredSVs_Sniffles.tsv"

df <- read.table(input_file,
                 header = TRUE,
                 sep = "\t",
                 stringsAsFactors = FALSE)

# all 5 wild type samples
samples_5 <- c(
  "1-bc2065","2-bc2066","3-bc2067","4-bc2068","5-bc2069"
)
df_filtered <- df %>%
  filter(SAMPLE %in% samples_5) %>%
  group_by(CHROM, START, END, SVLEN, SVTYPE) %>%
  mutate(n_samples = n_distinct(SAMPLE)) %>%
  ungroup() %>%
  filter(n_samples < length(samples_5)) %>%
  select(-n_samples)
write.table(df_filtered,
            file = output_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
cat("Done! Output saved to:\n", output_file, "\n")
