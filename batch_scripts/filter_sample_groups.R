library(dplyr)
library(tidyr)
library(stringr)

input_vcf <- "/work/pi_nhowlett_uri_edu/jessie/New-All-20-Bam/pbsv_variants/all_samples.pbsv.vcf"
output_file <- "/work/pi_nhowlett_uri_edu/jessie/New-All-20-Bam/pbsv_variants/Filtered.tsv"

# -----------------------------
# STEP 1: Read VCF
# -----------------------------

# Read header line separately
vcf_header <- readLines(input_vcf)
header_line <- vcf_header[grep("^#CHROM", vcf_header)]

col_names <- strsplit(header_line, "\t")[[1]]

# Read VCF data (skip metadata lines)
df <- read.table(input_vcf,
                 header = FALSE,
                 sep = "\t",
                 comment.char = "#",
                 stringsAsFactors = FALSE)

colnames(df) <- col_names

# -----------------------------
# STEP 2: Extract needed fields
# -----------------------------

df <- df %>%
  rename(CHROM = `#CHROM`,
         START = POS)

# Extract END and SVTYPE from INFO
df$END <- str_extract(df$INFO, "END=\\d+") %>% str_remove("END=")
df$SVTYPE <- str_extract(df$INFO, "SVTYPE=[^;]+") %>% str_remove("SVTYPE=")

df$END <- as.numeric(df$END)
df$START <- as.numeric(df$START)

# -----------------------------
# STEP 3: Convert to long format
# -----------------------------

sample_cols <- col_names[10:length(col_names)]

df_long <- df %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "SAMPLE",
               values_to = "GENOTYPE")

# Keep only variants present (not 0/0)
df_long <- df_long %>%
  filter(str_detect(GENOTYPE, "0/1|1/1"))

# -----------------------------
# STEP 4: Define sample groups
# -----------------------------

n2_nt  <- c("1-bc2065","2-bc2066","3-bc2067","4-bc2068","5-bc2069")
fcd2_nt <- c("6-bc2070","7-bc2071","8-bc2072","9-bc2073","10-bc2074")

n2_hu  <- c("11-bc2075","12-bc2076","13-bc2077","14-bc2078","15-bc2079")
fcd2_hu <- c("16-bc2080","17-bc2081","18-bc2082","19-bc2083","20-bc2084")

# -----------------------------
# STEP 5: Filtering function (same logic)
# -----------------------------

remove_shared_nt <- function(data, nt_samples, hu_samples) {
  
  shared_nt <- data %>%
    filter(SAMPLE %in% nt_samples) %>%
    group_by(CHROM, START, END, SVTYPE) %>%
    summarise(n_samples = n_distinct(SAMPLE), .groups = "drop") %>%
    filter(n_samples == length(nt_samples)) %>%
    select(-n_samples)
  
  cleaned <- data %>%
    filter(SAMPLE %in% c(nt_samples, hu_samples)) %>%
    anti_join(shared_nt,
              by = c("CHROM","START","END","SVTYPE"))
  
  return(cleaned)
}

# -----------------------------
# STEP 6: Apply filtering
# -----------------------------

df_n2_clean   <- remove_shared_nt(df_long, n2_nt, n2_hu)
df_fcd2_clean <- remove_shared_nt(df_long, fcd2_nt, fcd2_hu)

df_final <- bind_rows(df_n2_clean, df_fcd2_clean)

# -----------------------------
# STEP 7: Save output
# -----------------------------

write.table(df_final,
            file = output_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat("Done! Baseline SVs (5/5 NT) removed per genotype.\n")