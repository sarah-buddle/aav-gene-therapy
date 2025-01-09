library(tidyverse)
library(GenomicAlignments)

# Import alignment bam file
bam_raw <- readGAlignments(paste0(outdir, "/", sample, "_mapped.bam"), use.names = TRUE)

# Function to calculate longest match/mismatch segment from cigar string
longest_match_mismatch <- function(cigar) {
  matches <- str_extract_all(cigar, "\\d+[MX]")[[1]]
  lengths <- as.numeric(str_extract(matches, "\\d+"))
  max(lengths, na.rm = TRUE)
}

# Clean bam file
bam <- bam_raw %>%   
  as.data.frame(.) %>% 
  rownames_to_column(var = "qname") %>% 
  arrange(qname) %>% 
  separate_wider_delim(qname, names = c("qname", "n_alignment"), too_few = "align_start", delim = ".") %>% 
  mutate(n_alignment = ifelse(is.na(n_alignment), 1, as.numeric(n_alignment) + 1)) %>% 
  mutate(length_cigar = nchar(cigar)) %>% 
  mutate(longest_match_mismatch = sapply(cigar, longest_match_mismatch))
  
# Apply filters
filtered <- bam %>% 
  group_by(qname) %>% 
  summarise(total_aligned = sum(width), max_qwidth = max(qwidth), max_mm = max(longest_match_mismatch)) %>% 
  mutate(percent_aligned = total_aligned * 100 / max_qwidth) %>% 
  filter(percent_aligned >= 80 & max_qwidth >= 1000 & max_mm >= 100) %>% 
  select(qname)

# Output list of filtered reads
write_delim(filtered, file = paste0(outdir, "/", sample, "_mapped_filtered.txt"),
            col_names = FALSE, delim = "\t", quote = "none", escape = "none")
