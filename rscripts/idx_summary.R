library(tidyverse)

idx_path <- "/fh/fast/warren_h/ebv_enktl/output/bowtie_align/bam/idxstats"
idx_files <- list.files(path = idx_path, 
    pattern = "*.txt",
    full.names = TRUE)


idx_reader <- function(idx_path) {
    sample_name <- str_extract(basename(idx_path), "EBV\\d+")
    idx_table <- read_tsv(idx_path,
        col_names = c("chrom", "length", "aligned_reads", "unaligned_reads")) %>%
        mutate(sample = sample_name)
    human_reads <- idx_table %>% 
        filter(chrom != "chrEBV") %>% 
        pull(aligned_reads) %>% 
        sum()
    ebv_reads <- idx_table %>% 
        filter(chrom == "chrEBV") %>% 
        pull(aligned_reads) %>% 
        sum()
    total_reads <- idx_table %>% 
        pull(aligned_reads) %>% 
        sum()
    idx_sum <- tibble(sample = sample_name, 
        human = human_reads, 
        ebv = ebv_reads,
        total = total_reads) 
    return(idx_sum)
}

idx_tables <- idx_files %>% 
    purrr::map(idx_reader) %>% 
    bind_rows() %>% 
    mutate(ebv_percent = (ebv/total) * 100 )


