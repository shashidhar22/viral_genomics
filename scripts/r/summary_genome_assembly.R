library(tidyverse)
library(dplyr)
library(purrr)
library(qwraps2)
library(gtsummary)

# for explanation on gtsummary tables see https://cran.r-project.org/web/packages/gtsummary/vignettes/tbl_summary.html

#Summary statistics genome assembly with genome during assembly
gen_ref_quast_output <- Sys.glob("./output_data/genome_assembly/quast/*/report.tsv")

comb_quast_outputs <- function(x) {
  read_tsv(x) %>% pivot_wider(names_from = Assembly, values_from = str_extract(x, "EBV\\d+_\\w+")) %>% mutate_if(is.character,as.numeric) %>% mutate(sample = str_extract(x, "EBV\\d+_\\w+"), .before = "# contigs (>= 0 bp)") 
} 

comb_gen_ref_quast_summary <- purrr::map(gen_ref_quast_output, comb_quast_outputs) %>% bind_rows() 
gen_ref_quast_summary_cond <-comb_gen_ref_quast_summary[, c("N50", "L50", "# contigs", "Genome fraction (%)")]
gen_ref_quast_summary_cond_table <- gen_ref_quast_summary_cond %>% tbl_summary() %>% modify_caption("**Assembly With Reference Genome Summary**")
gen_ref_quast_summary_cond_table

#Summary statistics genome assembly without reference genome during assembly
gen_no_ref_quast_output <- Sys.glob("./output_data/genome_assembly_no_ref/quast/*/report.tsv")

comb_gen_no_ref_quast_summary <- purrr::map(gen_no_ref_quast_output, comb_quast_outputs) %>% bind_rows() 
gen_no_ref_quast_output_summary_cond <-comb_gen_no_ref_quast_summary[, c("N50", "L50", "# contigs", "Genome fraction (%)")]
gen_no_ref_quast_summary_cond_table <- gen_no_ref_quast_output_summary_cond %>% tbl_summary() %>% modify_caption("**Assembly Without Reference Genome Summary**")
gen_no_ref_quast_summary_cond_table


  