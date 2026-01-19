indir1 = "~/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs/unfiltered_vcfs"
file1 = dir(indir1, pattern = "DRUP01030028T", full.names = TRUE)

vcf <- ICAMS::ReadVCFs(file1, filter.status = "PASS")

list <-
  AnnotateIDVCF(
    vcf[[1]],
    ref.genome = "hg19",
    explain_indels = 1
  )

library(stringr)

list[[1]] %>%
  dplyr::filter(nchar(ins_or_del_seq) == 2) %>%
  dplyr::filter(ins_or_del == "d") %>%
  dplyr::count(short_visual) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::mutate(vis2 = gsub("\\}", "", gsub("\\{", "", short_visual))) %>%
  dplyr::mutate(dinuc = str_extract(vis2, "(?<=<).{2}(?=>)")) %>%
  dplyr::mutate(hasmh = grepl("\\{", short_visual)) -> foo

foo %>%
  group_by(dinuc, hasmh) %>%
  dplyr::summarize(total_n = sum(n), .groups = "drop") %>%
  dplyr::arrange(dinuc) -> bar
