indir1 = "~/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs/unfiltered_vcfs"
indir2 = "~/MEGA/important_mut_sig_data/pcawg_indel_vcfs"

file1 = dir(indir1, pattern = "DRUP01030028T", full.names = TRUE)

vcf <- ICAMS::ReadVCFs(file1, filter.status = "PASS")
library(ICAMS)
F <-
  AnnotateIDVCF(
    vcf[[1]],
    ref.genome = "hg19",
    explain_indels = 1
  )
# doesn't work need the stupid uuids
file2 = dir(indir1, pattern = "CPCT02100089T", full.names = TRUE)
# file3 = dir(indir2, pattern = "SP23925", full.names = TRUE)
# file4 = dir(indir2, pattern = "SP48008", full.names = TRUE)

library(stringr)

F[[1]] %>%
  dplyr::filter(nchar(ins_or_del_seq) == 2) %>%
  dplyr::filter(ins_or_del == "d") %>%
  dplyr::count(short_visual) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::mutate(vis2 = gsub("\\}", "", gsub("\\{", "", short_visual))) %>%
  dplyr::mutate(dinuc = str_extract(vis2, "(?<=<).{2}(?=>)")) %>%
  dplyr::mutate(hasmh = grepl("\\{", short_visual)) -> foo

write.csv("foo", "del2_in_DRUP01030028T.csv")

F[[1]] %>%
  dplyr::filter(nchar(ins_or_del_seq) == 2) %>%
  dplyr::filter(ins_or_del == "d") -> zz

foo %>%
  group_by(dinuc, hasmh) %>%
  dplyr::summarize(total_n = sum(n), .groups = "drop") %>%
  dplyr::arrange(dinuc) -> bar

View(dplyr::filter(foo, dinuc %in% c("AG", "TC")))

vcf <- ICAMS::ReadVCFs(file2, filter.status = "PASS")
library(ICAMS)
list <-
  AnnotateIDVCF(
    vcf[[1]],
    ref.genome = "hg19",
    explain_indels = 1
  )

list[[1]] %>%
  dplyr::filter(nchar(ins_or_del_seq) > 1) %>%
  dplyr::filter(ins_or_del == "d") %>%
  dplyr::count(short_visual) %>%
  dplyr::arrange(desc(n)) -> foo4

write.csv("foo4", "del_in_CPCT02100089T.csv")
