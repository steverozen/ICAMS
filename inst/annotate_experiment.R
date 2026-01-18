indir1 = "~/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs/unfiltered_vcfs"
file1 = dir(indir1, pattern = "DRUP01030028T", full.names = TRUE)

vcf <- ICAMS::ReadVCFs(file1, filter.status = "PASS")

list <-
  AnnotateIDVCF(
    vcf[[1]],
    ref.genome = "hg19",
    explain_indels = 1
  )
