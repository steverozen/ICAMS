zz = ReadVCF(
  "~/MEGA/important_mut_sig_data/pcawg_indel_vcfs/388a8875-c3f5-494e-8456-28be8d3626e1.consensus.indel.vcf.gz",
  filter.status = NULL
)
file.remove("/tmp/out.out")
Rprof("/tmp/out.out")
uu = AnnotateIDVCF(zz, ref.genome = "hg19")
Rprof(NULL)
summaryRprof("/tmp/out.out")

# 0980e7fd-051d-45e9-9ca6-2baf073da4e8.consensus.indel.vcf.gz
#161K

zz = ReadVCF(
  "~/MEGA/important_mut_sig_data/pcawg_indel_vcfs/0980e7fd-051d-45e9-9ca6-2baf073da4e8.consensus.indel.vcf.gz",
  filter.status = NULL
)
file.remove("/tmp/out.out")
Rprof("/tmp/out.out", line.profiling = TRUE)
uu = AnnotateIDVCF(zz, ref.genome = "hg19")
Rprof(NULL)
summaryRprof("/tmp/out.out", lines = "show")
