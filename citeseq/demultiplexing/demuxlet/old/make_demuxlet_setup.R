library(dplyr)
library(readr)
library(tibble)

citeseq_libs <- 
    c("pilot2" = "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/CITEseq_pilot_2/SN0257788/broad/hptmp/curtism/bwh10x/KW10170_Maria/220617_10X_KW10170_bcl/cellranger-6.1.1/GRCh38/BRI-1743_hashing/outs",
      "pilot2_reseq" = "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/CITEseq_pilot_2/SN0264956/broad/hptmp/sgurajal/bwh10x/KW10275_mgutierrez/220909_10X_KW10275_bcl/cellranger-6.1.1/GRCh38/BRI-1743_hashing/outs",
      "1984" = "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq/SN0268787/broad/hptmp/curtism/bwh10x/KW10456_Maria/221101_10X_KW10456_bcl/cellranger-6.1.1/GRCh38/BRI-1984/outs",
      "1988" = "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq/SN0273471/broad/hptmp/sgurajal/bwh10x/KW10598_mgutierrez/221211_10X_KW10598_bcl/cellranger-6.1.1/GRCh38/BRI-1988_hashing/outs",
      "1990" = "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq/SN0273471/broad/hptmp/sgurajal/bwh10x/KW10598_mgutierrez/221211_10X_KW10598_bcl/cellranger-6.1.1/GRCh38/BRI-1990_hashing/outs")  

vcfs <- 
    c("pilot2" = "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/citeseq/data/samples_23andMe_allsnps.vcf.gz",
     "pilot2_reseq" = "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/citeseq/data/samples_23andMe_allsnps.vcf.gz", 
     "1984" = "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/citeseq/data/allchr.mgb.generegions.vcf.gz",
     "1988" = "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/citeseq/data/allchr.mgb.generegions.vcf.gz",
     "1990" = "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/citeseq/data/allchr.mgb.generegions.vcf.gz")

setup <- 
    left_join(enframe(citeseq_libs, name = "lib_id", value = "scdata"),
	      enframe(vcfs, name = "lib_id", value = "vcf"))

write_tsv(setup, "demuxlet_setup.tsv", col_names = FALSE)
