bcftools filter \
--include 'MAF>=0.05' \
--regions-file /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/hg19knowngeneUCSCgenes_nochr.bed \
--IndelGap 0 \
-O v \
--output /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/allTestChromsRSid.reffix_MAF_0.05_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/allTestChromsRSid.reffix.vcf.gz

bcftools view -O v \
-s TEST_1,TEST_6,TEST_12,TEST_14,TEST_40,TEST_42,TEST_45,TEST_46 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/181010_lane3_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/allTestChromsRSid.reffix_MAF_0.05_exons.vcf.gz

bcftools view -O v \
-s TEST_15,TEST_17,TEST_23,TEST_25,TEST_28,TEST_32,TEST_39,TEST_43 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/181011_lane3_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/allTestChromsRSid.reffix_MAF_0.05_exons.vcf.gz

bcftools view -O v \
-s TEST_1,TEST_12,TEST_18,TEST_25,TEST_32,TEST_42,TEST_45,TEST_46 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/181017_lane3_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/allTestChromsRSid.reffix_MAF_0.05_exons.vcf.gz

bcftools view -O v \
-s TEST_17,TEST_25,TEST_28,TEST_40,TEST_47,TEST_50,TEST_53,TEST_55 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/181105_lane3_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/allTestChromsRSid.reffix_MAF_0.05_exons.vcf.gz

bcftools view -O v \
-s TEST_6,TEST_14,TEST_17,TEST_28,TEST_43,TEST_52,TEST_53,TEST_56 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/181106_lane3_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/allTestChromsRSid.reffix_MAF_0.05_exons.vcf.gz

bcftools view -O v \
-s TEST_1,TEST_6,TEST_18,TEST_23,TEST_45,TEST_47,TEST_55,TEST_56 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/181121_lane3_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/allTestChromsRSid.reffix_MAF_0.05_exons.vcf.gz

bcftools view -O v \
-s TEST_15,TEST_18,TEST_23,TEST_32,TEST_39,TEST_46,TEST_51,TEST_52 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/181122_lane3_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/allTestChromsRSid.reffix_MAF_0.05_exons.vcf.gz

bcftools view -O v \
-s TEST_15,TEST_42,TEST_43,TEST_47,TEST_50,TEST_51,TEST_55,TEST_56 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/181213_lane4_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/allTestChromsRSid.reffix_MAF_0.05_exons.vcf.gz

bcftools view -O v \
-s TEST_12,TEST_14,TEST_39,TEST_40,TEST_50,TEST_51,TEST_52,TEST_53 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/181214_lane3_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/allTestChromsRSid.reffix_MAF_0.05_exons.vcf.gz

bcftools filter \
--include 'MAF>=0.05' \
--regions-file /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/hg19knowngeneUCSCgenes_nochr.bed \
--IndelGap 0 \
-O v \
--output /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_v3_merged_20200624/stemi_v3_all_samples_20200624.RSid_MAF_0.05_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_v3_merged_20200624/stemi_v3_all_samples_20200624.RSid.vcf.gz

bcftools view -O v \
-s TEST_58,TEST_69,TEST_88,TEST_61,TEST_66,TEST_79,TEST_64,TEST_71 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/191125_lane1_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_v3_merged_20200624/stemi_v3_all_samples_20200624.RSid_MAF_0.05_exons.vcf

bcftools view -O v \
-s TEST_60,TEST_62,TEST_68,TEST_70,TEST_77,TEST_65,TEST_73,TEST_81 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/191125_lane2_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_v3_merged_20200624/stemi_v3_all_samples_20200624.RSid_MAF_0.05_exons.vcf

bcftools view -O v \
-s TEST_64,TEST_70,TEST_77,TEST_58,TEST_81,TEST_62,TEST_66,TEST_79 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/191126_lane1_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_v3_merged_20200624/stemi_v3_all_samples_20200624.RSid_MAF_0.05_exons.vcf

bcftools view -O v \
-s TEST_61,TEST_68,TEST_71,TEST_65,TEST_73,TEST_88,TEST_60,TEST_69 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/191126_lane2_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_v3_merged_20200624/stemi_v3_all_samples_20200624.RSid_MAF_0.05_exons.vcf

bcftools view -O v \
-s TEST_66,TEST_73,TEST_79,TEST_62,TEST_69,TEST_61,TEST_68,TEST_88 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/191209_lane1_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_v3_merged_20200624/stemi_v3_all_samples_20200624.RSid_MAF_0.05_exons.vcf

bcftools view -O v \
-s TEST_65,TEST_81,TEST_60,TEST_64,TEST_71,TEST_58,TEST_70,TEST_77 \
-o /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/191209_lane2_mmaf005_exons.vcf \
/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_v3_merged_20200624/stemi_v3_all_samples_20200624.RSid_MAF_0.05_exons.vcf
