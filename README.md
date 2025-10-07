<h1 align="center"> Cytokine expression profile in the human brain of older adults </h1><br>

> This repository includes code and plots. Exploratory analysis and intermediate processing files are too large for this repository.

<p align="center">
 <img src="https://github.com/rushalz/Cytokines_public_release/blob/main/Fig1_Graph_abstract_Cytokines_2025.png">
</p>

**Figure 1:** Study design overview.

Exploratory analysis can be found [here](https://rushalz.github.io/Cytokines_public_release/exploratory_analysis/exploratory_plots_bulk.html)

List of cytokines with family annotations can be found [here](https://rushalz.github.io/Cytokines_public_release/cytokines_families.html)

Cytokines expression per dataset [here](https://rushalz.github.io/Cytokines_public_release/upset_plot_expr_genes.html)

Results:

1. Association analysis:
	
	- Association analysis using ROSMAP data (Figure 2 and 4c-d) are [here](https://rushalz.github.io/Cytokines_public_release/association_analysis/signif_heatmap_correl.html)
	- Regression analysis with cytokine genes grouped by family: [bulk data is here](https://rushalz.github.io/Cytokines_public_release/association_analysis/family_bulk_LR_by_pheno.html) and [single-nuclei data are here](https://rushalz.github.io/Cytokines_public_release/association_analysis/family_sn_LR_by_pheno.html)
	- Replication analysis using MSBB bulk RNA-Seq data is [here](https://rushalz.github.io/Cytokines_public_release/association_analysis/lr_expr_bulk_MSBB.html)
	- Replication analysis using Mayo Clinic bulk RNA-Seq data is [here](https://rushalz.github.io/Cytokines_public_release/association_analysis/lr_expr_bulk_Mayo.html)
	
2. Biological sources of variance for the cytokines:

	- Variance Partition analysis using bulk RNAseq data are [here](https://rushalz.github.io/Cytokines_public_release/biological_drivers/vp_cytokines_bulk.html)
	- Variance Partition analysis using single-nuclei RNAseq data are [here](https://rushalz.github.io/Cytokines_public_release/biological_drivers/vp_cytokines_sn.html)
	
3. Genetic risk on cytokines expression:

	- AD genetic risk (AD-PRS and APOE) and their associations with cytokine brain expression in specific cell types can be found [here](https://rushalz.github.io/Cytokines_public_release/APOE_and_PRS_results/PRS_analysis_and_figures_review1.html)
	- Results of linear and/or logistic regressions between cytokines and [APOE ε4](https://rushalz.github.io/Cytokines_public_release/APOE_and_PRS_results/signif_heatmap_APOE_tests_review) and [AD-PRS](https://rushalz.github.io/Cytokines_public_release/APOE_and_PRS_results/signif_heatmap_PRS_tests_review)
	
4. Gene set enrichment analysis (GSEA):

	- Bulk RNAseq [here](https://rushalz.github.io/Cytokines_public_release/gsea/gsea_bulk_RNAseq.html)
	- Single-nuclei RNAseq [here](https://rushalz.github.io/Cytokines_public_release/gsea/gsea_sn_RNAseq.html)