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
	
	- Linear and/or logistic regressions using [bulk RNAseq data are here](https://rushalz.github.io/Cytokines_public_release/association_analysis/signif_heatmap_correl_bulkRNA.html)
	- Linear and/or logistic regressions using [single-nuclei RNAseq data are here](https://rushalz.github.io/Cytokines_public_release/association_analysis/signif_heatmap_correl_snRNA.html)
	- Scatter plot showing correlation among IL15 expression in the bulk RNAseq dataset and AD traits can be found [here](https://rushalz.github.io/Cytokines_public_release/association_analysis/scatter_plots_IL15_bulk.html)
	- Scatter plot showing correlation among IL15 microglial expression and AD pathological traits can be found [here](https://rushalz.github.io/Cytokines_public_release/association_analysis/scatter_plot_IL15_mic.html)
	- Boxplot showing IL15 expression in bulk RNAseq data by diagnosis can be found [here](https://rushalz.github.io/Cytokines_public_release/association_analysis/boxplot_IL15_bulk.html)
	- Boxplot showing IL15 expression in microglia by diagnosis can be found [here](https://rushalz.github.io/Cytokines_public_release/association_analysis/boxplot_IL15_mic.html)
	- Regression analysis with cytokine genes grouped by family: bulk data is [here](https://rushalz.github.io/Cytokines_public_release/association_analysis/family_bulk_LR_by_pheno.html) and single-nuclei data are [here](https://rushalz.github.io/Cytokines_public_release/association_analysis/family_sn_LR_by_pheno.html)
	
2. Biological sources of variance for the cytokines:

	- Variance Partition analysis using bulk RNAseq data are [here](https://rushalz.github.io/Cytokines_public_release/biological_drivers/vp_cytokines_bulk.html)
	- Variance Partition analysis using single-nuclei RNAseq data are [here](https://rushalz.github.io/Cytokines_public_release/biological_drivers/vp_cytokines_sn.html)
	
3. Genetic risk on cytokines expression:

	- AD-PRS analysis using single-nuclei RNAseq data can be found [here](https://rushalz.github.io/Cytokines_public_release/prs_analysis/prs_sn_association.html)
	
4. Gene set enrichment analysis (GSEA):

	- Bulk RNAseq [here](https://rushalz.github.io/Cytokines_public_release/gsea/gsea_bulk_RNAseq.html)
	- Single-nuclei RNAseq [here](https://rushalz.github.io/Cytokines_public_release/gsea/gsea_sn_RNAseq.html)