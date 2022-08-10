# GreaterGlider
A pipeline for sequence simulation (INDELible), ML tree reconstruction (IQTree), and site-concordance branch support measurement.

GreaterGlider.pl is a Perl written software pipeline for linux operating systems, allowing sequence simulation (INDELible) combined with ML tree reconstruction (IQTree) and subsequent conducted branch support estimation using the site-concordance (sCF) measure (IQTree2). To enable all three pipeline stages (simulation, ML recosntruction, and branch support measure) and thus a re-analysis of the original sCF study of Kück et al 2022, the pipeline consists a reduced source code of the more comprehensive, module based Appetite-pipeline (originally used in Kück et al. 2022). Processing greaterGlider.pl in default mode (perl greaterGlider <enter>) conducts a re-analysis of the ten tree-setups /T1 to T10) as published in Kück et al. 2022 with alpha=0.5. For other alpha values used by Kück et al. 2022, change the respective parameter in the script.

!IMPORTANT!: To use greaterGlider.pl following software must be pre-installed:
	:: INDELible v1.03 
		... Fletcher W. and Yang Z. 2009. INDELible: a flexible simulator of biological sequence evolution. Mol. Biol. and Evol. 2009 26(8):1879-1888
	:: IQTree v1.6.12 
		... Nguyen et al. 2015. IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. Mol. Biol. and Evol. 2015            32(1):268-274
	:: IQTree v2.2.0
		... Minh et al. 2020. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Mol. Biol. and Evol. 2015 37(5):1530-1534
	:: pdfunite
		... part of the TeX Live package (https://tug.org/texlive/) 
