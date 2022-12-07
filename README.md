# Raw data, processed data, and scripts for *Anti-diarrheal drug loperamide induces microbial dysbiosis in larval zebrafish via targeted bacterial inhibition*


This repository contains the scripts, processed sequencing artifacts, and the Rmd script files to reproduce the figures in the manuscript. The raw sequences generated for this study can be found in the NCBI Short Read Archive under BioProject no. PRJNA908751. The corresponding accession numbers for each 16S rRNA amplicon sample are detailed in [`PRJNA908751_NCBI_info.xlsx`](16SampliconAnalysis/metadata/PRJNA908751_NCBI_info.xlsx).

### To cite this work:
Stevick, R.J., Ghigo, J-M. & Pérez-Pascual, D. ....

### This repository has been archived on Zenodo. Access or cite the most recent release:  

----------------------------------------------------------------------  

# Contents

- [LoperamideStraininfo.xlsx](LoperamideStrainInfo.xlsx) - metadata for the strains used in the study
- [figures/](figures/) - general figures

## [16SampliconAnalysis](/16SampliconAnalysis)
This folder contains scripts used to perform QIIME2 analysis on the 16S rRNA amplicon data, the QIIME2 output artifacts, and the Rmd script to reproduce the figures in the manuscript.  
- [scripts/](16SampliconAnalysis/scripts) - QIIME2 bash script and Rmd files
- [output/](16SampliconAnalysis/output) - QIIME2 artifacts and processed files
- [metadata/](16SampliconAnalysis/metadata)
- [figures/](16SampliconAnalysis/figures)

## [inVitroAnalysis](/inVitroAnalysis)
This folder contains all the raw data files and the scripts to reproduce the figures and statistics in the manuscript.  
- [Growth curves](inVitroAnalysis/GrowthCurves)
- [Survival in water](inVitroAnalysis/WaterSurvival)  

## [inVivoAnalysis](/inVivoAnalysis)
This folder contains all the raw data files and the scripts to reproduce the figures and statistics in the manuscript.  
- [Mono-reconv CFUs](inVivoAnalysis/Mono)
- [Mix-reconv CFUs](inVivoAnalysis/Mix)
- [Fish growth and development](inVivoAnalysis/Growth_Development)

-------------------------------------------------------------------------

# Figure index


|  Fig | Description | Code | Figure |
| :---: | :------ | :--- | :----: |
| 1 | Experimental schematic | Made using BioRender | <img src="figures/Figure1_LoperamideExperimentalSchematic.png" width="300"> |
| 2 | Conventional microbiota changes | [Loperamide16Samplicon_analysis.Rmd](16SampliconAnalysis/scripts/Loperamide16Samplicon_analysis.Rmd#L605) | <img src="16SampliconAnalysis/figures/Figure2_LoperamideBetaAll_LimmaGenus.png" width="200"> |
| 3 | Bacterial survival | [WaterSurvivalCFUs_loperamide.Rmd](inVitroAnalysis/WaterSurvival/WaterSurvivalCFUs_loperamide_figure3.Rmd) | <img src="inVitroAnalysis/WaterSurvival/Figure3_WaterInVitro_LoperamideCFUs.png" width="200"> |
| 4 | Mono-colonization | [ZebrafishMonoCFUs_loperamide.Rmd](inVivoAnalysis/Mono/ZebrafishMonoCFUs_loperamide.Rmd#L116) | <img src="inVivoAnalysis/Mono/Figure4_LoperamideMonoColonization_withControl.png" width="200"> |
| 5 | Mix-colonization | [ZebrafishMixCFUs_loperamide.Rmd](inVivoAnalysis/Mix/ZebrafishMixCFUs_loperamide.Rmd#L102) | <img src="inVivoAnalysis/Mix/Figure5_MixReconvPlots.png" width=200> |
| 6 | Alpha-diversity comparison | [ZebrafishMixCFUs_loperamide.Rmd](inVivoAnalysis/Mix/ZebrafishMixCFUs_loperamide.Rmd#L180) | <img src="inVivoAnalysis/Mix/Figure6_AlphaDiversity_Conv_Mix5.png" width=200> |
|||||
| S1 | QC and rarefaction curves | [Loperamide16Samplicon_analysis.Rmd](16SampliconAnalysis/scripts/Loperamide16Samplicon_analysis.Rmd#L279) | <img src="16SampliconAnalysis/figures/FigureS1_16SLoperamideQC_rarecurve_reads.png" width="200"> |
| S2 | Sequencing controls  | [Loperamide16Samplicon_analysis.Rmd](16SampliconAnalysis/scripts/Loperamide16Samplicon_analysis.Rmd#L375)  | <img src="16SampliconAnalysis/figures/FigureS2_LoperamideQC_mocknegativecontrols.png" height="100"> |
| S3 | Phylum barplot  | [Loperamide16Samplicon_analysis.Rmd](16SampliconAnalysis/scripts/Loperamide16Samplicon_analysis.Rmd#L463) | <img src="16SampliconAnalysis/figures/FigureS3_LoperamideBarsPhylum.png" width="200"> |
| S4 | Genus barplot and venn diagram | [Loperamide16Samplicon_analysis.Rmd](16SampliconAnalysis/scripts/Loperamide16Samplicon_analysis.Rmd#L507) | <img src="16SampliconAnalysis/figures/FigureS4_LoperamideBarsGenusUpset.png" width="100"> |
| S5 | All samples beta-diversity | [Loperamide16Samplicon_analysis.Rmd](16SampliconAnalysis/scripts/Loperamide16Samplicon_analysis.Rmd#L740) | <img src="16SampliconAnalysis/figures/FigureS5_16Sloperamide_BetaDaysAndWithinDiv.png" width="100"> |
| S6 | Zebrafish growth | [LoperamideFishMeasurements.Rmd](inVivoAnalysis/Growth_Development/LoperamideFishMeasurents.Rmd) | <img src="inVivoAnalysis/Growth_Development/FigureS6_fishgrowth_loperamide.png" width="200"> |
| S7 | Bacterial growth | [LoperamideGrowthCurves.Rmd](inVitroAnalysis/GrowthCurves/LoperamideGrowthCurves_S10.Rmd) | <img src="inVitroAnalysis/GrowthCurves/FigureS7_GrowthCurvesLoperamide.png" width="200"> |
| S8 | Compare Mono vs Water | [ZebrafishMonoCFUs_loperamide.Rmd](inVivoAnalysis/Mono/ZebrafishMonoCFUs_loperamide.Rmd#L198) | <img src="inVivoAnalysis/Mono/FigureS8_WaterSurvivalFishMono.png" width="200"> |
| S9 | Compare Mix vs Mono | [ZebrafishMixCFUs_loperamide.Rmd](inVivoAnalysis/Mix/ZebrafishMixCFUs_loperamide.Rmd#L235) |<img src="inVivoAnalysis/Mix/FigureS9_CompareMonoMixReconv.png" width=200> |



-------------------------------------------------------------------------

![summary](figures/LoperamideGraphicalAbstract.png)
