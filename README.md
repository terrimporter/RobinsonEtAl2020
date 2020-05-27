# README

This repository contains the dataflow and scripts used to process the COI metabarcode reads in the paper Robinson et al., 2020 (submitted to Freshwater Science).

Infiles and scripts can be downloaded from https://github.com/terrimporter/RobinsonEtAl2020/releases .

## Infiles

cat.csv contains the ESV ids, sample names, and taxonomic assignments.  
Values_for_VennDiagrams.xlsx contains the values used in the Venn diagrams.
Sites.csv contains the coordinates needed to plot sites on the map.

## R Scripts

Fig1_Richness.R calcualtes ESV richness.  Uses cat.csv as an infile.

Fig2_NMDS.R creates NMDS plots.  Uses cat.csv as an infile.

Fig3_Venn.R draws Venn diagrams.  Uses values from Values_for_VennDiagrams.xlsx

Fig4_AntifreezeVsEthanol.R compares Ethanol vs Antifreeze xy-plot showing orders and genera.  Uses cat.csv as an infile.

Fig5_EPTheatmap.R compares Ethanol vs Antifreeze EPT families as heatmaps.  Uses cat.csv as an infile.

## Supplementary data files

1. Denoised, chimera-filtered ESVs (FASTA) for BR5, F230R, and ml-jg.  

2. Denoised, chimera-filtered ESVs, pseudogene-filtered ORFs (FASTA) for BR5, F230R, and ml-jg.

3. ESV x sample tables for BR5, F230R, and ml-jg.  These are based on the denoised, chimera-, and pseudogene-filtered ESVs.

4. FigS1_Map.R plots sites on the map.  Uses Sites.csv as an infile.

## References

Robinson, C.V., Porter, T.M., Wright, M.T., Hajibabaei, M. (2020) Propylene glycol-based antifreeze as an effective preservative for bulk-benthos DNA metabarcoding of macroinvertebrates.  Submitted to Freshwater Science.

## Acknowledgements

I would like to acknowledge funding from the Canadian government from the Genomics Research and Development Initiative (GRDI), Metagenomics-Based Ecosystem Biomonitoring (Ecobiomics) project.

Last updated: May 26, 2020
