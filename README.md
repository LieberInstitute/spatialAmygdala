

# Multiscale spatial transcriptomics resolves the cellular and molecular architecture of the human amygdala

## Overview

Welcome to the `spatialAmygdala` project! In this study, we generated
spatially-resolved transcriptomics (SRT) data from postmortem tissue sections of the human amygdala across  adult neurotypical donors. SRT data was generated using
[10x Genomics **Visium**](https://www.10xgenomics.com/products/spatial-gene-expression), [10x Genomics **Visium HD**](https://www.10xgenomics.com/products/visium-hd-spatial-gene-expression), and [10x Genomics **Xenium**](https://www.10xgenomics.com/platforms/xenium). 

Thank you for your interest in our work!

## Study design

![](code/Fig1_final.png)

TODO — full caption. Follow the dACC pattern: one summary sentence, then a
lettered walkthrough of each panel.

Experimental design to generate paired single-nucleus RNA-sequencing
(snRNA-seq) and spatially-resolved transcriptomics (SRT) data in the human
amygdala. **(A)** TODO — anatomical context; which nuclei/subregions, at what
level, from which donors. **(B)** TODO — tissue block, H&E validation, and
which assays were run from each block. **(C)** TODO — if you have a panel
showing the Xenium/Visium HD arm. (This figure was created with
[Biorender](https://biorender.com))

## Interactive Websites

All of these interactive websites are powered by the [`Vitessce`](http://vitessce.io/) open source software:

We provide the following interactive websites, organized by dataset with
software labeled by emojis:

- Visium (n = TODO)
  * [Visium](TODO_URL)
    + Provides interactive spot-level visualization of full-transcriptome profiles, spatial domains, 
    and predicted cell type locations.
- Xenium (n = TODO)
  * [Xenium](TODO_URL)
    + Provides interactive cell segmentations with a targeted gene panel, spatial domain, and transferred cell type labels.
- Visium HD (n = TODO)
  * [Visium HD](TODO_URL)
    + Provides interactive cell segmentations with full-transcriptome profliling, spatial domains, and 
        transferred cell type labels.
- Br9280 (n = TODO)
  * [Visium, VisiumHD, and Xenium](TODO_URL)
    + Provides interactive visualize of the same donor across all three technologies.


## Data Access

All data, including raw FASTQ files and `SpaceRanger` + `XeniumRanger`
processed data outputs, can be accessed via Gene Expression Omnibus (GEO) under
accessions [TODO_GSE](TODO_URL).

TODO — if the Zarr stores backing the Vitessce apps are hosted separately
(e.g. data.libd.org), add a line pointing there.

## Contact

We value feedback and public questions, as they allow other users to learn from the answers.
If you have any questions, please ask them at
[LieberInstitute/spatialAmygdala/issues](https://github.com/LieberInstitute/spatialAmygdala/issues)
and refrain from emailing us. Thank you again for your interest in our work!

## How to Cite

TODO_AUTHOR, et al. TODO_TITLE. BioRxiv (TODO_YEAR) doi:TODO_DOI.

## Internal

JHPCE location: `/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/`