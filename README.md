

# Multiscale spatial transcriptomics resolves the cellular and molecular architecture of the human amygdala

## Overview

Welcome to the `spatialAmygdala` project! In this study, we generated
spatially-resolved transcriptomics (SRT) data from postmortem tissue sections of the human amygdala across  adult neurotypical donors. SRT data was generated using
[10x Genomics **Visium**](https://www.10xgenomics.com/products/spatial-gene-expression), [10x Genomics **Visium HD**](https://www.10xgenomics.com/products/visium-hd-spatial-gene-expression), and [10x Genomics **Xenium**](https://www.10xgenomics.com/platforms/xenium). 

Thank you for your interest in our work!

## Study design

![](code/Fig1_final.png)

**Postmortem human brain tissue collection and spatial transcriptomics data generation.** (A) The human amygdala (pink) in midsagittal (left) and coronal (right) views; the dashed line marks the coronal level and the dashed box the dissected brain block. (B) Dissected coronal block with the amygdala (AMY) outlined (left) and a schematic of amygdala subnuclei at the corresponding A-P level (right). (C) UpSet plot of donors across the three SRT platforms: Visium (n=7), Xenium (n=4), and Visium HD (n=2), with Br9280 profiled on all three. (D) The Visium slide carries four 6.5 mm² capture arrays of ~5,000 barcoded 55 µm spots, capturing the whole transcriptome across multiple cells per spot. (E) Visium workflow for representative donor Br9280: AMY tissue is scored into six 6.5 mm² strips matched to the capture arrays (Visium Targeting), then H&E and SRT data are reassembled to represent the whole tissue (Visium Stitched) and spatial domains are defined. (F) The Xenium slide carries a single large capture array with probes for 366 genes, allowing subcellular resolution. (G) Xenium workflow for Br9280: unscored AMY tissue is placed on the array (Xenium Targeting), and subcellular resolution resolves cell types and spatial domains within the sample. (H) The Visium HD slide carries two 6.5 mm² capture arrays with a continuous lawn of 2 µm barcodes, resolving ~18,000 genes subcellularly. (I) Visium HD workflow for Br9280: AMY tissue is scored into two 6.5 mm² strips (Visium HD Targeting), and detailed spatial domains are determined for selected subnuclei: IA, CeA, and MeA.
[Biorender](https://biorender.com))

## Interactive Websites

All of these interactive websites are powered by the [`Vitessce`](http://vitessce.io/) open source software:

We provide the following interactive websites, organized by dataset with
software labeled by emojis:

- [Visium (n = 7)](TODO_URL)
    + Provides interactive spot-level visualization of full-transcriptome profiles, spatial domains, 
    and predicted cell type locations.
- [Xenium (n = 4)](TODO_URL)
    + Provides interactive cell segmentations with a targeted gene panel, spatial domain, and transferred cell type labels.
- [Visium HD (n = 5)](TODO_URL)
    + Provides interactive cell segmentations with full-transcriptome profliling, spatial domains, and 
        transferred cell type labels.
- [Br9280: Visium, VisiumHD, and Xenium](TODO_URL)
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

Totty, M. S., Bach, S. V., Valentine, M. R., Tippani, M., Maguire, S. E., Del Rosario Alvia, I., Miller, R. A., Kleinman, J. E., Maynard, K. R., Page, S. C., Hyde, T. M., Hicks, S. C., & Martinowich, K. Multiscale spatial transcriptomics resolves the cellular and molecular architecture of the human amygdala. bioRxiv (TODO_YEAR) doi:TODO_DOI.

## Internal

JHPCE location: `/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/`