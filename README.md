# master_thesis
Master thesis in Systems Biology on the topic: Prediction of cell type-specific hypothalamic neuronal response patterns to fasting, feeding, and refeeding through 3D-multi-omics data integration.

A. **Spatial cell type prediction**
For spatial cell type prediction using ISH data the feature selection procedure is as follows: (1) ISH filter applied, (2) MRx3 applied, (3) MRMR applied

Then the feature selection can be inputed into either MISS or Tangram algorithm.

For spacial cell type prediction using spatial transcriptomics data, cell2location or tangram script can be followed.

B. **scRNA-seq analysis**
To conduct scRNA dat analysis scRNA-seq script can be followed togtehr with the kNN script for label transfer using Hypomap, scRNA-seq reference.

C. **cFos LSFM analysis**
To conduct the cFos LSFM analysis: (1) preprocessing applied, (2) cfos LSFM downscaling applied, (3) cfos LSFM integration with spatila predictions applied.
