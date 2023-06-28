#IMPORTS
library('reticulate')
library('abind')
library('freesurferformats')
library('tiff')
library('moments')
np = import("numpy")

#PROCESS CFOS DATA
nic_res = np$load('nicnew_fedref_rychannelpad_reshape200.npy')
nic_melt <- reshape2::melt(nic_res)
annot = get_ccf_grid_annotation()
annot <- reshape2::melt(annot)
colnames(annot)[4] = c('annot')
nic_melt_annot = left_join(nic_melt, annot, by = c('Var1', 'Var2', 'Var3'))

mba_ontology = cocoframer::get_mba_ontology()
mba_ontology_flatten = cocoframer::flatten_mba_ontology(mba_ontology)
mba_ontology_flatten = rename(mba_ontology_flatten, annot = id)
nic_melt_annot = left_join(nic_melt_annot, mba_ontology_flatten, by = c('annot'))

mba_ontology_flatten = cocoframer::flatten_mba_ontology(mba_ontology)
h_children = filter_mba_ontology_children(mba_ontology_flatten, parent_acronym = "HY", 
                                          include_parent = TRUE)
h_children = h_children[-c(28),]
list_children = h_children$name
nic_annot_hyp = nic_melt_annot[nic_melt_annot$name %in% list_children,]

#CELL TYPE EXPRESSION LEVEL DISTRIBUTIONS FROM PPREDICTION
annot_mrmr = np$load("/beegfs/home/pmatyskova/project/ish_annot_hypnoSFO.npy", allow_pickle=TRUE)
annot_mrmr = annot_mrmr[[1]]
annot_mrmr = as.data.frame(annot_mrmr)

open_annot = function(name_d, e_annot){
  add_col_names = function(matrix) {
    colnames(matrix) = cell_types
    matrix = as.data.frame(matrix)
    return(matrix)
  }
  annotate_d = function(matrix, annot) {
    annot$num = rownames(annot)
    matrix$num = rownames(matrix)
    matrix_ann = left_join(annot, matrix, by = 'num')
    return(matrix_ann)
  }
  d_matrix <- np$load(name_d)
  d_matrix <- add_col_names(d_matrix)
  d_matrix <- annotate_d(d_matrix, e_annot)
}

dd = data.table::fread('/beegfs/scratch/bruening_scratch/lsteuernagel/data/petra_data/hypoMap_avg_signatures/hypomap_avg_signatures_C185_rna.txt',data.table = FALSE)
cell_types = colnames(dd[,-c(1)])
d_matrix = open_annot("d_miss_hm185cor_minres_mrx31300_filled0.npy", annot_mrmr)

d_matrix_plot = d_matrix[,-c(4:15)]
d_matrix_plot = d_matrix_plot[, colSums(d_matrix_plot != 0) > 0]
cell_type_list = colnames(d_matrix_plot)[4:112]

#cell type first threshold
kurt_list = c()
skew_list = c()
quant_list = c()
itter = 1
for (i in cell_type_list){
  kurt = kurtosis(d_matrix_plot[,itter+3])
  kurt_list[itter] = kurt
  skew = skewness(d_matrix_plot[,itter+3])
  skew_list[itter] = skew
  quant = quantile(d_matrix_plot[,itter+3], probs = c(.98))
  quant_list[itter] = quant
  itter = itter + 1
}
cellt_kurtskew = data.frame(cell_type  = cell_type_list, kurtosis = kurt_list, 
                            skewness = skew_list, quantile98 = quant_list)

kurt_filt = filter(cellt_kurtskew, kurtosis>25)
d_matrix_kurtfilt = d_matrix_plot[,colnames(d_matrix_plot) %in% kurt_filt$cell_type]
d_matrix_kurtfilt = cbind(d_matrix_plot[1:3], d_matrix_kurtfilt)

for (i in 4:ncol(d_matrix_kurtfilt)){
  d_matrix_kurtfilt[i][d_matrix_kurtfilt[i] < kurt_filt$quantile98[i-3]] = 0
}


#COMBINE SPACIAL CELL TYPE PREDICTION WITH CFOS
predict_celltype = function(cellt_predict, cfos, weighted_ct){
  if (weighted_ct == TRUE){
    cellt_predict[4:ncol(cellt_predict)] = sweep(cellt_predict[4:ncol(cellt_predict)], 2, 
                                                 colSums(cellt_predict[4:ncol(cellt_predict)]), `/`) #celltype normilized
  } else {
    cellt_predict[4:ncol(cellt_predict)][cellt_predict[4:ncol(cellt_predict)] > 0] = 1
  }
  
  cellt_predict$num = c(1:2031)
  cfos$num = c(1:2031)
  predic_list = c()
  celln_list = c()
  itter = 1
  for (i in 4:(ncol(cellt_predict)-1)){
    d_filt = filter(cellt_predict, cellt_predict[i] > 0)[c(1:3,i,ncol(cellt_predict))] #vars, cell type column, num merge column
    d_cfos = left_join(d_filt, cfos[c(4,9,ncol(cfos))], by = 'num')
    d_cfos$pred = d_cfos[4] * d_cfos$value #taking cFos values at the place with predicted cell types
    predic_list[itter] = sum(as.numeric(unlist(d_cfos$pred))) #taking a mean
    celln_list[itter] = colnames(d_cfos)[4]
    itter = itter + 1
  }
  
  cfos_cellt = data.frame(cell_type  = celln_list, cfos_mean = predic_list)
  return(cfos_cellt)  
}

mrx3_anoj = predict_celltype(d_matrix_kurtfilt, nic_annot_hyp, FALSE)
mrx3_anoj_norm = predict_celltype(d_matrix_kurtfilt, nic_annot_hyp, TRUE)
data5=mrx3_anoj[grepl("GLU|GABA", mrx3_anoj$cell_type),]
write.csv(data2, "nick_celltp_rychannel.csv", row.names = FALSE)

#compare with scRNA-seq analysis
singlecell_pred = data.table::fread('/beegfs/scratch/bruening_scratch/lsteuernagel/data/petra_data/nucseq_fasting_activated.txt',data.table = FALSE)

datafafei <- read_csv("/beegfs/home/pmatyskova/singularity/results/scrnaseq/list_fer_down_n.csv")
datafafei = filter(datafafei, min_pval < 0.1)
datafafei_joint <- merge(datafafei, data5, by = "cell_type")
