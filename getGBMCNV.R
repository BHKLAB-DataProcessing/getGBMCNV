input_dir <- "/pfs/downloadGBMData/"
cell_dir <- "/pfs/getGBMCellData/"
out_dir <- "/pfs/out/"

# input_dir <- "~/Documents/pfs/downloadGBMData/"
# cell_dir <- "~/Documents/pfs/getGBMCellData/"
# out_dir <- "~/Documents/pfs/getGBMCNV/" 
functions <- "https://github.com/BHKLAB-Pachyderm/getGBMCellData/raw/main/functions.R"

source(functions)
load(paste0(input_dir, "Ensembl.v99.annotation.RData"))
cell <- readRDS(paste0(cell_dir, "cell.rds"))

assay_cnv<- read.delim(paste(input_dir,"HGCC_DNA_copy_number_gene_level.txt", sep=""), header=T, sep="\t", stringsAsFactors = FALSE)
rownames(assay_cnv)<-assay_cnv$Row
assay_cnv<-assay_cnv[,-1]

feat_cnv<-fdata_builder(annotation=features_gene, assay=assay_cnv,ID_column="gene_name")
phen_cnv<-ph_data_builder(annotation=cell,assay=assay_cnv)
phen_cnv$Replicate[phen_cnv$cellid=="U10000" |phen_cnv$cellid=="U10001" ]<-NA

#Creating ExpressionSet 
assay_cnv<-assay_cnv[,rownames(phen_cnv)]#rearranging the colnames so it is similar to pheno data
cnv_eSet<- ExpressionSet(assayData = as.matrix(assay_cnv), phenoData = AnnotatedDataFrame(phen_cnv), featureData = AnnotatedDataFrame(feat_cnv)) 
print("CNV: done")

cnv_SE <- eSetToSE(cnv_eSet,annot_name="cnv")
saveRDS(cnv_SE, paste0(out_dir, "cnv_SE.rds"))
saveRDS(phen_cnv, paste0(out_dir, "phen_cnv.rds"))