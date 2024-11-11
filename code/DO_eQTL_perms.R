# LOD permutation
# Michael C. Saul
# michael.saul@jax.org

# Permuting LOD scores for with proximate markers as additive covariates
# Getting library for qtl2
library("qtl2")
library("parallel")

# Importing data
load("./RData/DO_str_2016_eQTL.RData")
cat("Data loaded.\n", sep = "")

# Getting information needed to run through the dataset.
DO_pheno = DO_str_2016_cross$pheno
DO_covar = DO_str_2016_cross$covar

set.seed(12345)
sample_pheno = sample(1:ncol(DO_pheno),replace=FALSE)

# Getting PBS array ID to position in index and number of cores
pbs_array_id = as.numeric(Sys.getenv("PBS_ARRAYID"))
n_cores = as.numeric(Sys.getenv("PBS_NUM_PPN"))

gene_i = colnames(DO_pheno)[sample_pheno[pbs_array_id]]
cat("Starting gene ", gene_i, ".\n", sep = "")
chr_i = eQTL_maRt[gene_i,"chromosome_name"]
start_i = eQTL_maRt[gene_i,"start_position"]
end_i = eQTL_maRt[gene_i,"end_position"]

pheno_col_i = which(colnames(DO_pheno) == gene_i)
addcovar_i = model.matrix(~ ngen + sex, data = DO_covar)[,-1]
qtl_cis_i = scan1(genoprobs = DO_str_2016_aprobs,
                  pheno = DO_pheno[,pheno_col_i, drop = FALSE],
                  kinship = DO_str_2016_kinship,
                  addcovar = addcovar_i,
                  cores = n_cores)
qtl_perm_i = scan1perm(genoprobs = DO_str_2016_aprobs,
                       pheno = DO_pheno[,pheno_col_i, drop = FALSE],
                       kinship = DO_str_2016_kinship,
                       addcovar = addcovar_i,
                       cores = n_cores,
                       n_perm = 1000)
qtl_perm = as.data.frame(qtl_perm_i[,1], stringsAsFactors = FALSE)
colnames(qtl_perm)[1] = gene_i

cat("Finishing gene ", gene_i,".\n", sep = "")


# Writing out eQTL data
cat("Writing out file for permutation scores.\n", sep = "")
saveRDS(qtl_perm, paste("./RDS/DO_perms_eQTL_info_arrayid_",pbs_array_id,".RDS",sep=""),
        compress = "xz")