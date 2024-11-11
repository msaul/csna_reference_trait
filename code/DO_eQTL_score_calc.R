# LOD score calculation
# Michael C. Saul
# michael.saul@jax.org

# Calculating LOD scores for eQTL with and without proximate markers as additive covariates
# Getting library for qtl2
library("qtl2")

# Importing data
load("./RData/DO_str_2016_eQTL.RData")
cat("Data loaded.\n", sep = "")

# Getting information needed to run through the dataset.
chrs = c(1:19,"X")
eQTL_maRt$cis_marker = rep(NA, times = nrow(eQTL_maRt))
eQTL_maRt$dist_tss = rep(NA, times = nrow(eQTL_maRt))
DO_pheno = DO_str_2016_cross$pheno
DO_covar = DO_str_2016_cross$covar
DO_covar$cis_marker = rep(NA, times = nrow(DO_covar))

# Getting PBS array ID to position in index and number of cores
pbs_array_id = as.numeric(Sys.getenv("PBS_ARRAYID"))
n_array_ids = 100
n_cores = as.numeric(Sys.getenv("PBS_NUM_PPN"))
cat("Running array ID ", pbs_array_id, " of ", n_array_ids, " with ", n_cores, " cores.\n", sep = "")
n_data = ncol(DO_pheno)
n_per = n_data %/% n_array_ids
n_remainder = n_data %% n_array_ids
remainder = n_remainder / n_array_ids
addvar = (pbs_array_id * remainder) %/% 1
addvar_minus_1 = ((pbs_array_id - 1) * remainder) %/% 1
diff_advar = addvar - addvar_minus_1
start_pbs = ((pbs_array_id - 1) * n_per) + addvar_minus_1 + 1
stop_pbs = (pbs_array_id * n_per) + addvar
cat("Calculating LOD scores for columns ", start_pbs, " through ", stop_pbs, ".\n\n", sep = "")
col_total = stop_pbs - start_pbs

# Looping through the columns needed for this specific dataset.
for (i in row.names(eQTL_maRt)[start_pbs:stop_pbs]) {
  cat("Starting gene ", i, ". ", sep = "")
  chr_i = eQTL_maRt[i,"chromosome_name"]
  start_i = eQTL_maRt[i,"start_position"]
  end_i = eQTL_maRt[i,"end_position"]
  
  pheno_col_i = which(colnames(DO_pheno) == i)
  addcovar_i = model.matrix(~ ngen, data = DO_covar)[,-1]
  qtl_cis_i = scan1(genoprobs = DO_str_2016_aprobs,
                    pheno = DO_pheno[,pheno_col_i, drop = FALSE],
                    kinship = DO_str_2016_kinship,
                    addcovar = addcovar_i,
                    cores = n_cores)
  
  if (exists("qtl_cis")) {
    qtl_cis = cbind(qtl_cis,
                    qtl_cis_i)
  } else {
    qtl_cis = qtl_cis_i
  }
  
  if (!(chr_i %in% chrs)) {
    next
  } else {
    # Getting the IDs of all markers within plus or minus 5 Mb of the gene body.
    markers_chr_i = DO_str_2016_cross$pmap[[chr_i]]
    markers_chr_i = markers_chr_i[which(markers_chr_i < ((end_i + 5e6) / 1e6) &
                                          markers_chr_i > ((start_i - 5e6) / 1e6))]
    if (length(markers_chr_i) == 0) {
      next
    } else {
      cis_lods_i = as.data.frame(qtl_cis_i)[names(markers_chr_i),1]
      
      # Using the marker with the highest LOD score within 5 Mb of the gene body as the cis-eQTL
      max_marker_i = names(markers_chr_i)[which(cis_lods_i == max(cis_lods_i))][1]
      # There's an error with a marker that contains only one genotype. Fixing that with an ifelse statement.
      max_marker_i = ifelse(max_marker_i == "UNC24016748", "UNCHS038149", max_marker_i)
      DO_covar$cis_marker = as.data.frame(pull_markers(DO_str_2016_cross,max_marker_i)$geno)[row.names(DO_covar),1]
      DO_covar$cis_marker = as.character(DO_covar$cis_marker)
      eQTL_maRt[i,"cis_marker"] = max_marker_i
      eQTL_maRt[i,"dist_tss"] = markers_chr_i[max_marker_i] - (eQTL_maRt[i,"marker_start"] / 1e6) 
      
      # Running the model with the cis marker as an additive covariate (categorical variable)
      addcovar_ciscovar_i = model.matrix(~ ngen + cis_marker, data = DO_covar)[,-1]
      qtl_nocis_i = scan1(genoprobs = DO_str_2016_aprobs,
                          pheno = DO_pheno[,pheno_col_i, drop = FALSE],
                          kinship = DO_str_2016_kinship,
                          addcovar = addcovar_ciscovar_i,
                          cores = n_cores)
      if (exists("qtl_nocis")) {
        qtl_nocis = cbind(qtl_nocis,
                          qtl_nocis_i)
      } else {
        qtl_nocis = qtl_nocis_i
      }
    }
  }
  DO_covar$cis_marker = rep(NA, times = nrow(DO_covar))
  cat("Finishing gene ", i, ".\n", sep = "")
}

# Writing out eQTL data
cat("Writing out files for LOD scores.\n", sep = "")
saveRDS(qtl_nocis, file = paste("./RDS/DO_lod_scores_nocis_arrayid_",pbs_array_id,".RDS",sep=""),
        compress = "xz")
saveRDS(qtl_cis, file = paste("./RDS/DO_lod_scores_cis_arrayid_",pbs_array_id,".RDS",sep=""),
        compress = "xz")
saveRDS(eQTL_maRt[start_pbs:stop_pbs,], paste("./RDS/DO_lod_scores_eQTL_info_arrayid_",pbs_array_id,".RDS",sep=""),
        compress = "xz")