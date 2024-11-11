# eQTL analysis script
# Michael C. Saul
# michael.saul@jax.org

load("~/Box/projects/DO_RTG/data/DO_str_2016_eQTL.RData")
DO_pheno = DO_str_2016_cross$pheno
DO_covar = DO_str_2016_cross$covar
addcovar_mat = model.matrix(~ ngen + sex, data = DO_covar)[,-1]
addcovar_i = addcovar_mat
library("qtl2")

plot_eQTL_str_2016 = function(gene_i, n_cores = 2, label = gene_i, color = "darkblue") {
  pheno_col_i = which(colnames(DO_pheno) == gene_i)
  qtl_cis_i = scan1(genoprobs = DO_str_2016_aprobs,
                    pheno = DO_pheno[,pheno_col_i, drop = FALSE],
                    kinship = DO_str_2016_kinship,
                    addcovar = addcovar_i,
                    cores = n_cores)
  return(plot_scan1(qtl_cis_i, map = DO_str_2016_cross$pmap, col = color, main = label))
}


okabe_ito = c("#F0E442", "#999999", "#E69F00", "#0072B2",
              "#56B4E9", "#009E73","#D55E00", "#CC79A7")
broman_do = c("#FFDC00", "#888888", "#F08080", "#0064C9",
              "#7FDBFF", "#2ECC40", "#FF4136", "#B10DC9")

plot_BLUP_str_2016 = function(gene_i,
                              chr_i,
                              label = gene_i,
                              n_cores = 2,
                              legendpos = "bottomleft",
                              founders_labels = c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"),
                              colors = okabe_ito) {
  pheno_col_i = which(colnames(DO_pheno) == gene_i)
  blup_i = scan1blup(DO_str_2016_aprobs[,chr_i],
                     pheno = DO_pheno[,pheno_col_i, drop = FALSE],
                     kinship = DO_str_2016_kinship[[chr_i]],
                     addcovar = addcovar_i,
                     cores = n_cores)
  blup_i = blup_i[,LETTERS[1:8]]
  colnames(blup_i) = founders_labels
  return(plot_coef(blup_i, map = DO_str_2016_cross$pmap, legend = legendpos,
                   main = label, bgcolor = "white", col = colors))
}



