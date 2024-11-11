# Pleiotropy IVSA and SENS Draft Analysis
# Michael C. Saul
# michael.saul [at] jax.org
# 2023-06-02

# This formal analysis of pleiotropy will determine whether there
# is pleiotropic regulation of the chr7 peaks for both IVSA and
# cocaine locomotor sensitization.

# First, getting information for the pleiotropy analysis
chromosome = "7"
peak_position = 87.650445
pheno_col1 = "RTG_a3"
pheno_col2 = "saline_novelty_cancor"

# Next, loading necessary libraries
library("qtl2")
library("qtl2pleio")
library("ggplot2")
library("dplyr")
library("parallel")

# Next getting cluster config options
n_cores = 29
cluster_n = makeCluster(n_cores, type = "PSOCK")

# Getting fit1_pvl for parallel
fit1_pvl2 = function (indices, start_snp, probs, addcovar, inv_S, S, pheno)
{
  indices <- as.numeric(indices)
  if (is.na(indices[length(indices)]))
    indices <- indices[-length(indices)]
  X_list <- qtl2pleio::prep_X_list(indices = indices, start_snp = start_snp,
                                   probs = probs, covariates = addcovar)
  X <- gemma2::stagger_mats(X_list)
  Bhat <- qtl2pleio::rcpp_calc_Bhat2(X = X, Sigma_inv = inv_S, Y = as.numeric(as.vector(as.matrix(pheno))))
  mymu <- as.vector(X %*% Bhat)
  out <- qtl2pleio::rcpp_log_dmvnorm2(inv_S = inv_S, mu = mymu, x = as.numeric(as.vector(as.matrix(pheno))),
                                      S = S)
  return(as.numeric(out))
}

# Getting pvl scan functions for use with cluster
# replacing mclapply with parLapply
scan_pvl_clean_2 = function (pheno, probs, addcovar,
                             Sigma_inv, Sigma, start_snp,
                             mytab, n_snp, cores = 1)
{
  list_result <- parallel::parLapply(cl = cores, X = as.data.frame(t(mytab)),
                                     fun = qtl2pleio::fit1_pvl, addcovar = addcovar,
                                     probs = probs, inv_S = Sigma_inv,
                                     S = Sigma, start_snp = start_snp, pheno = pheno)
  mytab$loglik <- unlist(list_result)
  marker_id <- dimnames(probs)[[3]][start_snp:(start_snp +
                                                 n_snp - 1)]
  mytab2 <- tibble::as_tibble(apply(FUN = function(x) marker_id[x],
                                    X = mytab[, -ncol(mytab)], MARGIN = 2))
  mytab2$log10lik <- mytab$loglik/log(10)
  return(mytab2)
}

scan_pvl_2 = function (probs, pheno, kinship = NULL,
                       addcovar = NULL, start_snp = 1,
                       n_snp, max_iter = 10000,
                       max_prec = 1/1e+08, cores = 1)
{
  inputs <- qtl2pleio:::process_inputs(probs = probs, pheno = pheno, addcovar = addcovar,
                                       kinship = kinship, max_iter = max_iter, max_prec = max_prec)
  d_size <- ncol(inputs$pheno)
  mytab <- qtl2pleio::prep_mytab(d_size = d_size, n_snp = n_snp)
  out <- scan_pvl_clean_2(mytab = mytab, addcovar = inputs$addcovar,
                          probs = inputs$probs, Sigma_inv = inputs$Sigma_inv, Sigma = inputs$Sigma,
                          start_snp = start_snp, pheno = inputs$pheno, n_snp = n_snp,
                          cores = cores)
  return(out)
}

boot_pvl_2 = function (probs, pheno, addcovar = NULL,
                       kinship = NULL, start_snp = 1,
                       n_snp, pleio_peak_index, nboot = 1,
                       max_iter = 10000, max_prec = 1/1e+08,
                       cores = 1)
{
  inputs <- qtl2pleio:::process_inputs(probs = probs, pheno = pheno, addcovar = addcovar,
                                       kinship = kinship, max_iter = max_iter, max_prec = max_prec)
  X1 <- inputs$probs[, , pleio_peak_index]
  if (!is.null(inputs$addcovar)) {
    Xpre <- cbind(X1, inputs$addcovar)
  }
  else {
    Xpre <- X1
  }
  d_size <- ncol(inputs$pheno)
  Xlist <- vector(length = d_size)
  Xlist <- lapply(Xlist, FUN = function(x) {
    x <- Xpre
    return(x)
  })
  X <- gemma2::stagger_mats(Xlist)
  Bcol <- rcpp_calc_Bhat2(X = X, Sigma_inv = inputs$Sigma_inv,
                          Y = as.vector(as.matrix(inputs$pheno)))
  B <- matrix(data = Bcol, nrow = ncol(Xpre), ncol = d_size,
              byrow = FALSE)
  Ysimlist <- list()
  for (i in 1:nboot) {
    foo <- sim1(X = X, B = B, Sigma = inputs$Sigma)
    Ysim <- matrix(foo, ncol = d_size, byrow = FALSE)
    rownames(Ysim) <- rownames(inputs$pheno)
    colnames(Ysim) <- paste0("t", 1:d_size)
    Ysimlist[[i]] <- Ysim
  }
  mytab <- qtl2pleio::prep_mytab(d_size = d_size, n_snp = n_snp)
  scan_out <- parallel::parLapply(cl = cores, X = Ysimlist, fun = scan_pvl_clean_2,
                                  probs = inputs$probs, addcovar = inputs$addcovar,
                                  Sigma_inv = inputs$Sigma_inv, Sigma = inputs$Sigma,
                                  start_snp = start_snp,
                                  mytab = mytab, n_snp = n_snp)
  lrt <- parallel::parLapply(cl = cores, X = scan_out, fun = function(x) {
    x %>% calc_profile_lods() %>% dplyr::select(profile_lod) %>%
      max()
  })
  return(unlist(lrt))
}


# Next, getting the data
DO_IVSA = readRDS("/projects/chesler-lab/USERS/saulm/RTG/DO_RTG_project_df_for_mapping.RDS")
DO_IVSA$Sex = ifelse(DO_IVSA$Sex == "Female", 0, 1)
DO_SENS = readRDS("/projects/chesler-lab/USERS/saulm/example/DO_novelty_saline_cocaine_cancor.RDS")
DO_pmap = readRDS("/projects/chesler-lab/USERS/saulm/misc/cleanup_chesler_lab/DO_pmap.RDS")

load("/projects/csna/csna_workflow/data/Jackson_Lab_11_batches/gm_DO2816_qc.RData")
load("/projects/csna/csna_workflow/data/Jackson_Lab_11_batches/apr_DO2816.RData")

# Getting phenotypes
pheno = cbind(DO_IVSA,
              DO_SENS[row.names(DO_IVSA),c("saline_novelty_cancor",
                                           "cocaine_novelty_cancor")])

geno.id <- as.character(do.call(rbind.data.frame, strsplit(ind_ids(gm_DO2816_qc), "_"))[,6])
#replace new ids in apr
nodup.geno.id <- data.frame(geno.id = geno.id,
                            nodup.geno.id = NA)
nodup.geno.id <- nodup.geno.id %>% mutate(nodup.geno.id = make.unique(as.character(geno.id)))
nodup.geno.id <- nodup.geno.id$nodup.geno.id
names(nodup.geno.id) <- ind_ids(gm_DO2816_qc)
apr <- replace_ids(apr,nodup.geno.id)

#overlap id
overlap.id <- intersect(pheno$Mouse_ID, nodup.geno.id)

#subset pheno
pheno <- pheno[pheno$Mouse_ID %in% overlap.id,]
dim(pheno)

apr <- apr[overlap.id,]
str(apr)
gm_qc <- replace_ids(gm_DO2816_qc,nodup.geno.id)[overlap.id,]
gm_qc
all.equal(as.character(pheno$Mouse_ID), ind_ids(gm_qc))
all.equal(rownames(apr$`1`), ind_ids(gm_qc))

k = calc_kinship(probs = apr,
                 type = "loco",
                 use_allele_probs = TRUE, cores = 29)

# Getting additive covariant matrix
addcovar = model.matrix(~ Generation, data = pheno)[,-1]

# Getting SNP IDs
pmap_chr = DO_pmap[[chromosome]]
min_chr = min(which(pmap_chr > (peak_position - 10) & pmap_chr < (peak_position + 10)))
n_snps = length(which(pmap_chr > (peak_position - 10) & pmap_chr < (peak_position + 10)))

# Saving input data
# save(list = ls(), file = "~/projects/example/DO_IVSA_SENS_pleio_input.RData")

# Calculating pleiotropy of the two phenotypes
pleio_tib <- scan_pvl_2(probs        = apr[[chromosome]],
                        pheno        = as.matrix(pheno[,c(pheno_col1,pheno_col2)]),
                        kinship      = k[[chromosome]],
                        addcovar     = addcovar,
                        cores        = cluster_n,
                        start_snp    = min_chr,
                        n_snp        = n_snps)
saveRDS(pleio_tib,"~/projects/example/DO_IVSA_SALINE_pleio_pvl.RDS")
pleio_tib_lods = calc_profile_lods(pleio_tib) %>%
  add_pmap(pmap = pmap_chr)
pleio_tib_lods$trait = gsub("^tr1$",pheno_col1,pleio_tib_lods$trait)
pleio_tib_lods$trait = gsub("^tr2$",pheno_col2,pleio_tib_lods$trait)
pleio_tib_lods$trait = gsub("^pleiotropy$","Pleiotropy",pleio_tib_lods$trait)
saveRDS(pleio_tib_lods,"~/projects/example/DO_IVSA_SALINE_pleio_lods.RDS")

# Running bootstraps
pleio_lrt   = calc_lrt_tib(pleio_tib)
pleio_index = find_pleio_peak_tib(pleio_tib, start_snp = min_chr)
# pleio_boot  = boot_pvl_2(probs = apr[[chromosome]],
#                          pheno = as.matrix(pheno[,c(pheno_col1,pheno_col2)]),
#                          pleio_peak_index = pleio_index,
#                          kinship = k[[chromosome]],
#                          nboot = 100,
#                          start_snp = min_chr,
#                          n_snp = n_snps,
#                          cores = cluster_n)
# save(list = c("pleio_lrt","pleio_index","pleio_boot"),
#      file = "~/projects/example/DO_IVSA_SALINE_pleio_boot.RData")

pleio_plot = ggplot(data = pleio_tib_lods, aes(x        = marker_position,
                                               y        = profile_lod,
                                               color    = trait)) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_line(color = "#FFFFFF"),
        legend.position = "bottom",
        panel.spacing = unit(0, "mm")) +
  xlab(paste0("chr",chromosome," Position (Mb)")) +
  ylab("LOD") +
  ggtitle("Pleiotropy Analysis: IVSA vs Saline Sensi")

