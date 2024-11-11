library("qtl2")
library("tidyverse")

pheno <- readRDS(file = "/projects/chesler-lab/USERS/saulm/example/DO_novelty_saline_cocaine_cancor.RDS")
# pheno$Sex[pheno$Sex == "Male"] <- 1
# pheno$Sex[pheno$Sex == "Female"] <- 0

# geno data ---------------------------------------------------------------
load("/projects/csna/csna_workflow/data/Jackson_Lab_11_batches/gm_DO2816_qc.RData")
load("/projects/csna/csna_workflow/data/Jackson_Lab_11_batches/apr_DO2816.RData")

#new geno id
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

#subset apr and gm
apr <- apr[overlap.id,]
str(apr)
gm_qc <- replace_ids(gm_DO2816_qc,nodup.geno.id)[overlap.id,]
gm_qc
all.equal(as.character(pheno$Mouse_ID), ind_ids(gm_qc))
all.equal(rownames(apr$`1`), ind_ids(gm_qc))

#kinship
k = calc_kinship(probs = apr, type = "loco", use_allele_probs = TRUE, cores = 4)

#pheno list
pheno.list <- colnames(pheno)[4:5]

# qtlmapping  -------------------------------------------------------------
query_variants <- create_variant_query_func("/projects/csna/csna_workflow/data/cc_variants.sqlite")
query_genes <- create_gene_query_func("/projects/csna/csna_workflow/data/mouse_genes_mgi.sqlite")

sigqtl.chr = "7"
pheno$Generation <- as.factor(pheno$Generation)
addcovar.i = model.matrix(~Generation, data = pheno)[,-1,drop=F]

load("/projects/chesler-lab/USERS/saulm/RTG/DO_SENSI_COCAINE/DO_SENSI_COCAINE.qtl.RData")
m2.sigqtl.peak = find_peaks(m2.qtl.out$cocaine_novelty_cancor, gm_qc$pmap, 3, prob = 0.95)
phe.i1 = pheno[,pheno.list[1]]
phe.i2 = pheno[,pheno.list[2]]
names(phe.i1) = row.names(pheno)
names(phe.i2) = row.names(pheno)


chr7_saline_snps <- scan1snps(genoprobs = apr,
                               map = gm_qc$pmap,
                               pheno = phe.i1,
                               kinship = k[[sigqtl.chr]],
                               addcovar = addcovar.i,
                               query_func=query_variants,
                               chr=sigqtl.chr,
                               start=m2.sigqtl.peak[which(m2.sigqtl.peak$chr == 7),"ci_lo"],
                               end=m2.sigqtl.peak[which(m2.sigqtl.peak$chr == 7),"ci_hi"],
                               keep_all_snps=TRUE)

chr7_cocaine_snps <- scan1snps(genoprobs = apr,
                               map = gm_qc$pmap,
                               pheno = phe.i2,
                               kinship = k[[sigqtl.chr]],
                               addcovar = addcovar.i,
                               query_func=query_variants,
                               chr=sigqtl.chr,
                               start=m2.sigqtl.peak[which(m2.sigqtl.peak$chr == 7),"ci_lo"],
                               end=m2.sigqtl.peak[which(m2.sigqtl.peak$chr == 7),"ci_hi"],
                               keep_all_snps=TRUE)

save(list = c("chr7_saline_snps","chr7_cocaine_snps"), file = "./data/SENSI_chr7_snps.RData")
