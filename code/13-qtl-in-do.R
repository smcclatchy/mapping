## ---- include=FALSE------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("13-")


## ----message=FALSE, echo=FALSE-------------------------------------------
library(qtl2)


## ----load_data, message=FALSE, results='hide'----------------------------
load("../data/qtl2_demo.Rdata")
sessionInfo()


## ----log_transform-------------------------------------------------------
pheno$log.MN.RET = log(pheno$prop.bm.MN.RET)


## ----hist_log_transform--------------------------------------------------
hist(pheno$log.MN.RET)


## ----dim_probs-----------------------------------------------------------
dim(probs[[1]])


## ----geno_plot, fig.width=8, fig.height=6--------------------------------
image(1:500, 1:ncol(probs[[1]]), t(probs[[1]][1,8:1,1:500]), breaks = 0:100/100,
      col = grey(99:0/100), axes = F, xlab = "Markers", ylab = "Founders",
      main = "Founder Allele Contributions for Sample 1")
abline(h = 0:8 + 0.5, col = "grey70")
usr = par("usr")
rect(usr[1], usr[3], usr[2], usr[4])
axis(side = 1, at = 0:5 * 100, labels = 0:5 * 100)
axis(side = 2, at = 1:8, labels = LETTERS[8:1], las = 1, tick = F)


## ----kinship, message=FALSE, results='hide'------------------------------
K = calc_kinship(probs = probs, type = "loco", use_allele_probs = TRUE)


## ----kinship_probs, fig.width=8, fig.height=8----------------------------
image(1:nrow(K[[1]]), 1:ncol(K[[1]]), K[[1]][,ncol(K[[1]]):1], xlab = "Samples", 
      ylab = "Samples", yaxt = "n", main = "Kinship between samples", 
      breaks = 0:100/100, col = heat.colors(length(0:100) - 1))


## ----covariates----------------------------------------------------------
addcovar = model.matrix(~Study, data = pheno)[,-1]


## ----select_100ppm-------------------------------------------------------
c100 = which(pheno$Conc == 100)


## ----QTL, warning=FALSE, error=FALSE-------------------------------------
index = which(colnames(pheno) == "prop.bm.MN.RET")
qtl = scan1(genoprobs = probs, pheno = pheno[c100,index, drop = FALSE], kinship = K, addcovar = addcovar)


## ----qtl_plot, fig.width=8, fig.height=6, warning=FALSE------------------
plot_scan1(x = qtl, map = map, main = "Proportion of Micro-nucleated Bone Marrow Reticulocytes")


## ----set_pheno_index_to_log----------------------------------------------
index = which(colnames(pheno) == "log.MN.RET")


## ----perms, message=FALSE, results='hide', warning=FALSE-----------------
perms = scan1perm(genoprobs = probs, pheno = pheno[c100,index, drop = FALSE], addcovar = addcovar, n_perm = 100)


## ----qtl_plot_thr, fig.width=8, fig.height=6, warning=FALSE--------------
plot(x = qtl, map = map,  main = "Proportion of Micro-nucleated Bone Marrow Reticulocytes")
thr = summary(perms)
abline(h = thr, col = "red", lwd = 2)


## ----find_peaks----------------------------------------------------------
find_peaks(scan1_output = qtl, map = map, threshold = thr)


## ----interval------------------------------------------------------------
find_peaks(scan1_output = qtl, map = map, threshold = thr, prob = 0.95)


## ----coef----------------------------------------------------------------
chr = 10
coef10 = scan1blup(genoprobs = probs[,chr], pheno = pheno[c100,index, drop = FALSE], kinship = K[[chr]], addcovar = addcovar)


## ----coef_plot, fig.width=8, fig.height=6--------------------------------
plot_coefCC(x = coef10, map = map, scan1_output = qtl, main = "Proportion of Micro-nucleated Bone Marrow Reticulocytes")


## ----assoc_map-----------------------------------------------------------
chr = 10
start = 30
end = 36
query_func = create_variant_query_func("../data/cc_variants.sqlite")
assoc = scan1snps(genoprobs = probs[,chr], map = map, pheno = pheno[c100,index,drop = FALSE], kinship = K, addcovar = addcovar, query_func = query_func, chr = chr, start = start, end = end, keep_all_snps = TRUE)


## ----assoc_fig, fig.width=9, fig.height=6--------------------------------
plot_snpasso(scan1output = assoc$lod, snpinfo = assoc$snpinfo, main = "Proportion of Micro-nucleated Bone Marrow Reticulocytes")


## ----get_genes-----------------------------------------------------------
query_genes = create_gene_query_func(dbfile = "../data/mouse_genes.sqlite", filter = "source='MGI'")
genes = query_genes(chr, start, end)
head(genes)


## ----plot_assoc2,fig.width=12,fig.height=12------------------------------
plot_snpasso(assoc$lod, assoc$snpinfo, main = "Proportion of Micro-nucleated Bone Marrow Reticulocytes", genes = genes)

