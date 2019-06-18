## ---- include=FALSE------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("12-")


## ----read_DOex_data------------------------------------------------------
DOex <- read_cross2(file = "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex/DOex.zip")


## ----DOex_calc_genoprob--------------------------------------------------
pr <- calc_genoprob(DOex, error_prob=0.002)
apr <- genoprob_to_alleleprob(pr)


## ----DOex_calc_kinship---------------------------------------------------
k <- calc_kinship(apr, "loco")


## ----DOex_create_sex_covar-----------------------------------------------
sex <- (DOex$covar$Sex == "male")*1
names(sex) <- rownames(DOex$covar)


## ----DOex_set_names------------------------------------------------------
sex <- setNames( (DOex$covar$Sex == "male")*1, rownames(DOex$covar) )


## ----DOex_scan1_pg-------------------------------------------------------
out <- scan1(apr, DOex$pheno, k, sex)


## ----plot_DOex_scan------------------------------------------------------
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot(out, DOex$gmap)


## ----DOex_effects_c2-----------------------------------------------------
coef_c2 <- scan1coef(apr[,"2"], DOex$pheno, k[["2"]], sex)


## ----plot_DOex_effects---------------------------------------------------
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_coefCC(coef_c2, DOex$gmap["2"], bgcolor="gray95", legend="bottomleft")


## ----plot_DOex_lod_curve-------------------------------------------------
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_coefCC(coef_c2, DOex$gmap["2"], scan1_output=out, bgcolor="gray95", legend="bottomleft")


## ----query_variants------------------------------------------------------
query_variants <- create_variant_query_func("../data/cc_variants.sqlite")


## ----grab_chr2_variants--------------------------------------------------
variants_2_97.5 <- query_variants(2, 97, 98)


## ----query_MGI-----------------------------------------------------------
query_genes <- create_gene_query_func("../data/mouse_genes_mgi.sqlite")


## ----query_chr2_genes----------------------------------------------------
genes_2_97.5 <- query_genes(2, 97, 98)


## ----DOex_find_peak_in_Mbp-----------------------------------------------
peak_Mbp <- max(out, DOex$pmap)$pos


## ----DOex_find_peak_on_chr2----------------------------------------------
variants <- query_variants(2, peak_Mbp - 1, peak_Mbp + 1)


## ----make_lookup_table---------------------------------------------------
out_snps <- scan1snps(pr, DOex$pmap, DOex$pheno, k[["2"]], sex, query_func=query_variants,
                      chr=2, start=peak_Mbp-1, end=peak_Mbp+1, keep_all_snps=TRUE)


## ----plot_snp_asso-------------------------------------------------------
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_snpasso(out_snps$lod, out_snps$snpinfo)


## ----plot_snp_asso_wplot-------------------------------------------------
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot(out_snps$lod, out_snps$snpinfo)


## ----id_and_plot_genes---------------------------------------------------
genes <- query_genes(2, peak_Mbp - 1, peak_Mbp + 1)
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot(out_snps$lod, out_snps$snpinfo, drop_hilit=1.5, genes=genes)


## ----id_top_snps---------------------------------------------------------
top <- top_snps(out_snps$lod, out_snps$snpinfo)
print(top[,c(1, 8:15, 20)], row.names=FALSE)


## ----gwas_scan-----------------------------------------------------------
out_gwas <- scan1snps(pr, DOex$pmap, DOex$pheno, k, sex, query_func=query_variants, cores=0)


## ----plot_gwas_scan------------------------------------------------------
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot(out_gwas$lod, out_gwas$snpinfo, altcol="green4", gap=0)

