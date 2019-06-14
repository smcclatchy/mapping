## ---- include=FALSE------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("06-")


## ----set_seed, echo=FALSE------------------------------------------------
set.seed(49237170)
RNGkind("Mersenne-Twister")


## ----scan1perm, eval=FALSE-----------------------------------------------
## operm <- scan1perm(genoprobs = pr, pheno = iron$pheno, Xcovar = Xcovar, n_perm = 1000)


## ----scan1perm_multicore, eval=FALSE-------------------------------------
## operm <- scan1perm(pr, iron$pheno, Xcovar=Xcovar, n_perm=1000, cores=0)


## ----create_summary, echo=FALSE------------------------------------------
# skipped the permutations above and am hard-coding the results here
operm_summary <- structure(c(3.46116054889584, 3.46335159869951), .Dim = 1:2, .Dimnames = list(
                 "0.05", c("liver", "spleen")), class = c("summary.scan1perm",
                 "matrix"), n_perm = structure(c(1000, 1000), .Dim = 1:2, .Dimnames = list(
                 NULL, c("liver", "spleen"))))


## ----summary_scan1perm, eval=FALSE---------------------------------------
## summary(operm)


## ----summary_scan1perm_run, echo=FALSE-----------------------------------
print(operm_summary)


## ----create_summary_B, echo=FALSE----------------------------------------
# skipped the permutations above and am hard-coding the results here
operm_summary_B <- structure(c(2.62593592893683, 3.46116054889584, 2.63661419257339,
                   3.46335159869951), .Dim = c(2L, 2L), .Dimnames = list(c("0.2",
                   "0.05"), c("liver", "spleen")), class = c("summary.scan1perm",
                   "matrix"), n_perm = structure(c(1000, 1000), .Dim = 1:2, .Dimnames = list(
                   NULL, c("liver", "spleen"))))


## ----summary_scan1perm_B, eval=FALSE-------------------------------------
## summary(operm, alpha=c(0.2, 0.05))


## ----summary_scan1perm_B_run, echo=FALSE---------------------------------
print(operm_summary_B)


## ----set_seed_again, echo=FALSE------------------------------------------
set.seed(49237170)


## ----scan1perm_Xsp, eval=FALSE-------------------------------------------
## operm2 <- scan1perm(pr, iron$pheno, Xcovar=Xcovar, n_perm=1000,
##                     perm_Xsp=TRUE, chr_lengths=chr_lengths(map))


## ----create_summary_C, echo=FALSE----------------------------------------
# skipped the permutations above and am hard-coding the results here
operm2_summary <- structure(list(A = structure(c(2.65418133901934, 3.42486248852004,
                  2.54301332030778, 3.22481615651871), .Dim = c(2L, 2L), .Dimnames = list(
                  c("0.2", "0.05"), c("liver", "spleen"))), X = structure(c(3.09791964868745,
                  3.89645628736863, 4.01757985708749, 5.17928309851641), .Dim = c(2L,
                  2L), .Dimnames = list(c("0.2", "0.05"), c("liver", "spleen")))), .Names = c("A",
                  "X"), class = c("summary.scan1perm", "list"), n_perm = structure(c(1000,
                  28243, 1000, 28243), .Dim = c(2L, 2L), .Dimnames = list(c("A",
                  "X"), c("liver", "spleen"))))


## ----summary_scan1perm_C, eval=FALSE-------------------------------------
## summary(operm2, alpha=c(0.2, 0.05))


## ----summary_scan1perm_C_run, echo=FALSE---------------------------------
print(operm2_summary)


## ----set_seed_yet_again, echo=FALSE--------------------------------------
set.seed(49237170)


## ----scan1perm_lmm, eval=FALSE-------------------------------------------
## operm3 <- scan1perm(pr, iron$pheno, kinship_loco, Xcovar=Xcovar, n_perm=1000,
##                     perm_Xsp=TRUE, chr_lengths=chr_lengths(map))


## ----create_summary_D, echo=FALSE----------------------------------------
# skipped the permutations above and am hard-coding the results here
operm3_summary <- structure(list(A = structure(c(2.64158056099834, 3.2863616632467,
                  2.62334670551258, 3.28542642483703), .Dim = c(2L, 2L), .Dimnames = list(
                  c("0.2", "0.05"), c("liver", "spleen"))), X = structure(c(3.13853435409231,
                  3.81699660922065, 4.36511444098698, 5.49658348000439), .Dim = c(2L,
                  2L), .Dimnames = list(c("0.2", "0.05"), c("liver", "spleen")))), .Names = c("A",
                  "X"), class = c("summary.scan1perm", "list"), n_perm = structure(c(1000,
                  28243, 1000, 28243), .Dim = c(2L, 2L), .Dimnames = list(c("A",
                  "X"), c("liver", "spleen"))))


## ----summary_scan1perm_D, eval=FALSE-------------------------------------
## summary(operm3, alpha=c(0.2, 0.05))


## ----summary_scan1perm_D_run, echo=FALSE---------------------------------
print(operm3_summary)


## ----scan1perm_binary, eval=FALSE----------------------------------------
## operm_bin <- scan1perm(pr, bin_pheno, Xcovar=Xcovar, model="binary",
##                        n_perm=1000, perm_Xsp=TRUE, chr_lengths=chr_lengths(map))


## ----scan1perm_binary_summary, echo=FALSE--------------------------------
operm_bin_summary <- structure(list(A = structure(c(2.59831201609688, 3.32949509346178,
                         2.62913892355286, 3.40933644523101), .Dim = c(2L, 2L), .Dimnames = list(
                         c("0.2", "0.05"), c("liver", "spleen"))), X = structure(c(3.16115845941841,
                         3.85862831021874, 3.0614902104213, 3.76657536814501), .Dim = c(2L,
                         2L), .Dimnames = list(c("0.2", "0.05"), c("liver", "spleen")))), .Names = c("A",
                         "X"), class = c("summary.scan1perm", "list"), n_perm = structure(c(1000,
                         28243, 1000, 28243), .Dim = c(2L, 2L), .Dimnames = list(c("A",
                         "X"), c("liver", "spleen"))))


## ----summary_scan1perm_binary, eval=FALSE--------------------------------
## summary(operm_bin, alpha=c(0.2, 0.05))


## ----summary_scan1perm_binary_run, echo=FALSE----------------------------
print(operm_bin_summary)


## ----permute_data, eval=FALSE--------------------------------------------
## shuffled_order <- sample(rownames(iron$pheno))
## pheno_permuted <- iron$pheno
## rownames(pheno_permuted) <- shuffled_order
## xcovar_permuted <- Xcovar
## rownames(xcovar_permuted) <- shuffled_order
## out_permuted <- scan1(genoprobs = pr, pheno = pheno_permuted, Xcovar = xcovar_permuted)
## plot(out_permuted, map)
## head(shuffled_order)

