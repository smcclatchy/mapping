# create and store 1000 permutations for iron data
# load these in for episode 06 Performing a permutation test

library(qtl2)

# read in cross data, insert pseudomarkers into map, calculate genotype
# probabilities, get X chromosome covariates, convert phenotype to binary
iron <- read_cross2( system.file("extdata", "iron.zip", package="qtl2") )
map <- insert_pseudomarkers(map=iron$gmap, step=1)
pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002)
Xcovar <- get_x_covar(iron)
bin_pheno <- apply(iron$pheno, 2, function(a) as.numeric(a > median(a)))
rownames(bin_pheno) <- rownames(iron$pheno)

# 1000 permutations
operm <- scan1perm(genoprobs = pr, 
                   pheno = iron$pheno, 
                   Xcovar = Xcovar, 
                   n_perm = 1000)
save(operm, file = "data/operm.Rdata")

# X chromosome permutations
operm2 <- scan1perm(genoprobs = pr, 
                    iron$pheno, 
                    Xcovar=Xcovar, 
                    n_perm=1000,
                    perm_Xsp=TRUE, 
                    chr_lengths=chr_lengths(map))
save(operm2, file = "data/operm2.Rdata")

# permutations for binary phenotype
operm_bin <- scan1perm(genoprobs = pr, 
                       bin_pheno, 
                       Xcovar=Xcovar, 
                       model="binary",
                       n_perm=1000,
                       perm_Xsp=TRUE, 
                       chr_lengths=chr_lengths(map))
save(operm_bin, file = "data/operm_bin.Rdata")