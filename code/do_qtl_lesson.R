# DO QTL mapping
# Sue McClatchy
# 2020-12-18

library(qtl2)
load('data/qtl2_demo.Rdata')
str(pheno)

# Challenge 1 Look at the Data
# 
# Determine the dimensions of pheno.
# 1). How many rows and columns does it have?
# 2). What are the names of the variables it contains?
# 3). What does the distribution of the prop.bm.MN.RET column look like? Is it normally distributed?

str(pheno)
hist(pheno$prop.bm.MN.RET)
pheno$log.MN.RET <- log(pheno$prop.bm.MN.RET)
hist(pheno$log.MN.RET)
dim(probs[[1]]) # samples x founders x markers

# plot genoprobs for 1st mouse on chr 1
plot_genoprob(probs, map, ind = 1, chr = 1)

# calculate kinship matrix
K <- calc_kinship(probs, type = "loco")

# heatmap
n_samples <- 50
heatmap(K[[1]][1:n_samples, 1:n_samples])

# create an additive covariate containing study cohort only; benzene 
# concentration would be another important covariate to include if we were 
# using all mice instead of a subset at one benzene concentration level
# also sex would be an important covariate except in this study, where all mice 
# were male
addcovar <- model.matrix(~Study, data = pheno[, -1])

# select only the 100 ppm exposure group
c100 <- which(pheno$Conc == 100)

# map untransformed phenotype for the 100ppm exposure group
index <- which(colnames(pheno) == "prop.bm.MN.RET")
qtl <- scan1(genoprobs = probs,
             pheno = pheno[c100, index, drop = FALSE],
             kinship = K,
             addcovar = addcovar)

# plot the genome scan
plot_scan1(qtl, map, 
           main = "Proportion of Micro-nucleated Bone Marrow Reticulocytes")

# use log-transformed reticulocytes
index <- which(colnames(pheno)=="log.MN.RET")

# redefine genome scan with log-transformed phenotype at 100ppm exposure level
qtl <- scan1(genoprobs = probs, pheno = pheno[c100, index, drop = FALSE], 
             kinship = K, addcovar = addcovar)

# run 100 permutations
perms <- scan1perm(genoprobs = probs, 
          pheno = pheno[c100,index, drop = FALSE],
          addcovar = addcovar,
          n_perm = 100)
