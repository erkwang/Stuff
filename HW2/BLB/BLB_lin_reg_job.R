#STA 250 Homework 2
#Problem 1
#Bag of Little Bootstraps

# Load packages:
library(BH)
library(bigmemory.sri)
library(bigmemory)
library(biganalytics)

# mini or full?
mini = FALSE
if (mini){
  rootfilename <- "blb_lin_reg_mini"
} else {
  rootfilename <- "blb_lin_reg_data"
}

# I/O specifications:
datapath <- "/home/pdbaines/data/"
outpath <- "~/STA250/Stuff/HW2/BLB/output/"

# Find r and s indices:
s = 5
r = 50

#============================== Setup for running on Gauss... ==============================#

args <- commandArgs(TRUE)

cat("Command-line arguments:\n")
print(args)

####
# sim_start ==> Lowest possible dataset number
###

###################
sim_start <- 1000
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
} else {
  # SLURM can use either 0- or 1-indexing...
  # Lets use 1-indexing here...
  #figure out which distinct subset it should use
  s_index = as.numeric(args[1]) - floor(as.numeric(args[1])/s)*s
  if (!s_index) s_index = 5
  r_index = ceiling(as.numeric(args[1])/s)
  sim_num <- sim_start + as.numeric(args[1])
  sim_seed <- (762*(s_index-1) + 121231)
  set.seed(sim_seed)
}

cat(paste("\nAnalyzing dataset number ",sim_num,"...\n\n",sep=""))
#==============================================================================================#

# Attach big.matrix:
full.data = attach.big.matrix(paste(datapath, rootfilename, ".desc", sep = ""))
# Remaining BLB specs:
gam = 0.7
b = floor(nrow(full.data)^gam)
# Extract the subset:
sample.index = sample(1:nrow(full.data), b, replace = FALSE)
y = full.data[sample.index,ncol(full.data)]
x = full.data[sample.index,-ncol(full.data)]
# Reset simulation seed:
set.seed(123*(sim_num - 1) + 121231)
# Bootstrap dataset:
boot.weights = as.numeric(rmultinom(1 , nrow(full.data), prob = rep(1, b)))
# Fit lm:
boot.fit = lm(y~-1+x, weights=boot.weights)
# Output file:
outfile = paste0(outpath,"coef_",sprintf("%02d",s_index),"_",sprintf("%02d",r_index),".txt", sep = "")
# Save estimates to file:
est = as.data.frame(coef(boot.fit))
write.table(est, file = outfile, row.names=FALSE, col.names = FALSE)




