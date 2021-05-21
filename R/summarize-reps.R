log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")



library(tidyverse)
library(hexbin)
library(rdist)
##CKMRpop::install_spip(Dir = system.file(package = "CKMRpop"))
##remotes::install_github("eriqande/CKMRpop", build_vignettes = TRUE)
library(CKMRpop)
#vignette("species_1_simulation", package = "CKMRpop")


#### these are the values that get changed for different runs ####
cohort_size <- as.integer(snakemake@wildcards$cohort_size)
SampleSize = as.integer(snakemake@wildcards$sample_size)


#### When summarizing all the runs, set it to the number of single-rep inputs ####
NReps <- length(snakemake@input)

SPD <- species_1_life_history

dummy = "2000CohortA"
dummy2 = paste(dummy,"sibs",sep="")

### define the scenario
survival = 0.7
alpha = 3  ## age at maturity
omega = 10 ## maximum age
adultlifespan = omega-alpha+1
phi = 1 ## ratio of Vk to kbar
##femalefecundity = c(0,0,alpha:omega)  ## fecundity proportional to age
femalefecundity = c(0,0,rep(1,adultlifespan))  ## constant fecundity
##femalefecundity = c(rep(0,9),1)  ## sweepstakes RS; only BOFFFs reproduce

samp_frac <- 2*SampleSize/cohort_size  ## twice as large as target for subsampling
SPD$`number-of-years` <- 56  # run the sim forward for 100 years
samp_start_year <- 51
samp_stop_year <- 55

SPD[[4]] = c(0,0,1,1,1,1,1,1,1,1)  ## prob of reproducing at each age
##SPD[[4]] = c(0,0,1,0,1,0,1,0,1,0)
##SPD[[4]] = c(0,0,rep(0.1,8))
SPD[[6]] = femalefecundity

a=as.numeric(Sys.time())
set.seed(a)

##Scenario B
SPD[[1]] = omega
SPD[[2]] = c(1,rep(survival,(SPD[[1]]-1)))
SPD[[3]] = SPD[[2]]

SPD[[5]] = SPD[[4]]

SPD[[7]] = SPD[[6]]
SPD[[8]] = "negbin"
SPD[[9]] = 1/phi
SPD[[10]] = SPD[[9]]
SPD[[11]] = -1
SPD[[12]] = 0.5

L <- leslie_from_spip(SPD, cohort_size)

# then we add those to the spip parameters
SPD$`initial-males` <- floor(L$stable_age_distro_fem)
SPD$`initial-females` <- floor(L$stable_age_distro_male)

# tell spip to use the cohort size
SPD$`cohort-size` <- paste("const", cohort_size, collapse = " ")
SPD$`fixed-cohort-size` <- ""  # define this flag and give it an empty string as an argument

sfspace = paste("0 ")
range = paste(samp_start_year,"-",samp_stop_year,sep="")
SPD$`discard-all` <- 0
SPD$`gtyp-ppn-fem-pre` <- paste(range, "0 ", samp_frac, paste(rep(sfspace, SPD$'max-age' - 2), collapse = ""))
SPD$`gtyp-ppn-male-pre` <- SPD$`gtyp-ppn-fem-pre`

Born = as.integer((samp_start_year:samp_stop_year)-2)
YearsSamp = length(Born)
BigGapHalf = array(data=0, dim = c(YearsSamp,YearsSamp,NReps))
colnames(BigGapHalf) = c("Year1","Year2","Year3","Year4","Year5")
BigGapFull = BigGapHalf
BigCohorts = matrix(NA,NReps,YearsSamp)
colnames(BigCohorts) = Born
BigNb = matrix(NA,NReps,YearsSamp)
BigNbdads = BigNb
BigNbmoms = BigNb
BigSibsSame = array(data = 0, dim = c(YearsSamp,4,NReps))
dimnames(BigSibsSame)[[2]] = c("cohort","FSP","PHSP","MHSP")
BigSibsSame[,1,] = Born




for(R in 1:NReps) {
  tmp <- read_rds(snakemake@input[[R]])
  BigGapHalf[,,R] = tmp$bsh
  BigGapFull[,,R] = tmp$bsf
}  # end for R

########################


ComparisonsWithin = SampleSize*(SampleSize-1)/2

#write.csv(BigNb,file="C:/users/lukes/dropbox/halflings/simulation/SibsB-Nb.csv")
#write.csv(BigNbmoms,file="C:/users/lukes/dropbox/halflings/simulation/SibsB-Nbmoms.csv")
#write.csv(BigNbdads,file="C:/users/lukes/dropbox/halflings/simulation/SibsB-Nbdads.csv")
#write.csv(BigGapHalf,file="C:/users/lukes/dropbox/halflings/simulation/SibsB-GapHalf.csv")
#write.csv(BigGapFull,file="C:/users/lukes/dropbox/halflings/simulation/SibsB-GapFull.csv")
#write.csv(BigSibsSame,file="C:/users/lukes/dropbox/halflings/simulation/SibsB-SibsSame.csv")

MeanGapsHalf = apply(BigGapHalf,1:2,mean)
MeanGapsFull = apply(BigGapFull,1:2,mean)

WangR = matrix(0,NReps,YearsSamp)
for (R in 1:NReps)  {
  matches = 1:YearsSamp
  for (k in 1:YearsSamp)  {
  matches[k] = 2*BigGapFull[k,k,R] + BigGapHalf[k,k,R] }
  WangR[R,] = matches
  }  ## end for R
WangNb = 4*ComparisonsWithin/WangR

TotWangR = matrix(0,NReps,YearsSamp-1)
colnames(TotWangR) = c("2Years","3Years","4Years","5Years")
for (R in 1:NReps)  {
for (j in 1:(YearsSamp-1))  {
  TotWangR[R,j] = sum(WangR[R,1:(1+j)])
  }  # end for j
  }  # end for R

ComboWangNb = matrix(0,NReps,YearsSamp-1)
for (j in 2:YearsSamp)  {
Top = 4*ComparisonsWithin*j
ComboWangNb[,j-1] = Top/TotWangR[,j-1]  }

WangData = matrix(NA,5,YearsSamp-1)
colnames(WangData) = c("2Years","3Years","4Years","5Years")
rownames(WangData) = c("5%","Median","95%","CV1","CV2")
for (j in 1:(YearsSamp-1))  {
WangData[1:3,j] = quantile(ComboWangNb[,j],probs = c(0.05, 0.5, 0.95))
WangData[4,j] = sd(ComboWangNb[,j])/mean(ComboWangNb[,j])
WangData[5,j] = sd(1/ComboWangNb[,j])/mean(1/ComboWangNb[,j])
}

AcrossHalfs = matrix(0,NReps,(YearsSamp-1))
colnames(AcrossHalfs) = c("2Years","3Years","4Years","5Years")
for (R in 1:NReps)  {
f=BigGapHalf[,,R]
f[!upper.tri(f)] <- 0
  for (j in 1:(YearsSamp-1))  {
    temp = f[1:(j+1),1:(j+1)]
  AcrossHalfs[R,j] = sum(temp)
  }  # end for j
  }  # end for R


gaps = 1:(YearsSamp-1)

ComparisonsBetween = matrix(2*SampleSize^2,YearsSamp-1,YearsSamp-1)
colnames(ComparisonsBetween) = c("Year2","Year3","Year4","Year5")
rownames(ComparisonsBetween) = c("Year1","Year2","Year3","Year4")
ComparisonsBetween[lower.tri(ComparisonsBetween)] <- NA
AdjustedComparisonsBetween = ComparisonsBetween

K = c(1,rep(survival,omega-1))
Lx = cumprod(K)
Nx = Lx*cohort_size
AdultN = sum(Nx[alpha:omega])
effectivesurvival = 1:4
for (j in 1:4)  {
effectivesurvival[j] = sum(Nx[(alpha+j):omega])/AdultN
}  # end for j

for (j in 1:(YearsSamp-1))  {
for (k in j:(YearsSamp-1)) {
  q = 1+k-j
  AdjustedComparisonsBetween[j,k] = ComparisonsBetween[j,k]*effectivesurvival[q]
  }}
TotalAdjustedComparisons = 1:(YearsSamp-1)
for (j in 1:(YearsSamp-1)) {
  temp = AdjustedComparisonsBetween[1:j,1:j]
  TotalAdjustedComparisons[j] = sum(temp,na.rm=T)
  }

X = t(2/AcrossHalfs)*TotalAdjustedComparisons
Nhat = t(X)

NData = matrix(NA,5,YearsSamp-1)
colnames(NData) = c("2Years","3Years","4Years","5Years")
rownames(NData) = c("5%","Median","95%","CV1","CV2")
for (j in 1:(YearsSamp-1))  {
NData[1:3,j] = quantile(Nhat[,j],probs = c(0.05, 0.5, 0.95))
NData[4,j] = sd(Nhat[,j])/mean(Nhat[,j])
NData[5,j] = sd(1/Nhat[,j])/mean(1/Nhat[,j])
}

MeanGapHalf = matrix(NA,YearsSamp,YearsSamp)
for (j in 1:YearsSamp)  {
for (k in 1:YearsSamp)  {
MeanGapHalf[j,k] = mean(BigGapHalf[j,k,])
}}

SplitSibs = matrix(NA,NReps,YearsSamp)
TotAcross = rep(NA,NReps)
for (R in 1:NReps)  {
temp = rep(0,YearsSamp)
A = BigGapHalf[,,R]
##A[lower.tri(A,diag=FALSE)] <- 0
for (j in 1:YearsSamp)  {
for (k in j:YearsSamp)  {
  X = 1 + k-j
  temp[X] = temp[X] + A[j,k]
  }}  ## end for j,k
SplitSibs[R,] = temp
TotAcross[R] = sum(temp[2:YearsSamp])
} ## end for R

BigSplit = cbind(SplitSibs,TotAcross)
MeanGaps = colMeans(BigSplit)
VarGaps = 1:(YearsSamp+1)
for (j in 1:(YearsSamp+1))  {
VarGaps[j] = var(BigSplit[,j])
}  # end for j

GapCompare = rowSums(ComparisonsBetween,na.rm=T)
ESibs = GapCompare*effectivesurvival*2/AdultN
GapData = cbind(MeanGaps,VarGaps,VarGaps/MeanGaps,sqrt(VarGaps)/MeanGaps)
colnames(GapData) = c("MeanSibs","Var","Phi","CV")
rownames(GapData) = c(0:(YearsSamp-1),"TotAcross")

### harmonic means
1/mean(1/WangNb)
1/mean(1/BigNb)
1/mean(1/Nhat)

sum(BigGapFull,na.rm=T)

t(NData)
t(WangData)

GapData
ESibs

cohort_size
SampleSize




# write everything out into a list
ret <- list(
  ### harmonic means
  hm_wang_nb = 1/mean(1/WangNb),
  hm_big_nb = 1/mean(1/BigNb),
  hm_hnat = 1/mean(1/Nhat),

  sum_big_gap_full = sum(BigGapFull,na.rm=T),

  tNData = t(NData),
  tWangData = t(WangData),

  GapCompare = GapCompare,
  ESibs = GapCompare*effectivesurvival*2/AdultN,

  GapData = GapData,

  cohort_size = cohort_size,
  SampleSize = SampleSize
)


write_rds(ret, file = snakemake@output[[1]])
