library(tidyverse)
library(hexbin)
library(rdist)
##CKMRpop::install_spip(Dir = system.file(package = "CKMRpop"))
##remotes::install_github("eriqande/CKMRpop", build_vignettes = TRUE)
library(CKMRpop)
#vignette("species_1_simulation", package = "CKMRpop")

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
cohort_size <- 2000
SampleSize = 100 ## fixed number to subsample from each cohort
samp_frac <- 2*SampleSize/cohort_size  ## twice as large as target for subsampling
SPD$`number-of-years` <- 56  # run the sim forward for 100 years
samp_start_year <- 51
samp_stop_year <- 55
NReps = 5
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

#########
GetSibs <- function(Pedigree2)  {
# remove duplicate rows in the pedigree
Pedigree2 = Pedigree2[!duplicated(Pedigree2), ]

# convert the parent IDs to unique integers (faster to compare)
unique_moms = sort(unique(Pedigree2$Mom))
unique_dads = sort(unique(Pedigree2$Dad))
mom_range = 1:length(unique_moms)
dad_range = 1:length(unique_dads)
names(mom_range) = unique_moms
names(dad_range) = unique_dads
Pedigree2$Mom_id = mom_range[Pedigree2$Mom]
Pedigree2$Dad_id = dad_range[Pedigree2$Dad]

# make n x n matrix (n=number of offspring) values are zero if the pair of offspring shares a parent, positive otherwise
mom_matrix = pdist(Pedigree2$Mom_id)
dad_matrix = pdist(Pedigree2$Dad_id)
# extract the inds sharing parents, don't double count, and sort
mom_matches = which(mom_matrix==0,arr.ind = T)
mom_matches = mom_matches[mom_matches[,1] < mom_matches[,2], ]
dad_matches = which(dad_matrix==0,arr.ind = T)
dad_matches = dad_matches[dad_matches[,1] < dad_matches[,2], ]
mom_matches = mom_matches[order(mom_matches[,1], mom_matches[,2]),]
dad_matches = dad_matches[order(dad_matches[,1], dad_matches[,2]),]

# convert to data.frame - i1, i2 gives the individual - row number in the original pedigree file
mom_df = data.frame(mom_matches)
names(mom_df) = c('i1', 'i2')
mom_df$parent = 'mom'
dad_df = data.frame(dad_matches)
names(dad_df) = c('i1', 'i2')
dad_df$parent = 'dad'

# merge the dfs of the mom and dad matches to find full sibs
sibs = merge(mom_df, dad_df, by = c('i1', 'i2'), all=T)
# arrays below are boolean indexes into the sibs df, telling us how the pair shares parents
share_both = !(is.na(sibs$parent.x) | is.na(sibs$parent.y))
share_mom = is.na(sibs$parent.x)
share_dad = is.na(sibs$parent.y)
# we can count the number of pairs in each category
##sum(share_both)
##sum(share_mom)
##sum(share_dad)

# construct the output file
ms = cbind(Pedigree2[sibs[share_mom,]$i1,][, c('sampleID', 'YearBorn')],
  Pedigree2[sibs[share_mom,]$i2,][, c('sampleID', 'YearBorn')])
names(ms) = c('ID1', 'Born1', 'ID2', 'Born2')
ms$Parent = 'mom'

ds = cbind(Pedigree2[sibs[share_mom,]$i1,][, c('sampleID', 'YearBorn')],
           Pedigree2[sibs[share_mom,]$i2,][, c('sampleID', 'YearBorn')])
names(ds) = c('ID1', 'Born1', 'ID2', 'Born2')
ds$Parent = 'dad'

if (sum(share_both)>0)  {
  fs = cbind(Pedigree2[sibs[share_both,]$i1,][, c('sampleID', 'YearBorn')],
    Pedigree2[sibs[share_both,]$i2,][, c('sampleID', 'YearBorn')])
  names(fs) = c('ID1', 'Born1', 'ID2', 'Born2')
  fs$Parent = 'both'
  all_sibs = rbind(fs, ms, ds)  }
else {  all_sibs = rbind(ms, ds)  }  

return(all_sibs)  }  # end function

######### start simulation

for (R in 1:NReps)  {

print(paste0("Replicate = ",R))
flush.console()  

BigSibsHalf = matrix(0,YearsSamp,YearsSamp)
BigSibsFull = BigSibsHalf

spip_dir <- run_spip(
  pars = SPD  
)
# now read that in and find relatives within the one-generation pedigree
slurped <- slurp_spip(spip_dir, 1)

# First, get the non-genotype info for each individual all together
non_geno_stuff <- slurped$samples %>%
  mutate(YearSampled = map_int(samp_years_list, 1))  %>% # this gets the sample year out of the samp_years_list
  select(ID, born_year, YearSampled, sex) %>%  # pick out column in the order desired
  left_join(slurped$pedigree %>% select(kid, ma, pa), by = c("ID" = "kid"))  %>%  # add mom and dad on there
  rename(
    sampleID = ID,
    YearBorn = born_year,
    Sex = sex,
    Mom = ma,
    Dad = pa
  )   # change the column names to what Robin wants
Pedigree = non_geno_stuff  
  
## get cohort size each year
a = table(Pedigree$YearBorn)
BigCohorts[R,] = a

## get total Nb each year

Nbmoms = 1:YearsSamp
Nbdads = 1:YearsSamp

for (j in 1:YearsSamp)  {
year = Born[j]
cohort = subset(Pedigree,Pedigree$YearBorn == year)
RSmoms = table(cohort$Mom)
SSmoms = sum(RSmoms^2)
  if(SSmoms > sum(RSmoms))  { Nbmoms[j] = (sum(RSmoms)-1)/(SSmoms/sum(RSmoms)-1) }
  else {Nbmoms[j] = 99999}
RSdads = table(cohort$Dad)
SSdads = sum(RSdads^2)
  if(SSdads > sum(RSdads))  {Nbdads[j] = (sum(RSdads)-1)/(SSdads/sum(RSdads)-1) }
  else {Nbdads[j] = 99999}
}

BigNbmoms[R,] = Nbmoms
BigNbdads[R,] = Nbdads
BigNb[R,] = 4*Nbmoms*Nbdads/(Nbmoms+Nbdads)

######### subsample cohorts of offspring
Pedigree2 = Pedigree[1,]

for (j in 1:YearsSamp)  {
  year = Born[j]
  cohort = subset(Pedigree,Pedigree$YearBorn == year)
  sampled = cohort[sample(nrow(cohort),SampleSize,replace=F),]
Pedigree2 = rbind(Pedigree2,sampled)
}  # end for j
Pedigree2 = Pedigree2[-1,]

sibs = GetSibs(Pedigree2)

## get age gaps between sibs
Gap = abs(sibs$Born2 - sibs$Born1)

MHSP = subset(sibs,sibs$Parent == "mom")
PHSP = subset(sibs,sibs$Parent == "dad")
Halfs = rbind(MHSP,PHSP)
FSP = subset(sibs,sibs$Parent == "both")

##get within cohort sibs
withinhalf = subset(Halfs,Halfs$Born1 == Halfs$Born2)
 for(j in 1:YearsSamp)  {
      bit = subset(withinhalf,withinhalf$Born1 == Born[j])
      BigSibsHalf[j,j] = nrow(bit) 
   }  # end for j   
withinfull = subset(FSP,FSP$Born1 == FSP$Born2)
    for(j in 1:YearsSamp)  {
      bit2 = subset(withinfull,withinfull$Born1 == Born[j])
      BigSibsFull[j,j] = nrow(bit2) 
      }## end for j
    
##get arcoss cohort sibs
acrosshalf = subset(Halfs,Halfs$Born1 != Halfs$Born2)
 for(j in 1:(YearsSamp-1))  {
 for(k in (j+1):YearsSamp)  {
    bit = subset(acrosshalf,acrosshalf$Born1 == Born[j] & acrosshalf$Born2 == Born[k])
    BigSibsHalf[j,k] = nrow(bit)
    }} # end for j,k
acrossfull = subset(FSP,FSP$Born1 != FSP$Born2)
 for(j in 1:(YearsSamp-1))  {
 for(k in (j+1):YearsSamp)  {
    bit = subset(acrossfull,acrossfull$Born1 == Born[j] & acrossfull$Born2 == Born[k])
    BigSibsFull[j,k] = nrow(bit)
    }} # end for j,k

BigGapHalf[,,R] = BigSibsHalf 
BigGapFull[,,R] = BigSibsFull   

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

TotSibs = rep(0,YearsSamp)
for (j in 1:YearsSamp)  {
for (k in j:YearsSamp)  {
  X = 1 + k-j
  TotSibs[X] = TotSibs[X] + MeanGapHalf[j,k]
  }}
  
##head(BigNb)
##head(WangNb)
#head(BigGapHalf)
#head(BigGapFull)
#BigSibsSame[,,1:3]
#colMeans(BigGapHalf)

### harmonic means
1/mean(1/WangNb)
1/mean(1/BigNb)
1/mean(1/Nhat)

sum(BigGapFull,na.rm=T)

t(NData)
t(WangData)

GapCompare = rowSums(ComparisonsBetween,na.rm=T)
ESibs = GapCompare*effectivesurvival*2/AdultN
TotSibs
ESibs

cohort_size
SampleSize


