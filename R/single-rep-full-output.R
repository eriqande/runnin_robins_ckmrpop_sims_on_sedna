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
cohort_size <- as.integer(snakemake@params$cohort_size)
SampleSize = as.integer(snakemake@params$SampleSize)
rep_num = as.integer(snakemake@params$rep_num)

# here are values for testing
cohort_size <- 2000
SampleSize <- 125
rep_num <- 1

my_seed = cohort_size * rep_num




#### Set some things for single runs ####
# We set NReps to 2 to make sure all the arrays get allocated and work
# the way they are supposed to.  Then we only cycle over R from 1 to NReps - 1
# (i.e. 1 to 1).
NReps <- 2

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

# eric reduces the memory over-allocation here
SPD$`alloc-extra` <- 2

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

for (R in 1:(NReps-1))  {

print(paste0("Replicate = ",R))
flush.console()

BigSibsHalf = matrix(0,YearsSamp,YearsSamp)
BigSibsFull = BigSibsHalf

set.seed(my_seed)
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

# Now we just save the BigSibsHalf and BigSibsFull for summarizing later.
write_rds(
  list(
    bsh = BigSibsHalf,
    bsf = BigSibsFull,
    slurped = slurped,
    Pedigree = Pedigree
  ), file = snakemake@output[[1]]
)
