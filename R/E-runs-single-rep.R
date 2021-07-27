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
DownSampleSize = as.integer(snakemake@params$DownSampleSize)
rep_num = as.integer(snakemake@params$rep_num)

outfile <- snakemake@output[[1]]

set.seed(cohort_size * SampleSize + rep_num)

# here are values for testing
#cohort_size <- 2000
#SampleSize <- 125  # roughly the number of indviduals to sample each year
#DownSampleSize <- 625  # the total number of individuals to be left with across all years, after downsampling


# set the population life history
SPD <- species_1_life_history

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


pair_nums <- lapply(
  rep_num,
  function(R)  {

    spip_dir <- run_spip(
      pars = SPD
    )

    # now read that in and find relatives within the one-generation pedigree
    slurped <- slurp_spip(spip_dir, 1)

    # get the average number of reproductively active adults
    mean_num_adults <- slurped$census_prekill %>%
      filter(age >= 3) %>%
      group_by(year) %>%
      summarise(Na_tot = sum(male) + sum(female)) %>%
      summarise(mean_Na = mean(Na_tot)) %>%
      pull(mean_Na)

    # get the related pairs
    crel <- compile_related_pairs(slurped$samples)

    # downsample to DownSample total samples, then to 1/sqrt(2) of that and 1/2 of that
    CNTS <- lapply(DownSampleSize * c(1, 1/sqrt(2), 1/2), function(ds) {
      dsamp <- downsample_pairs(slurped$samples, crel, ds)

      # retain only the cross-cohort half siblings
      xc_half_sibs_and_PO <- dsamp$ds_pairs %>%
        filter(
          (dom_relat == "Si" & born_year_1 != born_year_2 &  max_hit == 1) |
          dom_relat == "PO"
        )

      xc_half_sibs_and_PO %>%
        count(dom_relat,conn_comp) %>%
        rename(npairs_in_conn_comp = n) %>%
        count(dom_relat, npairs_in_conn_comp) %>%
        mutate(tot_pairs = n * npairs_in_conn_comp) %>%
        mutate(
          cohort_size = cohort_size,
          SampleSize = SampleSize,
          DownSampleTot = ds,
          mean_num_adults = mean_num_adults
        )

    }) %>%
      bind_rows()

    # return  a tibble with the rep_num
    CNTS %>%
      mutate(rep = R) %>%
      select(rep, everything())

  }) %>%
  bind_rows()


write_rds(pair_nums, file = outfile)

