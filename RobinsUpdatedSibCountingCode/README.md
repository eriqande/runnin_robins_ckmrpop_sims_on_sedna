README: CKMRpop runs with corrected sib counting
================
Eric C. Anderson

## GitHub Documents

Robin discovered a bug in some of his code and so wanted to re-run some
of the simulations. Everything to do so is contained in this directory,
`RobinsUpdatedSibCountingCode`.

He only needed the really large sizes re-done. I took his updated code
and inserted it into the scheme I had originally made for parallelizing
the runs so I could run it quickly on the cluster.

To launch this on the cluster, you must have Snakemake in your
environment. From the directory in which this README resides I did this:

``` sh
snakemake --use-envmodules --profile ../slurm_runnin_robin  --jobs 400
```

After this was done I tarred the summarized results up into
`RobinsUpdatedSibCountingCode/summarized-full-output.tar.gz`.

The version of the code that I used to run the simulations can be found
in commit `a865aeef323ff32bea1a7b694318827140a9c747`. Or you can grab it
by the tag `Phi=1`.

## Doing a set of runs with Phi=10

I had missed Robin’s request to run these at values of Phi equal to 1
and 10.

Rather than code those different runs up programmatically with different
wildcards/parameters, I am just going to re-run what I did, having
changed Phi to 10 in the R code and also changing the file output names
to reflect that.

That is at tag `Phi=10`.
