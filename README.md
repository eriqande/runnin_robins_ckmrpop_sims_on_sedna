README
================

Robin gave me a basic script to run multiple reps of spip simulations
and then summarize them. In order to use the full power of the cluster,
I am going to convert this into a system that has two parts:

1.  a script that does a single run of spip and processes the output and
    stores it in an RDS file.
2.  a script that processes all the reps done in item 1 (above) for a
    particular set of parameters, and spits out the results Robin wants.

We will do this within the context of Snakemake and just declare the
parameter values in simple Python lists, and expand them appropriately.

The variables that will become our wildcards/parameters will be:

  - `cohort_size`
  - `SampleSize`
  - `Rep`

We will set the seed for each run using the cohort\_size and the rep
number. Basically, we will just set it to the product of `cohort_size`
and `Rep`.

I have broken Robin’s original script into two separate ones:

  - `R/single_rep.R`
  - `R/summarize_reps.R`

The first runs just a single rep of `run_spip()` and the associated code
to get the numbers of different sibs, etc. It then writes that out to a
file.

The second uses all the saved files of ouput to summarize all the reps.

This is all working with Snakemake via the Snakefile.

## Checking

I ran Robin’s code, largely unmodified save for setting seeds for
reproducibility, in `sib-code-modified-by-eric.R`. I compared the
results for that run to what I got running everything in parallel via
the Snakefile. It was identical. Good.

## Running on SLURM

This should be relatively straightforward, but I want the resource
(memory) requirements to differ between the different cohort sizes. That
should not be difficult.

It sounded like Robin was running out of memory while spip was running.
He said it ran fine for cohort size 128000 for 20 reps, then
mysteriously failed. I suspect some garbage collection problems. I think
Robin’s machine has 32 Gb of RAM. Since I reset the over-alloc parameter
by a factor of 5, this means that cohort size 128000 ought to be OK at
6.4 Gb of RAM. So, let’s just round it to 8, or, actually, 9, which is a
little less than what SEDNA usually allocates for 2 CPUs. For every
smaller cohort size we should be fine with the default 4.7 Gb that it
gets with one core.

For other cohort sizes:

  - 512,000. Let’s shoot for 30 Gb which could fit three runs on a
    standard node.
  - 2,048,000. I will first try 94 Gb which should let it fit on a
    standard node. If that does not work, then I will go for the
    smallest number that works so I can fit as many as possible on the
    himem machines.

We want to restrict our total resource utilization to something
reasonable. Let’s say, we don’t want to use any more than 20 nodes.
That’s 400 cores, but we also don’t want to eat up any more memory
than what we have on 20 standard nodes: 94\*20 Gb = 1880000 Mb. So, we
just need to pass the memory constraint to keep the number of jobs down,
but I will add cpus in there too. We do that in the slurm profile. I’m
not sure if the default slurm profile maps cpus to the resources. Kind
of lame, but not an issue…

I set up a slurm profile here and I modified it to use appropriate
values.

## Report from initial 5 reps

The cohort\_size = 512000 runs took between 20 and 30 minutes. The
cohort\_size = 2048000 runs will take at least 4 times longer, and maybe
more. The bottleneck is tidyverse operations to make the files from
whence to build the C data structures for recursively traversing the
pedigree. Pretty lame that that is the part that takes the most time.

It appears that everything has sufficient memory.

Once the five remaining cohort\_size = 2048000 runs are done I will
launch everything and it should all be done within no more than a couple
of days.
