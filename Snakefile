# here are the different cohort and sample sizes and reps to use
csz = [500, 2000, 8000, 32000, 128000, 512000, 2048000, 128000, 512000, 2048000]
ssz = [63,   125,  251,   501,   1003,   2005,    4010,    709,   1418,    2836] 
reps = range(1,101)  # 1:100 

mem_mb_dict = {
	'128000': 9000,
	'512000': 30000,
	'2048000': 94000
}
# a function to return memory requirements given cohort_size wildcard
def mem_mb_needed(wc):
	if(int(wc.cohort_size) < 128000):
		return(4600)
	else:
		return(mem_mb_dict[wc.cohort_size])





rule all:
    input:
        expand("summarized/csz{cohort_size}_ssz{sample_size}.rds", zip, cohort_size = csz, sample_size = ssz)


rule simulate:
    output:
        "single-reps/{cohort_size}/{sample_size}/{rep}.rds"
    params:
    	cohort_size = "{cohort_size}",
    	SampleSize = "{sample_size}",
    	rep_num = "{rep}"
    resources:
    	mem_mb = lambda wc: mem_mb_needed(wc)
    envmodules:
    	"R/4.0.3"
    log:
    	"logs/single-reps/{cohort_size}/{sample_size}/{rep}-log.txt"
    script:
    	"R/single-rep.R"
        


rule summarize_reps:
    input:
        expand("single-reps/{{cohort_size}}/{{sample_size}}/{rep}.rds", rep = reps)
    output:
        "summarized/csz{cohort_size}_ssz{sample_size}.rds"
    envmodules:
    	"R/4.0.3"
    log:
        "logs/summarized/csz{cohort_size}_ssz{sample_size}-log.txt"
    script:
        "R/summarize-reps.R"