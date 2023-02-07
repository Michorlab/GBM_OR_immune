# module: call CNVs with cnvkit
_cnvkit_threads=32


def cnvkit_targets(wildcards):
    """ Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        r = config['runs'][run]
        tmr = r[1]
        ls.append("analysis/cnvkit/%s/%s_recalibrated.cnr" % (run,tmr))
        ls.append("analysis/cnvkit/%s/%s_recalibrated.cns" % (run,tmr))
    return ls

rule cnvkit_all:
    input:
        cnvkit_targets
    benchmark: "benchmarks/cnvkit/cnvkit_all.txt"

def cnvkitInputFn(wildcards):
    #run = config['runs'][wildcards.run]
    #tmr = run[1]
    tmr = wildcards.tmr
    #return tmr's recalibrated bam file
    ls = ["analysis/align/%s/%s_recalibrated.bam" % (tmr, tmr)]
    return ls


rule copynumber_ratios:
    """call copynumber ratios using cnvkit pipeline"""
    input:
        cnvkitInputFn
    output:
        cnr = "analysis/cnvkit/{run}/{tmr}_recalibrated.cnr"
    threads: _cnvkit_threads
    group: "cnvkit"
    params:
        fasta = config['genome_fasta'],
        refflat =  config['cnvkit'],
        output_dir=lambda wildcards: "analysis/cnvkit/%s" % wildcards.run,
    log: "analysis/logs/cnvkit/{run}/{tmr}.cnv.ratios.log"
    benchmark:
        "benchmarks/cnvkit/{run}/{tmr}.cnvkit.txt"
    shell:
        """cnvkit.py  batch {input} -m wgs  -f {params.fasta} --annotate {params.refflat}  -p 16 --scatter --diagram -d {params.output_dir}"""


rule copynumber_segment:
    """call copynumber segment from copynumber ratios"""
    input:
        cnr = "analysis/cnvkit/{run}/{tmr}_recalibrated.cnr"
    output:
        cns = "analysis/cnvkit/{run}/{tmr}_recalibrated.cns"
    group: "cnvkit"
    params:
        output_dir=lambda wildcards: "analysis/cnvkit/%s" % wildcards.run,
    log: "analysis/logs/cnvkit/{run}/{tmr}"
    benchmark:
       "benchmarks/cnvkit/{run}/{tmr}.segment.txt"
    shell:
        """cnvkit.py segment {input.cnr} -o {output.cns}"""
