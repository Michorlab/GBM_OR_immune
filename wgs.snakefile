#!/usr/bin/env python

import os
import sys
import subprocess
import pandas as pd
import json
import yaml

from string import Template


def getRunsCohorts(config):
    """parse metasheet for Run groupings"""
    ret = {}
    metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#', skipinitialspace=True)
    f = metadata.to_csv().split() #make it resemble an actual file with lines
    #SKIP the hdr
    cohorts = {}
    for l in f[1:]:
        tmp = l.strip().split(",")
        #print(tmp)
        ret[tmp[0]] = tmp[1:3]

        if len(tmp) > 3:
            chort=tmp[3] #it's the 4th col
            if chort in cohorts:
                cohorts[chort].append(tmp[0]) #add run to the cohort set
            else:
                #new cohort
                cohorts[chort] = [tmp[0]]

    #Find any cohorts that are singletons and print warning
    delete_these = []
    for c in cohorts:
        if len(cohorts[c]) < 2:
            print("WGS WARNING: cohort group %s ignored because it has too few runs in the set (at least 2 are required).  Please correct your metasheet" % c)
            delete_these.append(c)
    #REMOVE them
    for c in delete_these:
        del cohorts[c] #not valid cohort, remove

    #print(ret)
    config['runs'] = ret
    config['cohorts'] = cohorts if cohorts else None
    #print(config['cohorts'])
    return config

def addCondaPaths_Config(config):
    """ADDS the python2 paths to config"""
    #conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    conda_root=os.environ.get("CONDA_ROOT")
    config['conda_root'] = conda_root
    #config['wes_root'] = "%s/envs/wes" % conda_root # DEPRECATED
    config['vcf_root'] = "%s/envs/vcf" % conda_root


def loadRef(config):
    """Adds the static reference paths found in config['ref']
    NOTE: if the elm is already defined, then we DO NOT clobber the value
    """
    f = open(config['ref'])
    ref_info = yaml.safe_load(f)
    f.close()
    #print(ref_info[config['assembly']])
    for (k,v) in ref_info[config['assembly']].items():
        #NO CLOBBERING what is user-defined!
        if k not in config:
            config[k] = v

#---------  CONFIG set up  ---------------
configfile: "config.yaml"   # This makes snakemake load up yaml into config
config = getRunsCohorts(config)
addCondaPaths_Config(config)
loadRef(config)

#Finally check for 'remote_path' and 'expressions'
if 'remote_path' not in config:
    config['remote_path'] = ""
if 'expression_files' not in config:
    config['expression_files'] = {}
#-----------------------------------------


#------------------------------------------------------------------------------
# TARGETS
#------------------------------------------------------------------------------
def all_targets(wildcards):
    ls = []
    ls.extend(align_targets(wildcards))
    ls.extend(metrics_targets(wildcards))
    ls.extend(recalibration_targets(wildcards))
    ls.extend(cnvkit_targets(wildcards))
    return ls

rule target:
    input:
        targets=all_targets,
        benchmarks="benchmarks.tar.gz"
    message: "Compiling all outputs"

rule tar_benchmarks:
    input: all_targets
    output: "benchmarks.tar.gz"
    shell: "tar -c benchmarks | gzip > {output}"


include: "./modules/align.snakefile"     # common align rules
include: "./modules/metrics.snakefile"   # sentieon basic qc metrics module
include: "./modules/recalibration.snakefile" # sentieon BQSR
include: "./modules/copynumber.snakefile"
