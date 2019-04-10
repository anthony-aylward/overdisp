#===============================================================================
# overdisp.py
#===============================================================================

"""Estimate overdispersion parameters"""




# Imports ======================================================================

import json
import os
import subprocess
import tempfile




# Constants ====================================================================

ESTIMATE_LSSE_PARAMETERS_SCRIPT = (
    'library(chenimbalance);library(jsonlite);'
    'counts_frame<-read.table("{filename}",header=TRUE,stringsAsFactors=FALSE);'
    'total<-counts_frame[["coverage"]];'
    'allelic_ratio<-counts_frame[["ref_count"]]/total;'
    'lsse_parameters<-alleledb_beta_binomial('
        'total[!is.na(total)],'
        'allelic_ratio[!is.na(total)],'
        'r_by=0.025'
    ');'
    'cat(toJSON(lsse_parameters))'
)

ESTIMATE_NULL_PARAMETERS_SCRIPT = (
    'library(callimbalance);library(jsonlite);'
    'null_parameters<-estimate_null_parameters('
        'read.table("{filename}",header=TRUE,stringsAsFactors=FALSE),'
        'minimum_coverage={minimum_coverage},'
        'n_breaks={n_breaks},'
        'spline_order={spline_order},'
        'cores={processes}'
    ');'
    'cat(toJSON(null_parameters))'
)




# Functions ====================================================================

def estimate_lsse_parameters(counts_frame, temp_dir=None):
    with tempfile.NamedTemporaryFile(dir=temp_dir) as temp_counts:
        temp_counts_name = temp_counts.name
    counts_frame.to_csv(temp_counts_name, sep='\t', na_rep='NA', index=False)
    with subprocess.Popen(
        (
            'Rscript', '-e', ESTIMATE_LSSE_PARAMETERS_SCRIPT.format(
                filename=temp_counts_name
            )
        ),
        stdout=subprocess.PIPE
    ) as r:
        parameters = json.loads(r.communicate()[0])
    os.remove(temp_counts_name)
    return parameters


def estimate_null_parameters(
    counts_frame,
    minimum_coverage=10,
    n_breaks=11,
    spline_order=4,
    processes=1,
    temp_dir=None
):
    with tempfile.NamedTemporaryFile(dir=temp_dir) as temp_counts:
        temp_counts_name = temp_counts.name
    counts_frame.to_csv(temp_counts_name, sep='\t', na_rep='NA', index=False)
    with subprocess.Popen(
        (
            'Rscript', '-e', ESTIMATE_NULL_PARAMETERS_SCRIPT.format(
                filename=temp_counts_name,
                minimum_coverage=minimum_coverage,
                n_breaks=n_breaks,
                spline_order=spline_order,
                processes=processes
            )
        ),
        stdout=subprocess.PIPE
    ) as r:
        parameters = json.loads(r.communicate()[0])
    os.remove(temp_counts_name)
    return parameters
