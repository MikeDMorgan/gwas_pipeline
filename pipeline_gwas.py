##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""===========================
Pipeline template
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

.. Replace the documentation below with your own description of the
   pipeline's purpose

Overview
========

This pipeline computes the word frequencies in the configuration
files :file:``pipeline.ini` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_gwas.py config

Input files
-----------

None required except the pipeline configuration files.

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

.. Add any additional external requirements such as 3rd party software
   or R modules below:

Requirements:

* samtools >= 1.1

Pipeline output
===============

.. Describe output files of the pipeline here

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys
import os
import sqlite3
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import PipelineGWAS as gwas

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.
PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True))


# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.

# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

# ---------------------------------------------------
# load the UKBiobank phenotype data into an SQLite DB
# use this as the main accessor of phenotype data
# for the report and non-genetic analyses
@follows(mkdir("phenotypes.dir"))
@transform("%s/*.tab" % PARAMS['data_dir'],
           regex("%s/(.+).tab" % PARAMS['data_dir']),
           add_inputs(r"%s/\1.r" % PARAMS['data_dir']),
           r"phenotypes.dir/\1.tsv")
def formatPhenotypeData(infiles, outfile):
    '''
    Use the UKBiobank encoding dictionary/R script to
    set the factor levels for phenotype data
    '''

    pheno_file = infiles[0]
    r_file = infiles[1]

    statement = '''
    python /ifs/devel/projects/proj045/gwas_pipeline/pheno2pheno.py
    --task=set_factors
    --R-script=%(r_file)s
    %(pheno_file)s
    > %(outfile)s
    '''

    P.run()

@follows(formatPhenotypeData)
@transform(formatPhenotypeData,
           suffix(".tsv"),
           ".load")
def loadPhenotypes(infile, outfile):
    '''
    load all phenotype data in to an SQLite DB
    '''

    P.load(infile, outfile)

# ---------------------------------------------------
# Specific pipeline tasks
@follows(mkdir("plink.dir"),
         loadPhenotypes)
@collate("%s/*.bgen" % PARAMS['data_dir'],
         regex("%s/(.+)\.(.+)" % PARAMS['data_dir']),
         add_inputs("%s/*.sample" % PARAMS['data_dir']),
         r"plink.dir/\1.bed")
def convertToPlink(infiles, outfiles):
    '''
    Convert from Oxford binary (BGEN) format
    to Plink binary format.  One bed file
    per chromosome - keep the fam files the same
    '''

    job_threads = int(PARAMS['data_threads'])
    job_memory = "10G"
    infiles = ",".join([x for x in infiles[0]])

    log_out = outfiles.split(".")[0]
    out_pattern=outfiles.split(".")[0]

    statement = '''
    python /ifs/devel/projects/proj045/gwas_pipeline/geno2assoc.py
    --program=plink2
    --input-file-format=%(data_format)s
    --phenotypes-file=%(data_phenotypes)s
    --pheno=%(format_pheno)s
    --update-sample-attribute=gender
    --format-parameter=%(format_gender)s
    --method=format
    --format-method=change_format
    --reformat-type=plink_binary
    --output-file-pattern=%(out_pattern)s
    --log=%(log_out)s.log
    --threads=%(job_threads)i
    %(infiles)s
    '''

    P.run()

if PARAMS['candidate_region']:
    @follows(convertToPlink,
             mkdir("candidate.dir"))
    @transform(convertToPlink,
               regex("plink.dir/%s(.+).bed" % PARAMS['candidate_chromosome']),
               add_inputs([r"plink.dir/%s\1.fam" % PARAMS['candidate_chromosome'],
                           r"plink.dir/%s\1.bim" % PARAMS['candidate_chromosome']]),
               r"candidate.dir/%s\1-candidate_region.bed" % PARAMS['candidate_chromosome'])
    def getCandidateRegion(infiles, outfile):
        '''
        Pull out genotyping data on individuals over a
        candidate region for testing.
        '''

        bed_file = infiles[0]
        fam_file = infiles[1][0]
        bim_file = infiles[1][1]
        plink_files = ",".join([bed_file, fam_file, bim_file])

        region = ",".join(PARAMS['candidate_region'].split(":")[-1].split("-"))
        out_pattern = ".".join(outfile.split(".")[:-1])

        statement = '''
        python /ifs/devel/projects/proj045/gwas_pipeline/geno2assoc.py
        --program=plink2
        --input-file-format=plink_binary
        --method=format
        --restrict-chromosome=%(candidate_chromosome)s
        --snp-range=%(region)s
        --format-method=change_format
        --format-parameter=%(format_gender)s
        --update-sample-attribute=gender
        --reformat-type=plink_binary
        --output-file-pattern=%(out_pattern)s
        --log=%(outfile)s.log
        %(plink_files)s
        '''

        P.run()

    @follows(getCandidateRegion)
    @transform(getCandidateRegion,
               regex("candidate.dir/(.+).bed"),
               add_inputs([r"candidate.dir/\1.fam",
                           r"candidate.dir/\1.bim"]),
               r"candidate.dir/\1_assoc.assoc")
    def testCandidateRegion(infiles, outfile):
        '''
        Test the candidate region for association using
        Plink basic association - i.e. not a model-specific
        analysis
        '''

        bed_file = infiles[0]
        fam_file = infiles[1][0]
        bim_file = infiles[1][1]
        plink_files = ",".join([bed_file, fam_file, bim_file])

        out_pattern = ".".join(outfile.split(".")[:-1])

        job_threads = PARAMS['candidate_threads']
        job_memory = PARAMS['candidate_memory']

        statement = '''
        python /ifs/devel/projects/proj045/gwas_pipeline/geno2assoc.py
        --program=plink2
        --input-file-format=plink_binary
        --method=association
        --association-method=assoc
        --genotype-rate=0.01
        --indiv-missing=0.01
        --hardy-weinberg=0.0001
        --min-allele-frequency=0.001
        --output-file-pattern=%(out_pattern)s
        --threads=%(candidate_threads)s
        --log=%(outfile)s.log
        -v 5
        %(plink_files)s
        '''

        P.run()

    @follows(testCandidateRegion)
    @transform(getCandidateRegion,
               regex("candidate.dir/(.+).bed"),
               add_inputs([r"candidate.dir/\1.fam",
                           r"candidate.dir/\1.bim"]),
               r"candidate.dir/\1_conditional.%s" % PARAMS['conditional_model'])
    def conditionalTestRegion(infiles, outfile):
        '''
        Perform association analysis conditional on
        top SNP from previous analysis
        '''

        bed_file = infiles[0]
        fam_file = infiles[1][0]
        bim_file = infiles[1][1]
        plink_files = ",".join([bed_file, fam_file, bim_file])

        out_pattern = ".".join(outfile.split(".")[:-1])

        statement = '''
        python /ifs/devel/projects/proj045/gwas_pipeline/geno2assoc.py
        --program=plink2
        --input-file-format=plink_binary
        --method=association
        --association-method=logistic
        --genotype-rate=0.1
        --indiv-missing=0.1
        --hardy-weinberg=0.0001
        --conditional-snp=rs12924101
        --min-allele-frequency=0.01
        --output-file-pattern=%(out_pattern)s
        --threads=%(pca_threads)s
        --log=%(outfile)s.log
        -v 5
        %(plink_files)s
        '''
        
        P.run()
else:
    pass

@follows(convertToPlink,
         mkdir("grm.dir"))
@transform("plink.dir/*.*",
           regex("plink.dir/(.+).bed"),
           add_inputs([r"plink.dir/\1.fam",
                       r"plink.dir/\1.bim"]),
           r"grm.dir/\1.grm.N.bin")
def makeGRM(infiles, outfiles):
    '''
    Calculate the realised GRM per chromosome from
    Plink binary format files
    No filtering at this stage?
    '''

    job_threads = PARAMS['grm_threads']
    job_memory = "10G"

    bed_file = infiles[0]
    fam_file = infiles[1][0]
    bim_file = infiles[1][1]

    plink_files = ",".join([bed_file, fam_file, bim_file])
    out_pattern = ".".join(outfiles.split(".")[:-4])

    statement = '''
    python /ifs/devel/projects/proj045/gwas_pipeline/geno2assoc.py
    --program=gcta
    --input-file-format=plink_binary
    --method=matrix
    --method-compression=bin
    --matrix-form=grm
    --output-file-pattern=%(out_pattern)s
    --threads=%(grm_threads)s
    %(plink_files)s
    '''

    print out_pattern


@follows(makeGRM,
         mkdir("pca.dir"))
@transform("grm.dir/*.grm.bin",
           regex("grm.dir/(.+).grm.bin"),
           add_inputs([r"grm.dir/\1.grm.N.bin",
                       r"grm.dir/\1.grm.id"]),
           r"pca.dir/\1.eigenvec")
def runPCA(infiles, outfile):
    '''
    Run PCA on GRM(s)
    '''

    job_threads = PARAMS['pca_threads']
    n_file = infiles[0]
    bin_file = infiles[1][0]
    id_file = infiles[1][1]

    grm_files = ",".join([n_file, bin_file, id_file])

    out_pattern = outfiles.split(".")[-2]

    statement = '''
    python /ifs/devel/projects/proj045/gwas_pipeline/geno2assoc.py
    --program=gcta
    --input-file-format=GRM_binary
    --method=pca
    --principal-components=%(pca_components)i
    --threads=%(job_threads)s
    --output-file-pattern=%(out_pattern)s
    %(grm_files)s
    '''

    print grm_files
    
@follows(convertToPlink,
         mkdir("gwas.dir"))
@collate("plink.dir/chr16impv1*",
         regex("plink.dir/(.+).(...+)"),
         r"gwas.dir/test_assoc-\1.assoc")
def test_association(infiles, outfile):
    '''
    Run a test association on chromosome 16
    '''

    job_memory = "5G"
    job_threads = PARAMS['pca_threads']

    infiles = ",".join(infiles)

    statement = '''
    python /ifs/devel/projects/proj045/gwas_pipeline/geno2assoc.py
    --program=plink2
    --input-file-format=plink_binary
    --method=association
    --association-method=assoc
    --genotype-rate=0.01
    --indiv-missing=0.01
    --hardy-weinberg=0.0001
    --min-allele-frequency=0.01
    --output-file-pattern=gwas.dir/test_assoc
    --threads=%(pca_threads)s
    --log=%(outfile)s.log
    -v 5
    %(infiles)s
    '''

    P.run()

@follows(mkdir("dissect.dir"))
@collate("%s/*" % PARAMS['data_dir'],
         regex("%s/(.+).(.+)" % PARAMS['data_dir']),
         r"dissect.dir/test_out")
def test_dissect(infiles, outfile):
    '''
    Test DISSECT for parallelised analysis
    MUST have a separate phenotypes file with no header, but in plink format
    '''

    job_memory = "256G"
    infiles = ",".join(infiles)
    job_threads = 256
    job_queue = "mpi.q"
    cluster_parallel_environment = " mpi "
    statement = '''
    mpirun --np %(job_threads)s
    dissect
    --gwas
    --bfile data.dir/chr22impv1u
    --pheno data.dir/dissect.pheno
    --out dissect.dir/out_file
    '''

    P.run()

# ---------------------------------------------------
# Generic pipeline tasks
@follows()
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
