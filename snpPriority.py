'''
snpPriority.py - score SNPs based on their LD score and SE weighted effect sizes
===============================================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Score SNPs based on their LD score and SE weighted effect sizes from
association analysis.

Usage
-----

.. Example use case

Example::

   python snpPriority.py

Type::

   python snpPriority.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import PipelineGWAS as gwas
import re
import pandas as pd


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--score-method", dest="method", type="choice",
                      choices=["PICS", "LDscore", "ABF", "R2_rank"],
                      help="SNP scoring/prioritisation method to apply.")

    parser.add_option("--database", dest="database", type="string",
                      help="SQL database containing LD information "
                      "in table format. Expects columns SNP_A, "
                      "SNP_B, R2, BP_A and BP_B (Plink --r2 output)")

    parser.add_option("--table-name", dest="table", type="string",
                      help="name of the SQL table containing the LD"
                      "values")

    parser.add_option("--chromosome", dest="chromosome", type="string",
                      help="chromosome to subset the association results "
                      "file on")

    parser.add_option("--ld-threshold", dest="ld_threshold", type="float",
                      help="the threshold of LD above which variants will "
                      "be taken forward.")

    parser.add_option("--rank-threshold", dest="rank_threshold", type="float",
                      help="the threshold in terms of the top n% SNPs to "
                      "output based on the ranking metric. e.g. "
                      "--rank-threshold=0.01 is the top 1% SNPs")

    parser.add_option("--credible-interval", dest="interval", type="float",
                      help="The credible set interval size to generate the "
                      "credible set of SNPs")

    parser.add_option("--prior-variance", dest="prior_var", type="float",
                      help="the prior variance used to weight the SNP "
                      "variance")

    parser.add_option("--fine-map-window", dest="map_window", type="int",
                      help="the region size to included around the index "
                      "SNP as the fine-mapping region.")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    peek = pd.read_table(infile, nrows=5, sep="\s*", header=0)
    try:
        if len(peek["TEST"] != "ADD"):
            clean = False
        else:
            clean = True
    except KeyError:
        clean = True

    if options.method == "LDscore":
        snpscores = gwas.snpPriorityScore(gwas_results=infile,
                                          database=options.database,
                                          table_name=options.table,
                                          chromosome=options.chromosome,
                                          clean=clean)
    elif options.method == "PICS":
        # PICS scores expects the gwas results file to
        # only contain the region of interest, which
        # represents an independent association signal
        snpscores = gwas.PICSscore(gwas_results=infile,
                                   database=options.database,
                                   table_name=options.table,
                                   chromosome=options.chromosome,
                                   clean=clean,
                                   ld_threshold=options.ld_threshold)

        snpscores.columns = ["SNP", "PICS"]

    elif options.method == "R2_rank":
        # rank SNPs based on their LD with the lead
        # SNP, take the top n% SNPs
        snpscores = gwas.LdRank(gwas_results=infile,
                                database=options.database,
                                table_name=options.table,
                                chromosome=options.chromosome,
                                ld_threshold=options.ld_threshold,
                                top_snps=options.rank_threshold,
                                clean=clean)

    elif options.method == "ABF":
        snpscores = gwas.ABFScore(gwas_results=infile,
                                  region_size=options.map_window,
                                  chromosome=options.chromosome,
                                  set_interval=options.interval,
                                  prior_variance=options.prior_var,
                                  clean=clean)

    snpscores.to_csv(options.stdout, index_label="SNP",
                     sep="\t")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
