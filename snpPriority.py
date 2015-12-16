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

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    peek = pd.read_table(infile, nrows=5, sep=None, header=0)
    if len(peek["TEST"] != "ADD"):
        clean = False
    else:
        clean = True

    snpscores = gwas.snpPriorityScore(gwas_results=infile,
                                      database=options.database,
                                      table_name=options.table,
                                      chromosome=options.chromosome,
                                      clean=clean)

    snpscores.to_csv(options.stdout, index_label="SNP",
                     sep="\t")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
