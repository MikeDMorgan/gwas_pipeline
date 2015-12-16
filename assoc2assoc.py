'''
assoc2assoc.py filter, transform and process results from genome-wide analyses
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Format, manipulation and processing of genome-wide association results

Usage
-----

.. Example use case

Example::

   python assoc2assoc.py

Type::

   python assoc2assoc.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import PipelineGWAS as gwas


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--task", dest="task", type="choice",
                      choices=["get_hits"],
                      help="task to perform")

    parser.add_option("--p-threshold", dest="p_threshold", type="float",
                      help="threshold for association p-value, below "
                      "which results will be output")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # if the input is a list of files, split them
    infile = argv[-1]
    infiles = infile.split(",")
    if len(infiles) > 1:
        results = gwas.GWASResults(assoc_file=infiles)
    elif len(infiles) == 1:
        results = gwas.GWASResults(assoc_file=infile)
    else:
        raise IOError("no input files detected, please specifiy association "
                      "results files as the last command line argument")

    if options.task == "get_hits":
        for snp in results.getHits(float(options.p_threshold)):
            options.stdout.write("%s\n" % snp)
    else:
        pass

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
