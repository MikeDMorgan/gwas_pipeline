'''
assoc2plot.py - script for post-GWAS plotting
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Post-gwas plotting, for visualization and QC

Usage
-----

.. Example use case

Example::

   python assoc2plot.py

Type::

   python assoc2plot.py --help

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

    parser.add_option("--plot-type", dest="plot_type", type="choice",
                      choices=["manhattan", "qqplot"],
                      help="plot type to generate")

    parser.add_option("--resolution", dest="resolution", type="choice",
                      choices=["genome_wide", "chromosome",
                               "fine_map"],
                      help="the resolution of plotting, wether the plot "
                      "depicts the whole genome, a single chromosome or "
                      "a specific locus")

    parser.add_option("--save-path", dest="save_path", type="string",
                      help="path and filename to save image to")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    results = gwas.GWASResults(assoc_file=infile)
    if options.plot_type == "manhattan":
        results.plotManhattan(resolution=options.resolution,
                              save_path=options.save_path)

    else:
        pass


    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
