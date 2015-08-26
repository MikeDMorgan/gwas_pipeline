'''
pheno2pheno.py - format and manipulate phenotype files
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Project specific operations on phenotype files

Usage
-----

.. Example use case

Example::

   python pheno2pheno.py

Type::

   python pheno2pheno.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri
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

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--task", dest="task", type="string",
                      help="task to execute on phenotype file(s)")

    parser.add_option("--R-script", dest="r_script", type="string",
                      help="R script for table reformatting")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]
    if options.task == "set_factors":
        pandas2ri.activate()
        R('''source("%s")''' % options.r_script)
        R('''format <- format_phenotypes("%s")''' % infile)
        pheno_df = pandas2ri.ri2py_dataframe(R["format"])

        pheno_df.to_csv(options.stdout, sep="\t", index_col=None)
    else:
        pass
        

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
