#########################################################################
#########################################################################
# Classes for handling genome-wide association input and output files, ##
# analysis and qc programs, and post-hoc analyses                      ##
#########################################################################
#########################################################################

import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import numpy as np
import pandas as pd
import re
import random
import os
import rpy2.robjects as ro
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri as py2ri
from rpy2.robjects.packages import importr
# set matplotlib non-interactive backend to Agg to
# allow running on cluster
import matplotlib
#matplotlib.use("Qt4Agg")
import matplotlib.pyplot as plt
from ggplot import *
import collections

class FileGroup(object):
    '''
    An object for holding, formatting and processing files for genome-wide
    association analysis including compressed and binary files

    File types supported:
    * plink - .ped and .map files
    * plink binary - .bim, .fam. and .bed files
    * variant call format - .vcf and .bcf (including gzipped vcf)
    * Oxford format - .gen or .bgen with matched sample text file (must
                      be .sample)
    * GRM_binary - genetic relationship matrix calculated in an appropriate
      program in binary format.  File suffixes are *.grm.bin, *.grm.N.bin
      and *.grmid
    * GRM_gz - previously calcualted gzip compressed GRM, file suffixes
      are *.grm.gz and *.grm.id

    Phenotypes are assumed to be contained in the relevant files, if not
    then an additional phenotypes files can be included using the
    `phenotypes` argument.  Covariate files (if different from the phenotypes
    file) can also be included in the instantiation of a :FileGroup:
    object using the `covarite_files` argument.

    Only the `files` and `file_format` arguments are required.

    Genotype data are assumed to be raw genotype calls.  This can be modified
    using the `genotype_format` argument upon instantiation.  Values allowed
    are:
    * calls - standard bi-allelic genotype calls, i.e. AA, AB, BB
    * imputed_call - discrete genotype calls from imputed data,
                     essentially treated the same as ``calls``
    * genotype_prob - posterior probabilities for each genotype class,
                      i.e. 0.88 0.07 0.05 corresponding to homozygote
                      reference, heterozygote then homozygote rare allele.
    '''

    # Defaults for file formats
    ped_file = None
    map_file = None
    bim_file = None
    fam_file = None
    bed_file = None
    sample_file = None
    gen_file = None
    bgen_file = None
    vcf_file = None
    bcf_file = None

    def __init__(self, files, file_format, phenotypes=None,
                 genotype_format="calls", covariate_files=None):

        self.files = files
        self.file_format = file_format
        self.pheno_file = phenotypes
        self.genotype_format = genotype_format
        self.covariate_files = covariate_files
        self.set_file_prefix(files)

    def set_file_prefix(self, infiles):
        '''Get file prefixes from input files.  These are used across all
        file formats, e.g. myfile.bed, myfile.bim, myfile.fam name=myfile.
        Only use periods, '.' to denote file suffixes. use hyphens and
        underscores for separating file names.

        Set these to the appropriate attributes.
        '''

        file_prefixes = set()

        for f in infiles:
            # get all input file prefixes
            if len(f.split("/")) > 1:
                g = f.split("/")[-1]
                fdir = f.split("/")[:-1]
                fdir = "/".join(fdir)
                ffile = fdir + "/" + g.split(".")[0]
                file_prefixes.add(ffile)
            else:
                file_prefixes.add(f.split(".")[0])

        # if only prefix then use this for all data files
        if len(file_prefixes) == 1:
            self.name = [xf for xf in file_prefixes][0]
        else:
            # if there are multiple prefixes then use separate
            # flags for file inputs
            self.name = None

        # define file types by their suffix instead
        if self.file_format == "plink":
            self.ped_file = [pf for pf in infiles if re.search(".ped",
                                                               pf)][0]
            self.map_file = [mf for mf in infiles if re.search(".map",
                                                               mf)][0]

            # check files exist (i.e. are not the default None values)
            try:
                assert self.ped_file
            except AssertionError:
                raise ValueError(".ped file is missing, please "
                                 "specify")
            try:
                assert self.map_file
            except AssertionError:
                raise ValueError(".map file is missing, please "
                                 "specify")

        elif self.file_format == "plink_binary":
            self.fam_file = [ff for ff in infiles if re.search(".fam",
                                                               ff)][0]
            self.bim_file = [fb for fb in infiles if re.search(".bim",
                                                              fb)][0]
            self.bed_file = [bf for bf in infiles if re.search(".bed",
                                                               bf)][0]
            # check files exist (i.e. are not the default None values)
            try:
                assert self.fam_file
            except AssertionError:
                raise ValueError(".fam file is missing, please "
                                 "specify")
            try:
                assert self.bim_file
            except AssertionError:
                raise ValueError(".bim file is missing, please "
                                 "specify")
            try:
                assert self.bed_file
            except AssertionError:
                raise ValueError(".bed file is missing, please "
                                 "specify")

        elif self.file_format == "oxford":
            self.gen_file = [gf for gf in infiles if re.search(".gen",
                                                               gf)][0]
            self.sample_file = [sf for sf in infiles if re.search(".sample",
                                                                  sf)][0]
            # check files exist (i.e. are not the default None values)
            try:
                assert self.gen_file
            except AssertionError:
                raise ValueError(".gen file missing, please "
                                 "specify")
            try:
                assert self.sample_file
            except AssertionError:
                raise ValueError(".sample file missing, please "
                                 "specify")

        elif self.file_format == "oxford_binary":
            self.bgen_file = [bg for bg in infiles if re.search(".bgen",
                                                                bg)][0]
            self.sample_file = [sf for sf in infiles if re.search(".sample",
                                                                  sf)][0]
            # check files exist (i.e. are not the default None values)
            try:
                assert self.bgen_file
            except AssertionError:
                raise ValueError(".bgen file is missing, please "
                                 "specify")
            try:
                assert self.sample_file
            except AssertionError:
                raise ValueError(".sample file is missing, please "
                                 "specify")

        elif self.file_format == "vcf":
            self.vcf_file = [vf for vf in infiles if re.search(".vcf",
                                                               vf)][0]
            # check files exist (i.e. are not the default None values)
            try:
                assert self.vcf_file
            except AssertionError:
                raise ValueError(".vcf file is missing, please "
                                 "specify")

        elif self.file_format == "bcf":
            self.bcf_file = [bv for bv in infiles if re.search(".bcf",
                                                               bv)][0]
            # check files exist (i.e. are not the default None values)
            try:
                assert self.bcf_file
            except AssertionError:
                raise ValueError(".bcf file is missing, please "
                                 "specify")

        elif self.file_format == "GRM_binary":
            self.id_file = [ig for ig in infiles if re.search(".grm.id",
                                                              ig)][0]
            self.n_file = [gn for gn in infiles if re.search(".grm.N",
                                                             gn)][0]
            self.bin_file = [gb for gb in infiles if re.search(".grm.bin",
                                                               gb)][0]
            # check files exits
            try:
                assert self.id_file
            except AssertionError:
                raise ValueError("GRM ids file is missing, please "
                                 "specify")
            try:
                assert self.n_file
            except AssertionError:
                raise ValueError("grm.N file is missing, please "
                                 "specify")
            try:
                assert self.bin_file
            except AssertionError:
                raise ValueError("GRM binary file is missing, please "
                                 "specify")

        elif self.file_format == "GRM_plink":
            self.id_file = [ig for ig in infiles if re.search(".rel.id",
                                                              ig)][0]
            self.rel_file = [gn for gn in infiles if re.search(".rel.N.bin",
                                                             gn)][0]
            # check files exits
            try:
                assert self.id_file
            except AssertionError:
                raise ValueError("GRM ids file is missing, please "
                                 "specify")
            try:
                assert self.rel_file
            except AssertionError:
                raise ValueError("rel.N file is missing, please "
                                 "specify")


    def set_phenotype(self, pheno_file=None, pheno=1):
        '''
        Set the phenotype for a set of individuals
        using an external phenotypes file.

        Default is to use the (n+2)th column, designated
        as pheno 1.
        '''

        if type(pheno) == int:
            pheno = str(pheno)
        elif type(pheno) == str:
            pass
        else:
            raise AttributeError("Type of pheno unknown. "
                                 "Must be str or int.")
        self.pheno_file = pheno_file
        self.pheno = pheno


class GWASProgram(object):
    '''
    A base level object for programs designed to perform genome-wide
    association analysis and operate on genome-wide genotyping data.

    [INSERT PROPER DOCSTRING - see style guide]
    '''

    def __init__(self, executable=None, required_format=None):
        self.executable = executable
        self.require_format = required_format

    def program_call(self, infiles, outfile):
        '''build a statement to perform genome-wide
        analysis using infiles
        '''

        return ""

    def postprocess(self, infiles, outfile):
        '''collect and process output files from
        program - format for Result class objects'''

        return ""

    def build(self, infiles, outfile):
        '''run analysis program'''

        cmd_program = self.program_call(infile, outfile)
        cmd_postprocess = self.postprocess(infiles, outfile)

        if cmd_postprocess:
            cmd_postprocess = cmd_postprocess.strip().endswith(";")
            assert cmd_postprocess
        else:
            pass

        statement = " checkpoint; ".join((cmd_program,
                                          cmd_postprocess))

        return statement


class GCTA(GWASProgram):
    '''
    GCTA is designed for computing genetic relationship matrices, linear
    mixed model analyses and phenotype estimation/prediction.
    It can also perform SNP-wise GWAS.

    Files MUST be in Plink binary format
    '''

    def __init__(self, files, options=None, settings=None,
                 design=None):
        self.infiles = files
        self.options = options
        self.settings = settings
        self.design = design
        self.executable = "gcta64"
        self.statement = {}
        self.filters = []

    def program_call(self, infiles, outfile):
        '''build GCTA call statement on infiles'''

        statement = []
        statement.append(self.executable)

        if infiles.name:
            inputs = self._build_single_file_input(infiles,
                                                   infiles.file_format)
            statement.append(inputs)
        else:
            raise AttributeError("Files must be in binary plink format "
                                 "or as a GRM to use GCTA.  Please "
                                 "convert and try again.")

        if infiles.pheno_file:
            statement.append(" --pheno %s --mpheno %s " % (infiles.pheno_file,
                                                           infiles.pheno))
        else:
            pass

        self.statement["program"] = " ".join(statement)

    def _build_single_file_input(self, infiles, file_format):
        '''internal function only. Use it to construct the
        file input flags with --file, --bfile or --data
        '''

        statement = None

        if file_format == "plink":
            statement = " --file %s " % infiles.name
        elif file_format == "plink_binary":
            statement = " --bfile %s " % infiles.name
        elif file_format == "oxford" or file_format == "oxford_binary":
            statement = " --data %s" % infiles.name
        elif file_format == "GRM_binary" or file_format == "GRM_plink":
            statement = " --grm %s " % infiles.name
        else:
            raise AttributeError("file format is not defined or recognised."
                                 "Please define the input corectly when "
                                 "instantiating a FileGroup object")

        return statement

    def PCA(self, n_pcs="20"):
        '''
        Perform PCA analysis on previosly generated GRM, output the number n
        principal componets, default = 20
        '''

        self._run_tasks(pca=n_pcs)

    def filter_genotypes(self, filter_type, filter_value):
        '''
        * chromosome - exclude all variants not on the specified chromosome(s).
          [str/list]
        * autosome_number - for non-human species, the number of chromosomes to
          be considered autosomes
        * exclude_snps - text file list of variant IDs to exclude from analysis.
          [file]
        * extract - text file list of variant IDs to include in analysis,
          ignores all others. [file]
        * min_allele_frequency - only include SNPs with cohort/case allele
          frequency above this threshold. [float]
        * max_allele_frequency - include all SNPs with a MAF equal to or below
          this value. [float]
        '''

        if filter_type == "chromosome":
            self._construct_filters(chromosome=filter_value)
        elif filter_type == "autosome_number":
            self._construct_filters(autosome_number=filter_value)
        elif filter_type == "exclude_snps":
            self._construct_filters(exclude_snps=filter_value)
        elif filter_type == "extract":
            self._construct_filters(extract=filter_value)
        elif filter_type == "min_allele_frequency":
            self._construct_filters(min_allele_frequency=filter_value)
        elif filter_type == "max_allele_frequency":
            self._construct_filters(max_allele_frequency=filter_value)

    def _construct_filters(self, **kwargs):
        '''
        Add filter to each GCTA run.

        The filters accepted are defined below.  These are input as keyword
        arguments supported by this function.

        * min_allele_frequency - only include SNPs with cohort/case allele
          frequency above this threshold. [float]
        * max_allele_frequency - include all SNPs with a MAF equal to or below
          this value. [float]
        * keep - keep individuals with matching individual and family IDs.
          [file]
        * remove - remove all individuals with matching individual and family
          IDs. [file]
        * extract - text file list of variant IDs to include in analysis,
          ignores all others. [file]
        * exclude - text file list of variant IDs to exclude from analysis.
          [file]
        * chromosome - exclude all variants not on the specified chromosome(s).
          [str/list]
        * autosome - exclude all non-place and non-autosomal variants.
          [boolean]
        * covariates_file - specify the covariates file with family and
          individual IDs in the first two columns.  Covariates are in the
          (n+2)th column. Only used in conjunction with `covariate_filter`.
          [file]
        * covariate_filter - covariate columns value to filter on.  Can be
          used with non-numeric values to filter out individuals with
          covariate =/= `covariate_filter` value. [str/int/float]
        * covariate_column - column number to apply filtering to if more
          than one covariate in the file. [int]
        * update_gender - provide gender information in a separate text
          file. [file]
        * grm_threshold - remove one of a pair of individuals with
          estimated relatedness greater than this value.
        * ld_significance - p-value threshold for regression test
          of LD significance
        * genotype_call - GenCall score cut-off for calling raw
          genotypes into Plink PED format
        * meta_pval - p-value threshold cut-off for conditional
          and joint genome-wide analysis
        * cojo_window - distance in kb beyond wich SNPs this
          distance apart are assumed to be in linkage equilibrium
        * cojo_collinear - multiple regression R^2 on selected SNPs
          value above which the testing SNP will not be selected.
        * cojo_inflation - adjust COJO analysis test statistics
          for genomic control. [boolean]
        * reml_iterations - maximum number of iterations to use
          during reml analysis.  Default is 100. [int]
        '''

        statement = []

        # map of keyword arguments recognised to Plink2 filtering flags
        filter_map = {"min_allele_frequency": " --maf %s ",
                      "max_allele_frequency": " --max-maf %s ",
                      "keep": " --keep %s ",
                      "remove": " --remove %s ",
                      "extract": " --extract %s ",
                      "exclude": " --exclude %s ",
                      "chromosome": " --chr %s ",
                      "autosome": " --autosome ",
                      "autosome_number": " --autosome-num %s ",
                      "grm_threshold": " --grm-cutoff %s ",
                      "ld_significance": " --ls-sig %s ",
                      "genotype_call": " --gencall %s ",
                      "meta_pval": " --cojo-p %s ",
                      "cojo_window": " --cojo-wind %s ",
                      "cojo_collinear": " --cojo-collinear %s ",
                      "cojo_inflation": " --cojo-gc ",
                      "reml_iterations": " --reml-maxit %s "}

        # compile all filters together, checking for dependencies.
        # use a mapping dictionary to extract the relevant flags and
        # combinations to use.
        filters = []
        filter_dict = {}
        for key, value in kwargs.iteritems():
            filter_dict[key] = value

        for each in filter_dict.keys():
            try:
                assert filter_map[each]
                # check for data type <- behaviour is type dependent
                if type(filter_dict[each]) == 'bool':
                    filters.append(filter_map[each])
                else:
                    filter_val = filter_dict[each]
                    filters.append(filter_map[each] % filter_val)

            except KeyError:
                E.warn("%s filter not recognised, please see "
                       "documentation for allowed filters" % each)
                pass
                
        self.filters.append(" ".join(filters))
        self.statement["filters"] = " ".join(self.filters)

    def _run_tasks(self, parameter=None, **kwargs):
        '''
        The principal functions of GCTA revolve around GRM estimation
        and variance components analysis, such as REML estimation of
        heritability and variance components, BLUP and phenotype prediciton.

        It can also be used to do PCA and conditional and joint GWAS.

        Tasks
        -----
        * pca - perform principal components analysis on a GRM
        * greml - perform restricted maximum likelihood analysis
          for estimation of variance components
        * estimate_ld - estimate the linkage disequilibrium structure
          over the genomic regions specified
        * simulate_gwas - simulate genome-wide association data based
          on observed genotype data
        * cojo - conditional and joint genome-wide association
          analysis across SNPs and covariates
        * bivariate_reml - perform GREML on two traits, either both
          binary, both quantitative or one of each
        * lmm - perform a linear mixed model based association analysis
        '''

        statement = []

        # set up a dictionary of recognised tasks with key word argument
        # values as further dictionaries. Use the parameter argument
        # to pass arguments by value to string formatting

        # put all of the other tasks as options in the calling function

        task_map = {"pca": " --pca %s ",
                    "greml": {"standard": " --reml ",
                              "priors": " --reml --reml-priors %s ",
                              "reml_algorithm": " --reml --reml-alg %s ",
                              "unconstrained": " --reml --reml-no-constrain ",
                              "GxE": " --reml --gxe %s ",
                              "LRT": " --reml --reml-lrt %s ",
                              "BLUP_EBV": " --reml --reml-pred-rand ",
                              "snpBLUP": " --blup-snp "},
                    "estimate_ld": " --ld %s ",
                    "simulate_gwas": {"quantitative": " --simu-qt ",
                                      "case_control": " --simu-cc %s %s "},
                    "cojo": {"stepwise": " --cojo-file %s --cojo-slct ",
                             "no_selection": " --cojo-file %s --cojo-joint ",
                             "snp_conditional": " --cojo-file %s --cojo-cond %s "},
                    "bivariate_reml": {"standard": " --reml-bivar %s ",
                                       "no_residual": " --reml-bivar %s --reml-bivar-nocove ",
                                       "fixed_cor": " --reml-bivar %s --reml-bivar-lrt-rg %s "},
                    "lmm": {"standard": " --mlma ",
                            "loco": " --mlma-loco ",
                            "covar": " --mlma-no-adj-covar "}}

        for task, value in kwargs.iteritems():
            # check for PCA first as it is not nested in task_map
            if task == "pca":
                try:
                    state = task_map[task] % value
                    statement.append(state)
                except TypeError:
                    statement.append(task_map[task])
                statement.append
            # LD estimation is likewise not nested
            elif task == "estimate_ld":
                try:
                    state = task_map[task] % value
                    statement.append(state)
                except TypeError:
                    raise IOError("no SNP file list detected")
            elif task != "parameter":
                try:
                    # sub_task is a nested dictionary
                    sub_task = task_map[task]
                    try:
                        assert sub_task[value]
                        try:
                            # some tasks do not contain task values for the
                            # parameter argument - catch these with the TypeError
                            # exception
                            statement.append(sub_task[value] % parameter)

                            # the default for parameter is None, check this is appropriate
                            if not parameter:
                                E.warn("Parameter value is set to NoneType. "
                                       "Please check this is an appropriate value "
                                       "to pass for this task")
                            else:
                                pass
                        except TypeError:
                            statement.append(sub_task[value])
                    except KeyError:
                        raise KeyError("% Task not recognised, see docs for details of "
                                       "recognised tasks" % task)                    
                except KeyError:
                    raise KeyError("Task not recognised, see docs for details of "
                                   "recognised tasks")
            else:
                pass

        self.statement["tasks"] = " ".join(statement)

    def genetic_relationship_matrix(self, compression="binary", metric=None,
                                    shape="square", options=None):
        '''
        Calculate the estimated genetic relationship matrix from
        genotyping data

        * estimate_grm - estimate the realized genetic relationship
          matrix between individuals from genotyping data
        '''

        mapf = {"binary": " --make-grm-bin ",
                "gzip": " --make-grm-gz ",
                "no_compress": " --make-grm ",
                "X_chr": " --make-grm-chr ",
                "X_chr_gz": " --make-grm-gz ",
                "inbreeding": " --ibc "}

        if options == "X_chr":
            if compression == "gz":
                state = mapf["X_chr_gz"]
            else:
                state = mapf["X_chr"]
        elif options == "inbreding":
            state = mapf["inbreeding"]
        else:
            pass

        # check compression is compatible
        if compression == "gz":
            state = mapf["gzip"]
        elif compression == "bin":
            state = mapf["binary"]
        elif compression is None and not options:
            state = mapf["no_compress"]

        self.statement["matrix"] = state       

    def build_statement(self, infiles, outfile, threads=None,
                        memory=None, parallel=None):
        '''
        Build statement and execute from components
        '''

        statement = []
        exec_state = self.executable
        # calls to function add to the self.statement dictionary
        try:
            statement.append(self.statement["program"])
        except KeyError:
            raise AttributeError("Input files and format not detected")
            
        try:
            statement.append(self.statement["filters"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["tasks"])
        except KeyError:
            pass
        
        try:
            statement.append(self.statement["matrix"])
        except KeyError:
            pass

        if threads:
            statement.append(" --thread-num %i " % threads)
        else:
            pass

        # add output flag
        statement.append(" --out %s " % outfile)

        os.system(" ".join(statement))


class Plink2(GWASProgram):
    '''
    Run various Plink functions and analysis, including file processing, GRM
    calculation, PCA and other GWA tasks

    Require Plink v1.9 to be in the users PATH variable as ``plink2`` to
    distinguish it from Plink v1.07.
    '''

    def __init__(self, files, options=None,
                 settings=None, design=None):
        self.infiles = files
        self.options = options
        self.settings = settings
        self.design = design
        self.executable = "plink2"
        self.statement = {}
        self.filters = []

    def program_call(self, infiles, outfile):
        ''' build Plink call statement on infiles'''

        statement = []
        statement.append(self.executable)

        if infiles.name:
            inputs =self. _build_single_file_input(infiles,
                                                   infiles.file_format)
            statement.append(inputs)

        else:
            inputs = self._build_multiple_file_input(infiles,
                                                     infiles.file_format)
            statement.append(inputs)

        # check for the presence of an additional phenotypes file
        try:
            if infiles.pheno_file:
                statement.append(" --pheno %s --mpheno %s " % (infiles.pheno_file,
                                                               infiles.pheno))
            else:
                pass
        except AttributeError:
            pass

        self.statement["program"] = " ".join(statement)

    def hamming_matrix(self, shape, compression, options):
        '''
        Calculate genomic pair-wise distance matrix between
        individuals using Hamming distance across all variants
        '''

        # check shape is compatible
        if not shape:
            shape = "triangle"
        elif shape in ["square", "square0", "triangle"]:
            pass
        else:
            raise ValueError("matrix shape %s not recognised."
                             "Valid options are square, square0, "
                             "and triangle." % shape)

        # check compression is compatible
        if compression in ["gz", "bin", "bin4"]:
            pass
        else:
            raise ValueError("compression %s not recognised. Accepted "
                             "formats are gz, bin and bin4." % compression)

        if options:
            state = self._matrices(matrix_type="hamming", shape=shape,
                                   compression=compression, options=options)
        else:
            state = self._matrices(matrix_type="hamming", shape=shape,
                                   compression=compression)

        self.statement["matrix"] = state

    def ibs_matrix(self, shape, compression, options):
        '''
        Calculate genomic pair-wise similarity matrix between
        individuals using proportion of IBS alleles
       '''

        # check shape is compatible
        if shape in ["square", "square0", "triangle"]:
            pass
        else:
            raise ValueError("matrix shape %s not recognised."
                             "Valid options are square, square0, "
                             "and triangle." % shape)

        # check compression is compatible
        if compression in ["gz", "bin", "bin4"]:
            pass
        else:
            raise ValueError("compression %s not recognised. Accepted "
                             "formats are gz, bin and bin4." % compression)

        if options:
            state = self._matrices(matrix_type="ibs", shape=shape,
                                   compression=compression, options=options)
        else:
            state = self._matrices(matrix_type="ibs", shape=shape,
                                   compression=compression)

        self.statement["matrix"] = state

    def genome_matrix(self, shape, compression, options):
        '''
        Calculate genomic pair-wise distance matrix between
        individuals using 1 - proportion of IBS alleles
       '''

        # check shape is compatible
        if shape in ["square", "square0", "triangle"]:
            pass
        else:
            raise ValueError("matrix shape %s not recognised."
                             "Valid options are square, square0, "
                             "and triangle." % shape)

        # check compression is compatible
        if compression in ["gz", "bin", "bin4"]:
            pass
        else:
            raise ValueError("compression %s not recognised. Accepted "
                             "formats are gz, bin and bin4." % compression)

        if options:
            state = self._matrices(matrix_type="genomic", shape=shape,
                                   compression=compression, options=options)
        else:
            state = self._matrices(matrix_type="genomic", shape=shape,
                                   compression=compression)

        self.statement["matrix"] = state

    def genetic_relationship_matrix(self, shape, compression, metric,
                                    options=None):
        '''
        Calculate genomic pair-wise distance matrix between
        individuals using proportion of IBS alleles
        Requires the use of the Plink2 parallelisation to run with large
        cohorts of patients
        '''

        # check shape is compatible
        if shape in ["square", "square0", "triangle"]:
            pass
        else:
            raise ValueError("matrix shape %s not recognised."
                             "Valid options are square, square0, "
                             "and triangle." % shape)

        # check compression is compatible
        if compression in ["gz", "bin", "bin4"]:
            pass
        else:
            raise ValueError("compression %s not recognised. Accepted "
                             "formats are gz, bin and bin4." % compression)

        if metric in ["cov", "ibc2", "ibc3"]:
            state = self._matrices(matrix_type="grm", shape=shape,
                                   compression=compression, options=metric)
        else:
            E.info("%s metric not recognised.  Running with default Fhat1" % metric)
            state = self._matrices(matrix_type="grm", shape=shape,
                                   compression=compression)

        self.statement["matrix"] = state       

    def filter_genotypes(self, filter_type, filter_value):
        '''
        arguments supported by this function.

        * genotype_rate - exclude SNPs with a genotyping rate below this
          value. [float]        
        * min_allele_frequency - only include SNPs with cohort/case allele
          frequency above this threshold. [float]
        * max_allele_frequency - include all SNPs with a MAF equal to or below
          this value. [float]
        * exclude_snp - exclude this single variant
        * exclude_snps - text file list of variant IDs to exclude from analysis.
          [file]
        * chromosome - exclude all variants not on the specified chromosome(s).
          [str/list]
        * exclude_chromosome - exclude all variants on the specified
          chromosome(s). [str/list]
        * autosome - exclude all non-place and non-autosomal variants.
          [boolean]
        * pseudo_autosome - include the pseudo-autosomal region of chromosome
          X. [boolean]
        * ignore_indels - remove all indels/multi-character allele coding
          variants. [boolean]
        * snp_bp_range - (from, to) range in bp of variants to include in
          analysis. [tuple]

        '''

        if filter_type == "genotype_rate":
            self._construct_filters(genotype_rate=filter_value)
        elif filter_type == "hwe":
            self._construct_filters(hwe=filter_value)
        elif filter_type == "missingness":
            self._construct_filters(missingness=filter_value)
        elif filter_type == "min_allele_frequency":
            self._construct_filters(min_allele_frequency=filter_value)
        elif filter_type == "max_allele_frequency":
            self._construct_filters(max_allele_frequency=filter_value)
        elif filter_type == "exclude_snp":
            self._construct_filters(exclude_snp=filter_value)
        elif filter_type == "exclude":
            self._construct_filters(exclude=filter_value)
        elif filter_type == "extract":
            self._construct_filters(extract=filter_value)
        elif filter_type == "chromosome":
            self._construct_filters(chromosome=filter_value)
        elif filter_type  == "exclude_chromosome":
            self._constuct_filters(exclude_chromosome=filter_value)
        elif filter_type == "autosome":
            self._construct_filters(autosome=filter_value)
        elif filter_type == "pseudo_autosome":
            self._construct_filters(pseudo_autosome=filter_value)
        elif filter_type == "ignore_indels":
            self._construct_filters(ignore_indels=filter_value)
        elif filter_type == "snp_bp_range":
            self._construct_filters(snp_bp_range=filter_value)
        elif filter_type == "conditional_snp":
            self._construct_filters(conditional_snp=filter_value)

    def _build_multiple_file_input(self, infiles, file_format):
        '''
        internal function only.  Use it to construct
        the appropriate file input flags
        '''

        statement = None

        if file_format == "oxford":
            statement = " --gen %s --sample %s " % (infiles.gen_file,
                                                    infiles.sample_file)
        elif file_format == "oxford_binary":
            statement = " --bgen %s --sample %s " % (infiles.bgen_file,
                                                     infiles.sample_file)
        elif file_format == "plink":
            statement = " --ped %s --map %s " % (infiles.ped_file,
                                                 infiles.sample_file)
        elif file_format == "plink_binary":
            statement = " --bed %s --bim %s --fam %s " % (infiles.bed_file,
                                                          infiles.bim_file,
                                                          infiles.fam_file)
        elif file_format == "vcf":
            statement = " --vcf %s " % infiles.vcf_file
        elif file_format == "bcf":
            statement = " --bcf %s " % infiles.vcf_file
        else:
            raise AttributeError("file format is not defined.  Please "
                                 "define the input file formats when "
                                 "instantiating a FileGroup object")

        return statement

    def _build_single_file_input(self, infiles, file_format):
        '''internal function only. Use it to construct the
        file input flags with --file, --bfile or --data
        '''

        statement = None

        if file_format == "plink":
            statement = " --file %s " % infiles.name
        elif file_format == "plink_binary":
            statement = " --bfile %s " % infiles.name
        elif file_format == "oxford" or file_format == "oxford_binary":
            statement = " --data %s" % infiles.name
        elif file_format == "GRM_plink":
            statement = " --grm.bin  %s " %infiles.name
        else:
            raise AttributeError("file format is not defined or recognised."
                                 "Please define the input corectly when "
                                 "instantiating a FileGroup object")

        return statement

    def _construct_filters(self, **kwargs):
        '''
        Add filter to each plink run. [data type]

        The filters accepted are defined below.  These are input as keyword
        arguments supported by this function.

        * genotype_rate - exclude SNPs with a genotyping rate below this
          value. [float]
        * missingness - exclude individuals with total genotype missingness
          above this value. [float]
        * hwe - p-value threshold for excluding SNPs deviating from
          Hardy-Weinberg expectations. [float]
        * min_allele_frequency - only include SNPs with cohort/case allele
          frequency above this threshold. [float]
        * max_allele_frequency - include all SNPs with a MAF equal to or below
          this value. [float]
        * mendelian_error - filter out samples/trios exceeding the error
          threshold. [float]
        * keep - keep individuals with matching individual and family IDs.
          [file]
        * remove - remove all individuals with matching individual and family
          IDs. [file]
        * quality_score_file - vcf file with variants and quality scores.  Use
          `qual_score_column` and `var_id_col` to specify which columns
          correspond to the quality score and variant ID columns.
          [file] <int> <int>
        * min_qual_score - alters the lower bound of the quality score
          threshold; default is 0.[int]
        * max_qual_score - sets an upper limit on the quality scores;
          default is Inf. [int]
        * allow_no_sex - prevents phenotypes set to missing if there is no
          gender information. [boolean]
        * enforce_sex - force phenotype missing when using --make-bed, --recode
          or --write-covar. [boolean]
        * subset_filter - filter on a particular subset.  Choices are: cases,
          controls, males, females, founders, nonfounders. [str]
        * extract - text file list of variant IDs to include in analysis,
          ignores all others. [file]
        * exclude - text file list of variant IDs to exclude from analysis.
          [file]
        * chromosome - exclude all variants not on the specified chromosome(s).
          [str/list]
        * exclude_chromosome - exclude all variants on the specified
          chromosome(s). [str/list]
        * autosome - exclude all non-place and non-autosomal variants.
          [boolean]
        * pseudo_autosome - include the pseudo-autosomal region of chromosome
          X. [boolean]
        * ignore_indels - remove all indels/multi-character allele coding
          variants. [boolean]
        * snp_bp_range - (from, to) range in bp of variants to include in
          analysis. [tuple]
        * specific_snp - only load the variant specified. [str]
        * exclude_snp - exclude this single variant
        * window_size - alters behaviour of `specific_snp` and `exclude_snp`
          to include/exclude SNPs within +/- half of this distance (kb) are
          also included. [float]
        * range_resolution - sets the resolution of the (from, to) range.
          Either bp, kb or mb. If set it will take the values from
          `snp_bp_range`. [str/int/float]
        * covariates_file - specify the covariates file with family and
          individual IDs in the first two columns.  Covariates are in the
          (n+2)th column. Only used in conjunction with `covariate_filter`.
          [file]
        * covariate_filter - covariate columns value to filter on.  Can be
          used with non-numeric values to filter out individuals with
          covariate =/= `covariate_filter` value. [str/int/float]
        * covariate_column - column number to apply filtering to if more
          than one covariate in the file. [int]
        '''

        statement = []

        # map of keyword arguments recognised to Plink2 filtering flags
        filter_map = {"genotype_rate": " --geno %s ",
                      "missingness": "--mind %s ",
                      "hwe": " --hwe %s ",
                      "min_allele_frequency": " --maf %s ",
                      "max_allele_frequency": " --max-maf %s ",
                      "mendelian_error": " --me %s ",
                      "keep": " --keep %s ",
                      "remove": " --remove %s ",
                      "quality_score_file": " --qual-scores %s ",
                      "qual_score_column": " %s ",
                      "var_id_col": " %s ",
                      "min_qual_score": " --qual-threshold %s ",
                      "max_qual_score": " --qual-max-threshold %s ",
                      "allow_no_sex": " --allow-no-sex ",
                      "enforce_sex": " --must-have-sex ",
                      "subset_filter": " --filter-%s ",
                      "extract": " --extract %s ",
                      "exclude": " --exclude %s ",
                      "chromosome": " --chr %s ",
                      "exclude_chromosome": " --not-chr %s ",
                      "autosome": " --autosome ",
                      "pseudo_autosome": " --autosome-xy ",
                      "ignore_indels": " --snps-only no-DI ",
                      "snp_id_range": " --from %s --to %s ",
                      "specific_snp": " --snp %s ",
                      "window_size": " --window %s ",
                      "exclude_snp": " --exclude-snp %s ",
                      "snp_bp_range": "--from-bp %s --to-bp %s ",
                      "covariates_file": " --filter %s ",
                      "covariate_filter": " %s ",
                      "covariate_column": " --mfilter %s ",
                      "missing_phenotype": " --prune ",
                      "conditional_snp": " --condition %s "}

        # compile all filters together, checking for dependencies.
        # use a mapping dictionary to extract the relevant flags and
        # combinations to use.
        filters = []
        filter_dict = {}
        for key, value in kwargs.iteritems():
            filter_dict[key] = value

        # need to check for covariates and qual scores - these
        # are more complex.  Deal with these first and remove
        # from dictionary once complete.
        try:
            assert filter_dict["quality_score_file"]
            assert filter_dict["qual_score_column"]
            assert filter_dict["var_id_col"]

            quals = []
            qual_file = filter_dict["quality_score_file"]
            score_col = filter_dict["qual_score_column"]
            id_col = filter_dict["var_id_col"]

            quals.append(filter_map["quality_score_file"] % qual_file)
            quals.append(filter_map["qual_score_column"] % score_col)
            quals.append(filter_map["var_id_col"] % id_col)

            # remove from dictionary
            filter_dict.pop("qual_score_column", None)
            filter_dict.pop("var_id_col", None)

            filters.append(" ".join(quals))

        except KeyError:
            pass

        try:
            assert filter_dict["covariates_file"]
            assert filter_dict["covariate_filter"]

            covars = []
            covar_file = filter_dict["covariates_file"]
            covar_val = filter_dict["covariate_filter"]
            covars.append(filter_map["covariates_file"] % covar_file)
            covars.append(filter_map["covariate_filter"] % covar_val)

            # check to filter on specific column numnber, default is 3rd file
            # column, i.e. (n+2)th column
            try:
                assert filter_dict["covariate_column"]
                covar_col = filter_dict["covariate_column"]
                covars.append(filter_map["covariate_column"] % covar_col)
                filter_dict.pop("covariate_column", None)
            except KeyError:
                pass

            # remove from dictionary
            filter_dict.pop("covariates_file", None)
            filter_dict.pop("covariate_filter", None)

            filters.append(" ".join(covars))

        except KeyError:
            pass

        # range_resolution and snp_bp_range are used together
        try:
            assert filter_dict["snp_bp_range"]
            flags = filter_map["snp_bp_range"]
            from_pos = filter_dict["snp_bp_range"].split(",")[0]
            to_pos = filter_dict["snp_bp_range"].split(",")[1]
            filters.append(flags % (from_pos, to_pos))

            # remove so they are not duplicated - source of bugs
            filter_dict.pop("snp_bp_range", None)

        except KeyError:
            pass

        for each in filter_dict.keys():
            try:
                assert filter_map[each]
                # check for data type <- behaviour is type dependent
                if type(filter_dict[each]) == 'bool':
                    filters.append(filter_map[each])
                # handle multiple arguments in string format
                elif len(filter_dict[each].split(",")) > 1:
                    vals = tuple(filter_dict[each].split(","))
                    filters.append(filter_map[each] % vals)
                else:
                    filter_val = filter_dict[each]
                    filters.append(filter_map[each] % filter_val)

            except KeyError:
                E.warn("%s filter not recognised, please see "
                       "documentation for allowed filters" % each)
                pass
                
        self.filters.append(" ".join(filters))
        self.statement["filters"] = " ".join(self.filters)

    def _run_tasks(self, parameter=None, **kwargs):
        '''
        Plink2 is capable of much more than just running basic association
        analyses.

        These include file processing, reformating, filtering, data summaries,
        PCA, clustering, GRM calculation (slow and memory intense), etc.

        multiple tasks can be added by separate calls to this function.
        For instance, adding phenotype and gender information using the
        update_samples task whilst change the file format.
        

        Tasks
        -----

        * change_format - convert from input format to an alternative format
          after applying filters.
        * change_missing_values - alters the genotype or phenotype missing
          value into the value supplied.
        * update_variants - use this to fill in missing variant IDs, useful
          for data from exome or whole-genome sequencing that have
          non-standard IDs.
        * update_samples - update phenotype and sample information
        * flip_strands - flip the strand for alleles, swaps A for T and
          C for G.
        * flip_scan - use the LD-based scan to check SNPs have not had
          incorrect strand assignment. Particularly useful if cases and
          controls were genotyped separately, or the cohort was genotyped
          in different batches.
        * sort - sort files by individual and/or family IDs
        * merge - merge new filesets with reference fileset.
        * merge_mode - handling of missing values and overwriting values
        * find_duplicates - find and output duplicate variants based on bp position,
          or variant ID.  Useful to output for the --exclude filtering flag.
        '''

        statement = []

        # set up a dictionary of recognised tasks with key word argument
        # values as further dictionaries. Use the parameter argument
        # to pass arguments by value to string formatting

        task_map = {'change_format': {"plink_binary": " --make-bed ",
                                      "plink": " --recode ",
                                      "oxford": " --recode oxford ",
                                      "oxford_binary": " --recode oxford gen-gz "},
                    "change_missing_values": {"genotype": " --missing-genotype %s ",
                                              "phenotype": " --missing-phenotype %s "},
                    "update_variants": {"variant_ids": " --set-missing-var-ids %s ",
                                        "missing_id": " --mising-var-code %s ",
                                        "chromosome": " --update-chr %s ",
                                        "centimorgan": " --update-cm %s ",
                                        "name": " --update-name %s ",
                                        "alleles": " --update-alleles  %s ",
                                        "map": " --update-map %s "},
                    "update_samples": {"sample_ids": " --update-ids %s ",
                                       "parents": " --update-parents %s ",
                                       "gender": " --update-sex %s %s "},
                    "flip_strands": {"all_samples": " --flip %s ",
                                     "subset": " --flip-subset %s "},
                    "flip_scan": {"default": " --flip-scan verbose ",
                                  "window": "--flip-scan --flip-scan-window %s ",
                                  "kb": " --flip-scan  --flip-scan-window-kb %s ",
                                  "threshold": " --flip-scan  --flip-scan-threshold %s "},
                    "sort": {"none": " --indiv-sort %s ",
                             "natural": " --indiv-sort %s ",
                             "ascii": " --indiv-sort %s ",
                             "file": " --indiv-sort %s "},
                    "merge": {"plink": " --merge %s ",
                              "binary_plink": " --bmerge %s "},
                    "merge_mode": {"default": " --merge-mode 1 ",
                                   "orginal_missing": " --merge-mode 2 ",
                                   "new_nonmissing": " --merge-mode 3 ",
                                   "no_overwrite": " --merge-mode 4 ",
                                   "force": " --merge-mode 5 ",
                                   "report_all": " --merge-mode 6 ",
                                   "report_nonmissing": " --merge-mode 7"},
                    "find_duplicates": {"same_ref": " --list-duplicate-vars require-same-ref ",
                                        "id_match": " --list-duplicate-vars ids-only ",
                                        "suppress_first": " --list-duplicate-vars suppress-first"},
                    "pca": " --pca %s "}


        for task, value in kwargs.iteritems():
            # check for PCA first as it is not nested in task_map
            if task == "pca":
                try:
                    state = task_map[task] % value
                    statement.append(state)
                except TypeError:
                    statement.append(task_map[task])
                statement.append
            elif task != "parameter":
                try:
                    # sub_task is a nested dictionary
                    sub_task = task_map[task]
                    try:
                        assert sub_task[value]
                        try:
                            # gender has two string formats
                            if value == "gender":
                                gcol = 1
                                statement.append(sub_task[value] % (parameter,
                                                                    gcol))
                            else:
                                # some tasks do not contain task values for the
                                # parameter argument - catch these with the TypeError
                                # exception
                                statement.append(sub_task[value] % parameter)
                            # the default for parameter is None, check this is appropriate
                            if not parameter:
                                E.warn("Parameter value is set to NoneType. "
                                       "Please check this is an appropriate value "
                                       "to pass for this task")
                            else:
                                pass
                        except TypeError:
                            statement.append(sub_task[value])
                    except KeyError:
                        raise KeyError("No sub task found, see docs for details of "
                                       "recognised tasks")                    
                except KeyError:
                    raise KeyError("Task not recognised, see docs for details of "
                                   "recognised tasks")
            else:
                pass
        # handle multiple tasks for a single run
        try:
            curr_tasks = self.statement["tasks"]
            new_tasks = " ".join(statement)
            self.statement["tasks"] = " ".join([curr_tasks, new_tasks])
        except KeyError:
            self.statement["tasks"] = " ".join(statement)


    def _output_statistics(self, **kwargs):
        '''
        Summary statistics are written to specific files dictated by the
        type of statistic

        Statistics
        ----------
        * allele_frequency - writes out MAF to `plink`.frq, this can be
          modified with specific keywords.
        * missing_data - generates a report of data missingness, can be subset
          into within family and/or cluster reports
        * hardy_weinberg - calculates all HWE p-values using exact test
          statistics. For case/control studies reports are written for case,
          controls and combined.
        * mendel_errors - generates a Mendelian error report across all trios.
          There are 10 different codes responding to different Mendelian error
          scenarios.
        * inbreeding - calculate observed and expected homozygosity across
          individuals and F statistics.  If the sample size is small then a
          file of MAFs is required. Inbreeding coefficients can also be
          reported on request using inbreeding_coef.
        * gender_checker - checks gender assignment against X chromosome
          genotypes. Gender values can also be imputed based on genotype
          information using gender_impute.
        * wrights_fst - calculate Wright's Fst statistic given a set of
          subpopulations for each autosomal diploid variant.  Used in
          conjunction with the --within flag.
        '''

        stats_map = {"allele_frequency": " --freq %s ",
                     "missing_data": " --missing %s ",
                     "hardy_weinberg": " --hardy midp ",
                     "mendel_errors": " --mendel %s ",
                     "inbreeding": " --het %s ",
                     "inbreeding_coef": " --ibc ",
                     "gender_checker": " --check-sex ",
                     "gender_impute": " --impute-sex ",
                     "wrights_fst": " --fst --within %s "}

        statement = []
        for key, value in kwargs.iteritems():
            if value:
                try:
                    assert stats_map[key]
                    statement.append(stats_map[key] % value)
                except KeyError:
                    raise KeyError("statistic not recognised.  Please "
                                   "consult the documentation for allowed "
                                   "options.")
            else:
                try:
                    assert stats_map[key]
                    flag = stats_map[key].rstrip("%s ")
                    statement.append(flag)
                except KeyError:
                    raise KeyError("statistic not recognised.  Please "
                                   "consult the documentation for allowed "
                                   "options.")
        self.statement["stats"] = " ".join(statement)

    def run_association(self, association=None, model=None,
                        run_options=None,
                        permutation=False, n_perms=None,
                        random_seed=None, permutation_options=None,
                        covariates_file=None, covariates=None):
        '''
        Construct a statement for a plink2 association analysis.

        QC filters are constructed from input during instantiation.

        run options include redirecting logging output, using parallelisation,
        defining number of threads to use, etc

        The default association uses the --assoc flag.  Plink will check
        phenotype coding, if it is not case/control it assumes
        it is a continuous trait and uses linear regression.

        Alternative regression models that include covariates can be used,
        i.e. logistic and linear regression.

        key
        ***
        {CC} - applies to case/control analysis only
        {quant} - applies to quantitative trait only
        {CC/quant} - applies to both

        run_options
        -----------
        ``--assoc``:
            * `fisher | fisher-midp` - uses Fisher's exact test to calculate
            association p-values or applies Lancaster's mid-p adjustment. {CC}
            * `counts` - causes --assoc to report allele counts instead of
            frequencies. {CC}
            * `set-test` - implements and tests the significance of variant
            sets.  See documentation below. {CC/quant}
            * `qt-means` - generates a .qassoc.means file reporting trait means
            and standard deviations by genotype. {quant}
            * `lin` - reports the Lin et al (2006) statistic to be reported. If
            multiple testing adjustments and/or permutation is also used, they
            will be based on this statistic. {quant}

        ``--model``:
            * `fisher | fisher-midp | trend-only` - uses Fisher's exact test
            to calculate association p-values or applies Lancaster's mid-p
            adjustment. trend-only forces only a trend test to be performed.
            {CC}
            * `dom | rec | gen | trend` - use the specified test as the basis
            for the model permutation.  If none are defined the result with the
            smallest p-value is reported. {CC}
            * --cell - sets the minimum number of observations per cell in the
            2x3 contingency table.  The default is 0 with the Fisher and
            Fiser-midp test, otherwise 5. {CC}

        ``--linear/logistic``:
            * `set-test` - implements and tests the significance of variant
            sets.  See documentation below. {CC/quant}
            * `hide-covar` - removes the covariate specific sections from the
            results output. {CC/quant
            * `sex | no-x-sex` - `sex` adds sex as covariate to all models,
            whislt `no-x-sex` does not include gender into X-chromosome SNP
            models. {CC/quant}
            * `interaction` - adds in genotype X covariate interaction terms
            into the model. Can only be used with permutation is ``--tests``
            is also specified. {CC/quant}
            * `beta` - reports the beta coefficients instead of the OR in a
            logistic model. {CC}
            * `standard-beta` - standardizes the phenotype and all predictor
            variables to zero mean and unit variance prior to regressionm
            (separate for each variant analysed). {quant}
            * `intercept` - includes the intercept in the output results.
            {quant}

        model
        -----
        * `recessive` - `recessive` specifies the model assuming the A1 allele
          as recessive. {CC/quant}
        * `dominant` - `dominant` specifies the model assuming the A1 allele is
          dominant. {CC/quant}
        * `genotype` - `genotype` adds an additive effect/dominance deviation
          2df joint test with two genotype variables in the test (coded 0/1/2
          and 0/1/0). {CC/quant}
        * `trend` - forces a trend test to be performed. {CC/quant}
        * `hethom` - `hethom` uses 0/0/1 and 0/1/0 instead of the genotype
          coding. With permutation it will be based on the joint test instead
          of just the additive effects.  This can be overriden using the
          `--tests` flag. {CC/quant}
        * `no-snp` - `no-snp` defines a regression of phenotype on covariates
          without reference to genotype data, except where `--conditon{-list}`
          is specified. If used with permuation, test results will be reported
          for every covariate. {CC/quant}

        permutation
        -----------
        If permutation is True, run an adaptive Monte Carlo permutation test.
        If n_perms is set, this will run a max(T) permutation test with the n
        replications. A random seed will need to be provided.

        * `perm-count` - this alters the permutation output report to include
          counts instead of frequencies

        covariates
        ----------
        These should be provided in a separate file.  Specifying which
        covariates to include can be done as either a comma-separated list
        of covariate names or numbers. These numbers will correspond to the
        (n+2)th covariate file column as per the plink documentation.
        '''

        # model map maps common option effects onto specific syntax
        model_map = {"--logistic": {"recessive": "recssive",
                                    "dominant": "dominant",
                                    "genotype": "genotypic"},
                     "--linear": {"recessive": "recssive",
                                  "dominant": "dominant",
                                  "genotype": "genotypic"},
                     "--model": {"recessive": "rec",
                                 "dominant": "dom",
                                 "genotype": "gen"}}

        statement = []
        # construct analysis flags
        # add model, i.e. additive, recessive, dominant, etc.
        # see docstring for details.  Make sure correct modifier is used
        # with a mapping dictionary

        if association == "logistic":
            statement.append(" --logistic ")
            m_map = model_map["--logistic"]
            if model:
                statement.append(m_map[model])
            else:
                pass

        elif association == "linear":
            statement.append(" --linear ")
            m_map = model_map["--linear"]
            if model:
                statement.append(m_map[model])
            else:
                pass

        elif association == "model":
            statement.append(" --model ")
            m_map = model_map["--model"]
            statement.append(m_map[model])
        else:
            statement.append(" --assoc ")

        # add in run options.  These need to be in their correct
        # format already
        if run_options:
            modifiers = " ".join(run_options)
            statement.append(modifiers)
        else:
            pass

        # permutation should have a random seed set by the user.  Allow
        # this to set it's own seed if one not provided, but report it in
        # the log file
        if permutation:
            try:
                assert random_seed
            except AssertionError:
                rand_seed = random.randint(0, 100000000)
                E.warn("No seed is provided for the permutation test. "
                       "Setting seed to %s. Record this for future "
                       "replicability" % random_seed)
            if n_perms:
                statement.append(" mperm=%i --seed %s " % (n_perms,
                                                              random_seed))
            else:
                statement.append(" perm --seed %s " % (random_seed))
        else:
            pass

        # if using linear or logistic, covariates can be added into the model
        # to adjust for their effects - assumes fixed effects of covariates
        # mixed models are not yet implemented in Plink2.

        if covariates:
            if type(covariates[0]) == str:
                covar_names = ",".join([c for c in covariates])
                m_covar = " --covar-name %s " % covar_names

            elif type(covariates[0]) == int:
                covar_nums = ",".join([ci for ci in covariates])
                m_covar = " --covar-number %s " % covar_nums
            else:
                # if none are specified then don't adjust the model for any
                # and log a warning
                E.warn("Covariate header or numbers are not recognised."
                       "No covariates will be included in the model.  Please"
                       "specifiy them exactly")
                covariates = None
                covariates_file = None
        else:
            pass

        if covariates and covariates_file:
            statement.append(" --covar %s %s " % m_covar)
        elif covariates and not covaries_file:
            E.warn("No covariate file specified.  None included in model.")
        elif covariates_file and not covariates:
            E.warn("No covariates specified to include in the model."
                   "None included")
        else:
            pass

        self.statement["assoc"] = " ".join(statement)


    def PCA(self, n_pcs="20"):
        '''
        Perform PCA analysis on previosly generated GRM, output the number n
        principal componets, default = 20
        '''

        self._run_tasks(pca=n_pcs)


    def _dimension_reduction(self, **kwargs):
        '''
        Use PCA to perform dimensionality reduction on
        input samples.  A PCA  can be calculated using
        a subset of samples which can then be projected on
        to other samples.
        '''

        # FINISH ME!!!!

    def _matrices(self, matrix_type, shape="triangle", compression=None, options=None):
        '''
        Calculate a number of different distance matrices:
        realised genetic relationship matrix
        relationship covariance matrix
        identity by descent/state matrix
        hamming distance matrix

        * matrix_type - matrix to compute.  Can be either IBS, 1 - IBS,
          Hamming, GRM
          
        '''

        statement = []
        if matrix_type == "hamming":
            flag = " --distance "
        elif matrix_type == "ibs":
            flag = " --distance ibs "
        elif matrix_type == "genomic":
            flag = " --distance 1-ibs "
        elif matrix_type == "grm":
            flag = " --make-grm-bin "

        if options:
            statement.append(" ".join([flag, shape, compression, options]))
        elif matrix_type == "grm":
            statement.append(flag)
        else:
            statement.append(" ".join([flag, shape, compression]))

        return " ".join(statement)

    def _qc_methods(self, parameter=None, **kwargs):
        ''''
        Perform QC on genotyping data, SNP-wise and sample-wise.
        All arguments are passed as key word arguments, except
        cases detailed in `Parameters` where they are passed with
        the ``parameter`` argument.
        
        Methods
        -------
        * ld_prune - generate a list of SNPs in linkage equilibrium by
          pruning SNPs on either an LD statistic threshold, i.e. r^2,
          or use a variance inflation factor (VIF) threshold
        * heterozygosity - calculate average heterozygosity from each
          individual across a set of SNPs, threshold on individuals
          with deviation from expected proportions
        * ibd - calculate the genetic relationship of individuals to
          infer relatedness between individuals, threshold on given
          degree of relatedness, e.g. IBD > 0.0625, 3rd cousins
        * genetic_gender - estimate the gender of an individual
          from the X chromosome genotypes - correlate with reported
          gender and output discrepancies
        * ethnicity_pca - perform PCA using a subset of independent
          SNPs to infer genetic ancestry.  Compare and contrast this
          to individuals reported ancestry.  Report discrepancies
          and individuals  greater than a threshold distance away
          from a reference population.
        * homozygosity - identifies sets of runs of homozygosity
          within individuals.  These may be indicative of inbreeding,
          systematic genotyping errors or regions under selection.

        Parameters
        ----------
        Method parameters can also be passed through this function
        as keyword=value pairs.
        * ld_prune:
          `kb` - this modifier changes the window resolution to kb
          rather than bp.
          `r2` - the r^2 threshold above which SNPs are to be removed
          `vif` - the VIF threshold over which SNPs will be removed
          `window` - window size to calculate pair-wise LD over
          `step` - step size to advance window by
        '''

        qc_dict = {"ld_prune": {"R2": " --indep-pairwise %s %s %s ",
                                "VIF": " --indep %s %s %s "},
                   "heterozygosity": {"gz": " --het gz",
                                      "raw": " --het "},
                   "ibd": {"relatives": " --genome gz rel-check ",
                           "full": " --genome gz full ",
                           "norm": " --genome gz "},
                   "genetic_gender": "none",
                   "ethnicity_pca": "none",
                   "homozygosity": {"min_snp": " --homozyg-snp %s ",
                                    "min_kb": " --homozyg-kb %s ",
                                    "default": " --homozyg ",
                                    "density": " --homozyg-density ",
                                    "set_gap": " --homozyg-gap ",
                                    "snp_window": " --homozyg-window-snp %s ",
                                    "het_max": " --homozyg-het %s "}}

        task_dict = {}
        state = []

        # put everything in an accessible dictionary first
        for task, value in kwargs.iteritems():
            task_dict[task] = value

        # LD pruning can be passed multiple parameters,
        # handle this separately
        try:
            sub_task = task_dict["ld_prune"]
            ld_prune_task = qc_dict["ld_prune"]

            try:
                step = task_dict["step"]
            except KeyError:
                raise AttributeError("No step size found, please "
                                     "pass a step size to advance the "
                                     "window by")
            try:
                window = task_dict["window"]
                try:
                    task_dict["kb"]
                    window = "".join([window, "kb"])
                    task_dict.pop("kb", None)
                except KeyError:
                    pass

            except KeyError:                        
                raise AttributeError("No window size found.  Please input "
                                     "a window size to prune over")
            try:
                threshold = task_dict["threshold"]
            except KeyError:
                raise AttributeError("No threshold value, please input "
                                     "a value to LD prune SNPs on")

            # add in the kb if it is passed as an argument
            state.append(ld_prune_task[sub_task] % (window, step, threshold))

            task_dict.pop("threshold", None)
            task_dict.pop("ld_prune", None)
            task_dict.pop("window", None)
            task_dict.pop("step", None)

        except KeyError:
            pass

        for task,value in task_dict.iteritems():
            try:
                sub_task = qc_dict[task]
                try:
                    state.append(sub_task[value] % parameter)
                except ValueError:
                    state.append(sub_task[value])
            except KeyError:
                raise AttributeError("Task not found, please see "
                                     "documentation for available features")

        self.statement["QC"] = " ".join(state)

    def build_statement(self, infiles, outfile, threads=None,
                        memory="60G", parallel=None):
        '''
        Build statement and execute from components
        '''

        statement = []
        exec_state = self.executable
        # calls to function add to the self.statement dictionary
        try:
            statement.append(self.statement["program"])
        except KeyError:
            raise AttributeError("Input files and format not detected")
            
        try:
            statement.append(self.statement["QC"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["filters"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["tasks"])
        except KeyError:
            pass
        
        try:
            statement.append(self.statement["stats"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["assoc"])
        except KeyError:
            pass

        try:
            statement.append(self.statement["matrix"])
        except KeyError:
            pass

        if threads:
            statement.append(" --threads %i " % threads)
        else:
            pass

        if memory != "60G":
            memory = int(memory.strip("G")) * 1000
            statement.append(" --memory %i " % memory)
        else:
            statement.append(" --memory 60000 ")

        # add output flag
        # outfile needs to be complete path for Plink to save
        # results properly - check if it starts with '/',
        # if so is already a full path
        if os.path.isabs(outfile):
            statement.append(" --out %s " % outfile)
        else:
            outpath = "/".join([os.getcwd(), outfile])
            statement.append(" --out %s " % outpath)

        # parallelisation only really applies to GRM calculation
        # at the moment
        if parallel:
            statements = []
            cat_state = [" cat "]

            for i in range(1, parallel+1):
                p_state = statement[:] # copy list, assigning just makes a pointer
                p_state.append(" --parallel %i %i " % (i, parallel))
                statements.append(" ".join(p_state))
                cat_state.append(" %s.grm.bin.%i " % (outpath, i))
            cat_state.append(" > %s.grm.N.bin " % outpath)
            statements.append(" ".join(cat_state))
            os.system(";".join(statements))
        else:
            os.system(" ".join(statement))
        


class GWASResults(object):
    '''
    A class for handling the results from a GWA, used for plotting
    and post-analysis QC
    '''

    def __init__(self, assoc_file):
        self.infile = assoc_file
        # results is a pandas dataframe to operate on
        self.results = self.get_results(assoc_file)

    def get_results(self, association_file):
        '''
        Parse a GWA results file and assign the table
        as an attribute.
        '''

        # use Pandas for now - try something different later
        # SQLite DB maybe?
        # inconsistent number of white spaces between
        # fields means Pandas parsing breaks down
        # fields need to be the correct data type,
        # i.e. BP = int, P = float, SNP = str, etc

        l_count = 0
        with open(association_file, "r") as ifile:
            for line in ifile:
                parsed = line.split(" ")
                # remove multiple blank spaces
                for i in range(parsed.count('')):
                    parsed.remove('')
                # get rid of the newline
                parsed.remove("\n")
                if l_count == 0:
                    header = [iy for ix, iy in enumerate(parsed)]
                    head_idx = [ix for ix, iy in enumerate(parsed)]
                    map_dict = dict(zip(head_idx, header))
                    res_dict = dict(zip(header, [[] for each in header]))
                    l_count += 1
                else:
                    col_idx = [lx for lx, ly in enumerate(parsed)]
                    col = [ly for lx, ly in enumerate(parsed)]
                    for i in col_idx:
                        res_dict[map_dict[i]].append(col[i])
                    l_count += 1

        # substract one from the index for the header column
        df_idx = range(l_count-1)

        results_frame = pd.DataFrame(res_dict, index=df_idx)
        results_frame["BP"] = [int(bx) for bx in results_frame["BP"]]
        results_frame["P"] = [np.float64(fx) for fx in results_frame["P"]]
        try:
            results_frame["STAT"] = [np.float64(sx) for sx in results_frame["STAT"]]
        except KeyError:
            results_frame["CHISQ"] = [np.float64(sx) for sx in results_frame["CHISQ"]]
        try:
            results_frame["F_U"] = [np.float64(ux) for ux in results_frame["F_U"]]
        except KeyError:
            pass

        try:
            results_frame["F_A"] = [np.float64(ax) for ax in results_frame["F_A"]]
        except KeyError:
            pass

        try:
            results_frame["OR"] = [np.float64(ox) for ox in results_frame["OR"]]
        except KeyError:
            pass

        return results_frame

    def plotManhattan(self, save_path, resolution="chromosome"):
        '''
        Generate a basic manhattan plot of the association results
        Just deal with chromosome-by-chromosome for now.
        '''

        # use the python ggplot plotting package
        # need to calculate -log10P values separately
        self.results["log10P"] = np.log(self.results["P"])

        # manplot = ggplot(self.results, aes(x="BP",
        #                                    y="-log10P")) + \
        #     geom_point() + xlab("Chromosome position (bp)") + \
        #     ylab("-log10 P-value")

        # or using rpy2
        py2ri.activate()
        R('''suppressPackageStartupMessages(library(ggplot2))''')
        r_df = py2ri.py2ri_pandasdataframe(self.results)
        R.assign("assoc.df", r_df)
        R('''p <- ggplot(assoc.df, aes(x=BP, y=-log10P)) + geom_point() + '''
          '''theme_bw() + labs(x="Chromosome position (bp)", '''
          '''y="-log10 P-value")''')
        R('''png("%s")''' % save_path)
        R('''print(p)''')
        R('''dev.off()''')

    def plotQQ(self, save_path, resolution="chromosome"):
        '''
        Generate a QQ-plot of expected vs. observed
        test statistics
        '''

        pass


##########################################################
# unbound methods that work on files and data structures #
##########################################################

def plotMapPhenotype(data, coords, coord_id_col, lat_col,
                     long_col, save_path, xvar, var_type,
                     xlabels=None, level=None):
    '''
    Generate a map of the UK, with phenotype data overlaid
    '''
    
    # merge co-ordinate data with phenotype data
    merged_df = pd.merge(left=coords, right=data, left_on="f.eid",
                         right_on=coord_id_col, how='inner')

    # pheno column and set level of categorical variable
    if xlabels and var_type == "categorical":
        # convert to string type as a categorical variable
        # drop NA observations from the merged data frame
        na_mask = pd.isnull(merged_df.loc[:, xvar])
        merged_df = merged_df[~na_mask]

        rvar = merged_df.loc[:, xvar].copy()
        nvar = pd.Series(np.nan_to_num(rvar), dtype=str)
        var = [v for v in set(nvar)]
        var.sort()
        # recode the variables according to the input labels
        xlabs = xlabels.split(",")
        print xlabs, var
        lbls = [str(xlabs[ix]) for ix in range(len(var))]
        for xv in range(len(var)):
            nvar[nvar == var[xv]] = lbls[xv]
        merged_df.loc[:, "cat_var"] = nvar

    else:
        pass

    if level:
        lvar = merged_df.loc[:, "cat_var"].copy()
        mask = lvar.isin([level])
        lvar[mask] = 1
        lvar[~mask] = 0
        lvar = lvar.fillna(0)
        merged_df.loc[:, "dichot_var"] = lvar
    else:
        pass

    # push the df into the R env
    py2ri.activate()
    r_df = py2ri.py2ri_pandasdataframe(merged_df)
    R.assign("pheno.df", r_df)

    # setup the map and plot the points
    R('''suppressPackageStartupMessages(library(maps))''')
    R('''suppressPackageStartupMessages(library(mapdata))''')

    R('''uk_map <- map("worldHires", c("UK", "Isle of Wight",'''
      '''"Ireland", "Isle of Man", "Wales:Anglesey"), '''
      '''xlim=c(-11, 3), ylim=c(50, 60.9), plot=F)''')
    # colour by reference, or a colour for each discrete value
    if level:
        R('''red <- rep("#FF0000", '''
          '''times=length(pheno.df$dichot_var[pheno.df$dichot_var == 1]))''')
        R('''black <- rep("#000000", '''
          '''times=length(pheno.df$dichot_var[pheno.df$dichot_var == 0]))''')

        R('''png("%(save_path)s", width=540, height=540, res=90)''' % locals())
        R('''map(uk_map)''')

        R('''points((pheno.df[,"%(long_col)s"])[pheno.df$dichot_var == 0], '''
          '''(pheno.df[,"%(lat_col)s"])[pheno.df$dichot_var == 0], pch=".", col=black)''' % locals())

        R('''points((pheno.df[,"%(long_col)s"])[pheno.df$dichot_var == 1], '''
          '''(pheno.df[,"%(lat_col)s"])[pheno.df$dichot_var == 1], pch=".", col=red)''' % locals())

        R('''legend('topleft', legend=c("not-%(level)s", "%(level)s"),'''
          '''fill=c("#000000", "#FF0000"))''' % locals())
        R('''dev.off()''')
    else:
        R('''print(pheno.df$cat_var[pheno.df$cat_var == "nan"])''')
        R('''png("%(save_path)s", width=540, height=540, res=90)''' % locals())
        R('''map(uk_map)''')

        R('''points(pheno.df[,"%(long_col)s"], pheno.df[,"%(lat_col)s"], pch=".", '''
          '''col=factor(pheno.df$cat_var))''' % locals())

        R('''legend('topleft', legend=unique(pheno.df$cat_var),'''
          '''fill=unique(pheno.df$cat_var))''' % locals())
        R('''dev.off()''')
        
    
def plotPhenotype(data, plot_type, x, y=None, group=None,
                  save_path=None, labels=None, xlabels=None,
                  ylabels=None, glabels=None, var_type="continuous"):
    '''
    Generate plots of phenotypes using ggplot
    '''

    # change data format if necessary and convert nan/NA to missing
    if not y and var_type == "categorical":
        var = np.nan_to_num(data.loc[:, x].copy())
        data.loc[:, x] = pd.Series(var, dtype=str)
        if group:
            gvar = np.nan_to_num(data.loc[:, group].copy())
            data.loc[:, group] = pd.Series(gvar, dtype=str)
        else:
            pass

    elif not y and var_type == "integer":
        var = np.nan_to_num(data.loc[:, x].copy())
        data.loc[:, x] = pd.Series(var, dtype=np.int64)
        if group:
            gvar = np.nan_to_num(data.loc[:, group].copy())
            data.loc[:, group] = pd.Series(gvar, dtype=str)
        else:
            pass


    elif not y and var_type == "continuous":
        var = data.loc[:, x].copy()
        data.loc[:, x] = pd.Series(var, dtype=np.float64)
        if group:
            gvar = np.nan_to_num(data.loc[:, group].copy())
            data.loc[:, group] = pd.Series(gvar, dtype=str)
        else:
            pass

    elif y and var_type == "categorical":
        xvar = np.nan_to_num(data.loc[:, x].copy())
        yvar = np.nan_to_num(data.loc[:, y].copy())

        data.loc[:, x] = pd.Series(xvar, dtype=str)
        data.loc[:, y] = pd.Series(yvar, dtype=str)

        if group:
            gvar = np.nan_to_num(data.loc[:, group].copy())
            data.loc[:, group] = pd.Series(gvar, dtype=str)
        else:
            pass


    elif y and var_type == "integer":
        xvar = np.nan_to_num(data.loc[:, x].copy())
        yvar = np.nan_to_num(data.loc[:, y].copy())

        data.loc[:, x] = pd.Series(xvar, dtype=np.int64)
        data.loc[:, y] = pd.Series(yvar, dtype=np.int64)

        if group:
            gvar = np.nan_to_num(data.loc[:, group].copy())
            data.loc[:, group] = pd.Series(gvar, dtype=str)
        else:
            pass


    elif y and var_type == "continuous":
        # NAs and NaNs should be handled by ggplot
        xvar = data.loc[:, x].copy()
        yvar = data.loc[:, y].copy()

        data.loc[:, x] = pd.Series(xvar, dtype=np.float64)
        data.loc[:, y] = pd.Series(yvar, dtype=np.float64)

        if group:
            gvar = np.nan_to_num(data.loc[:, group].copy())
            data.loc[:, group] = pd.Series(gvar, dtype=str)
        else:
            pass


    R('''suppressPackageStartupMessages(library(ggplot2))''')
    # put the pandas dataframe in to R with rpy2
    py2ri.activate()

    r_df = py2ri.py2ri_pandasdataframe(data)
    R.assign("data_f", r_df)

    # plotting parameters, including grouping variables and labels
    # axis labels
    try:
        labs = labels.split(",")
    except AttributeError:
        labs = []

    # if variable labels have been provided then assume they are
    # categorical/factor variables.
    # assume variable labels are input in the correct order
    if xlabels:
        try:
            unique_obs = len(set(data.loc[:, x]))
            xfact = len(xlabels.split(","))
            if xfact == unique_obs:
                R('''lvls <- unique(data_f[,"%(x)s"])''' % locals())
                lbls = ro.StrVector([ri for ri in xlabels.split(",")])
                R.assign("lbls", lbls)
                R('''lvls <- lvls[order(lvls, decreasing=F)]''')
                R('''data_f[,"%(x)s"] <- ordered(data_f[,"%(x)s"], '''
                  '''levels=lvls, labels=lbls)''' % locals())
            else:
                E.warn("the number of labels does not match the "
                       "number of unique observations, labels not "
                       "used.")
                pass
        except AttributeError:
            xlabels = None
            
    else:
        pass

    if glabels:
        unique_obs = len(set(data.loc[:, group]))
        gfact = len(glabels.split(","))
        if gfact == unique_obs:
            R('''lvls <- unique(data_f[, "%(group)s"])''' % locals())
            lbls = ro.StrVector([rg for rg in glabels.split(",")])
            R.assign("lbls", lbls)
            R('''lvls <- lvls[order(lvls, decreasing=F)]''')
            R('''data_f[,"%(group)s"] <- ordered(data_f[,"%(group)s"], '''
              '''levels=lvls, labels=lbls)''' % locals())
        else:
            E.warn("the number of labels does not match the "
                   "number of unique observations, labels not "
                   "used.")
            pass

    # start constructing the plot
    # if X and Y are specified, assume Y is a variable to colour
    # observations by, unless group is also set.
    # If Y and group then colour by group and split by Y
    if y:
        R('''p <- ggplot(aes(x=%s, y=%s), data=data_f)''' % (x, y))

        if plot_type == "histogram":
            if group:
                R('''p <- p + geom_histogram(aes(colour=%(group)s)) + '''
                  '''facet_grid(. ~ %(y)s)''' % locals())
            else:
                R('''p <- p + geom_histogram(aes(colour=%(y)s))''' % locals())

        elif plot_type == "barplot":
            if group:
                R('''p <- p + geom_bar(aes(colour=%(group)s)) + '''
                  '''facet_grid(. ~ %(y)s)''' % locals())
            else:
                R('''p <- p + geom_bar(aes(colour=%(y)s))''' % locals())

        elif plot_type == "density":
            if group:
                R('''p <- p + geom_density(aes(colour=%(group)s)) + '''
                  '''facet_grid(. ~ %(y)s)''' % locals())
            else:
                R('''p <- p + geom_density(aes(colour=%(y)s))''' % locals())

        elif plot_type == "boxplot":
            if group:
                R('''p <- p + geom_boxplot(group=%(group)s,'''
                  '''aes(x=factor(%(x)s), y=%(y)s, fill=%(group)s))''' % locals())
            else:
                R('''p <- p + geom_boxplot(aes(colour=%(x)s))''' % locals())

        elif plot_type == "scatter":
            if group:
                R('''p <- p + geom_point(size=1, aes(colour=%(group)s))''' % locals())
            else:
                R('''p <- p + geom_point(size=1)''')

        if len(labs) == 1:
            xlab = labs[0]
            R('''p <- p + labs(x="%s")''' % xlab)
        elif len(labs) == 2:
            xlab = labs[0]
            ylab = labs[1]
            R('''p <- p + labs(x="%(xlab)s", y="%(ylab)s")''' % locals())
        elif len(labs) == 3:
            xlab = labs[0]
            ylab = labs[1]
            title = labs[2]
            R('''p <- p + labs(x="%(xlab)s", y="%(ylab)s", '''
              '''title="%(title)s")''' % locals())
        elif len(labs) == 4:
            xlab = labs[0]
            ylab = labs[1]
            glab = labs[2]
            title = labs[3]
            R('''p <- p + labs(x="%(xlab)s", y="%(ylab)s",'''
              '''title="%(title)s")''' % locals())
            # need to add in guide/legend title

    else:
        R('''p <- ggplot(data=data_f)''')  

        if plot_type == "histogram":
            if group:
                R('''p <- p + geom_histogram(aes(%(x)s)) + '''
                  '''facet_grid(. ~ %(group)s)''' % locals())
            else:
                R('''p <- p + geom_histogram(aes(%s))''' % x)

        elif plot_type == "barplot":
            if group:
                R(''' p <- p + geom_bar(aes(%(x)s)) + '''
                  '''facet_grid(. ~ %(group)s)''')
            else:
                R('''p <- p + geom_bar(aes(%s))''' % x)

        elif plot_type == "density":
            if group:
                R('''p <- p + geom_density(aes(%(x)s)) + '''
                  '''facet_grid(. ~ %(group)s)''' % locals())
            else:
                R('''p <- p + geom_density(aes(%s))''' % x)

        elif plot_type == "boxplot":
            if group:
                R('''p <- p + geom_boxplot(aes(y=%(x)s, '''
                  '''x=factor(%(group)s)))''' % locals())
            else:
                raise AttributeError("Y or group variable is missing")

        if len(labs) == 1:
            xlab = labs[0]
            R('''p <- p + labs(x="%s")''' % xlab)
        elif len(labs) == 2:
            xlab = labs[0]
            title = labs[1]
            R('''p <- p + labs(x="%(xlab)s", '''
              '''title="%(title)s")''' % locals())
        elif len(labs) == 3:
            if group:
                xlab = labs[0]
                glab = labs[1]
                title = labs[2]
                R('''p <- p + labs(x="%(glab)s", y="%(xlab)s",'''
                  '''title="%(title)s")''' % locals())
            else:
                E.warn("too many labels provided, assume first is X, "
                       "and second is plot title")
                xlab = labs[0]
                title = labs[1]
                R('''p <- p + labs(x="%(xlab)s", '''
                  '''title="%(title)s")''' % locals())
            
    # the default theme is bw
    R('''p <- p + theme_bw()''')

    R('''png("%(save_path)s")''' % locals())
    R('''print(p)''')
    R('''dev.off()''')


def countByVariantAllele(ped_file, map_file):
    '''
    Count the number of individuals carrying the variant allele
    for each SNP.
    Count the number of occurences of each allele with the variant
    allele of each other SNP.

    Requires ped file genotyping to be in format A1(minor)=1, A2=2
    '''

    # parse the ped file - get the variant column headers from
    # the map file - no headers with these files
    # variant order in the map file matters, use an ordered dict
    variants = collections.OrderedDict()
    with open(map_file, "r") as mfile:
        for snp in mfile.readlines():
            attrs = snp.split("\t")
            snpid = attrs[1]
            variants[snpid] = {"chr": attrs[0],
                               "pos": attrs[-1].strip("\n")}

    variant_ids = variants.keys()
    # store genotype matrix as an array
    # rows and columns are variant IDs
    homA1 = np.zeros((len(variant_ids), len(variant_ids)),
                     dtype=np.int64)

    homA2 = np.zeros((len(variant_ids), len(variant_ids)),
                     dtype=np.int64)

    het  = np.zeros((len(variant_ids), len(variant_ids)),
                     dtype=np.int64)

    tcount = 0
    with open(ped_file, "r") as pfile:
        for indiv in pfile.readlines():
            indiv = indiv.strip("\n")
            indiv_split = indiv.split(" ")
            fid = indiv_split[0]
            iid = indiv_split[1]
            mid = indiv_split[2]
            pid = indiv_split[3]
            gender = indiv_split[4]
            phen = indiv_split[5]
            alleles = indiv_split[6:]
            genos = ["".join([alleles[i],
                              alleles[i+1]]) for i in range(0, len(alleles), 2)]
            tcount += 1
            # get genotype counts
            for i in range(len(genos)):
                # missing genotypes are coded '00' in plink format
                if genos[i] == "00":
                    pass
                elif genos[i] == "11":
                    homA1[i, i] += 1
                elif genos[i] == "12":
                    het[i, i] += 1
                else:
                    homA2[i, i] += 1
    allele_counts = ((2 * homA2) + het)/float(2 * tcount)
    mafs = 1 - allele_counts.diagonal()
    maf_df = pd.DataFrame(zip(variant_ids, mafs), columns=["SNP", "MAF"],
                          index=[x for x,y in enumerate(variant_ids)])
    maf_df["A2_HOMS"] = (2 * homA1).diagonal()
    maf_df["A2_HETS"] = het.diagonal()
    maf_df.index = maf_df["SNP"]
    maf_df.drop(["SNP"], axis=1, inplace=True)

    E.info("allele frequencies calculated over %i SNPs and "
           "%i individuals" % (len(genos), tcount))

    return maf_df


def calcPenetrance(ped_file, map_file, mafs=None,
                   subset=None):
    '''
    Calculate the proportion of times an allele is observed
    in the phenotype subset vs it's allele frequency.
    This is the penetrance of the allele
    i.e. if observed in 100% of affected individuals and 0%
    of controls, then penetrance is 100%
    Generates a table of penetrances for each variants/allele
    and a plot of MAF vs # cases carrying the allele

    Generates a heatmap of compound heterozygotes, and homozygotes
    with penetrances.

    Outputs a table of SNPs, homozygote and heterozygote counts
    among subset individuals and proportion of subset individual
    phenotype explained by homozygotes and heterozygotes

    Requires alleles are coded A1(minor)=1, A2=2
    '''
    # check subset is set, if not then throw an error
    # cannot calculate penetrance without a phenotype
    if not subset:
        raise ValueError("Cannot calculate penetrance of alleles "
                         "without a phenotype to subset in")
    else:
        pass

    # parse the ped file - get the variant column headers from
    # the map file - no headers with these files
    # variant order in the map file matters, use an ordered dict
    variants = collections.OrderedDict()
    with open(map_file, "r") as mfile:
        for snp in mfile.readlines():
            attrs = snp.split("\t")
            snpid = attrs[1]
            variants[snpid] = {"chr": attrs[0],
                               "pos": attrs[-1].strip("\n")}

    variant_ids = variants.keys()
    case_mat = np.zeros((len(variant_ids), len(variant_ids)),
                        dtype=np.float64)

    all_mat = np.zeros((len(variant_ids), len(variant_ids)),
                       dtype=np.float64)

    tcount = 0
    ncases = 0
    # missing phenotype individuals must be ignored, else
    # they will cause the number of individuals explained
    # to be underestimated
    with open(ped_file, "r") as pfile:
        for indiv in pfile.readlines():
            indiv = indiv.strip("\n")
            indiv_split = indiv.split(" ")
            fid = indiv_split[0]
            iid = indiv_split[1]
            mid = indiv_split[2]
            pid = indiv_split[3]
            gender = int(indiv_split[4])
            phen = int(indiv_split[5])
            if phen != -9:
                if subset == "cases":
                    select = phen
                elif subset == "gender":
                    select = gender
                else:
                    select = None
                alleles = indiv_split[6:]
                genos = ["".join([alleles[i],
                                  alleles[i+1]]) for i in range(0, len(alleles), 2)]
                tcount += 1

                het = np.zeros(len(genos), dtype=np.float64)
                hom = np.zeros(len(genos), dtype=np.float64)

                for i in range(len(genos)):
                    # missing values are coded '00' in plink format
                    # A2 homs are coded '11' in plink format
                    if genos[i] == "11":
                        hom[i] += 1
                    elif genos[i] == "12":
                        het[i] += 1
                    else:
                        pass

                hom_mat = np.outer(hom, hom)
                het_mat = np.outer(het, het)
                homs = hom_mat.diagonal()
                het_mat[np.diag_indices(len(genos))] = homs

                gen_mat = het_mat

                # separate matrix for subset
                # reference is always level 2 for plink files,
                # either cases or females
                if select == 2:
                    case_mat += gen_mat
                    all_mat += gen_mat
                    ncases += 1
                else:
                    all_mat += gen_mat
            else:
                pass

    E.info("alleles counted over %i SNPs "
           "and %i individuals, of which %i are "
           "in the %s subset" % (len(genos), tcount, ncases, subset))

    penetrance = np.divide(case_mat, all_mat)
    # round for the sake of aesthetics
    penetrance = np.round(penetrance, decimals=5)
    pen_df = pd.DataFrame(penetrance, columns=variant_ids,
                          index=variant_ids)
    pen_df = pen_df.fillna(0.0)

    case_df = pd.DataFrame(case_mat, columns=variant_ids,
                           index=variant_ids)
    all_df = pd.DataFrame(all_mat, columns=variant_ids,
                          index=variant_ids)
    # plot heatmap of penetrances as percentages
    indf = pen_df * 100
    py2ri.activate()

    # only plot penetrances > 0%
    r_pen = py2ri.py2ri_pandasdataframe(indf)
    r_cases = py2ri.py2ri_pandasdataframe(case_df)
    r_all = py2ri.py2ri_pandasdataframe(all_df)
    R.assign("pen.df", r_pen)
    R.assign("case.df", r_cases)
    R.assign("all.df", r_all)
    R('''suppressPackageStartupMessages(library(gplots))''')
    R('''suppressPackageStartupMessages(library(RColorBrewer))''')

    # penetrances
    E.info("plotting penetrance matrix")
    R('''hmcol <- colorRampPalette(brewer.pal(9, "BuGn"))(100)''')
    R('''rowpen <- pen.df[rowSums(pen.df) > 0,]''')
    R('''colpen <- rowpen[,colSums(rowpen) > 0]''')
    R('''png("%s/penetrance-matrix.png", width=720, height=720)''' % os.getcwd())
    R('''heatmap.2(as.matrix(colpen), trace="none", col=hmcol,'''
      '''dendrogram="none", Colv=colnames(colpen), key=FALSE, '''
      '''Rowv=rownames(colpen), margins=c(10,10), cellnote=round(colpen),'''
      '''notecol="white")''')
    R('''dev.off()''')

    E.info("plotting case counts matrix")
    R('''rowcase <- case.df[rowSums(case.df) > 0,]''')
    R('''colcase <- rowcase[,colSums(rowcase) > 0]''')
    R('''png("%s/cases-matrix.png", width=720, height=720)''' % os.getcwd())
    R('''heatmap.2(as.matrix(colcase), trace="none", col=rep("#F0F8FF", 100),'''
      '''dendrogram="none", Colv=colnames(colcase), key=FALSE, '''
      '''colsep=seq(1:length(colnames(colcase))), '''
      '''rowsep=seq(1:length(rownames(colcase))),'''
      '''Rowv=rownames(colcase), margins=c(10,10), cellnote=round(colcase),'''
      '''notecol="black")''')
    R('''dev.off()''')

    E.info("plotting all individuals matrix")
    R('''rowall <- all.df[rownames(colcase),]''')
    R('''colall <- rowall[,colnames(colcase)]''')
    R('''png("%s/all-matrix.png", width=720, height=720)''' % os.getcwd())
    R('''heatmap.2(as.matrix(colall), trace="none", col=rep("#F0F8FF", 100),'''
      '''dendrogram="none", Colv=colnames(colall), key=FALSE, '''
      '''colsep=seq(1:length(colnames(colall))), '''
      '''rowsep=seq(1:length(rownames(colall))), '''
      '''Rowv=rownames(colall), margins=c(10,10), cellnote=round(colall),'''
      '''notecol="black")''')
    R('''dev.off()''')

    # plot MAF vs homozygosity
    maf_df = pd.read_table(mafs, sep="\t", header=0, index_col=0)
    plot_df = pd.DataFrame(columns=["MAF"],
                           index=maf_df.index)
    plot_df["MAF"] = maf_df["MAF"]

    homs = case_mat.diagonal()
    hom_series = pd.Series({x: y for x, y in zip(variant_ids,
                                                 homs)})
    plot_df["explained_by_homozygotes"] = hom_series
    plot_df["SNP"] = plot_df.index
    plot_df.index = [ix for ix, iy in enumerate(plot_df.index)]
    plotPenetrances(plotting_df=plot_df)

    out_df = summaryPenetrance(maf_df=maf_df,
                               case_counts=case_mat,
                               variants=variant_ids,
                               n_cases=ncases,
                               n_total=tcount)
    return out_df, pen_df


def summaryPenetrance(maf_df, case_counts,
                      variants, n_cases, n_total):
    '''
    Summarise genotype counts and proportion of cases explained
    by the observed homozygotes and compound heterozygotes.
    This is a function of the total population size and
    population allele frequency - does this assume 100%
    penetrance of each allele?
    '''

    # homozygous individuals  are on the
    # diagonal of the case_counts array
    homozyg_cases = case_counts.diagonal()
    homozyg_series = pd.Series({x: y for x, y in zip(variants,
                                                     homozyg_cases)})

    # heterozygotes are on the off-diagonal elements
    # get all off diagonal elements by setting diagonals to zero
    # matrix is diagonal symmetric
    np.fill_diagonal(case_counts, 0)

    het_counts = np.sum(case_counts, axis=0)
    het_series = pd.Series({x: y for x, y in zip(variants,
                                                  het_counts)})
    out_df = pd.DataFrame(columns=["homozygote_cases",
                                   "heterozygote_cases"],
                          index=maf_df.index)
    out_df["MAF"] = maf_df["MAF"]
    out_df["homozygote_cases"] = np.round(homozyg_series, 1)
    out_df["expected_cases"] = np.round(((out_df["MAF"] ** 2) * n_total), 3)
    out_df["heterozygote_cases"] = het_series
    out_df["hom_prop_explained"] = np.round(homozyg_series/float(n_cases), 3)
    out_df["het_prop_explained"] = np.round(het_series/float(n_cases), 3)

    return out_df


def plotPenetrances(plotting_df):
    '''
    Plot the proportion of cases/phenotype explained by 
    individuals carrying allele vs. population allele frequency.
   
    Generate final output summary table (should be in separate function)
    '''

    # only need to plot variants with MAF >= 0.01
    low_frq = plotting_df["MAF"] < 0.01
    hi_df = plotting_df[~low_frq]
    
    # get into R and use ggplot for MAF vs homozygosity amongs cases
    r_plot = py2ri.py2ri_pandasdataframe(hi_df)
    R.assign("hom.df", r_plot)

    R('''suppressPackageStartupMessages(library(ggplot2))''')
    R('''png("%s/penetrance-plot.png", height=720, width=720)''' % os.getcwd())
    R('''pen_p <- ggplot(hom.df, aes(x=explained_by_homozygotes, y=MAF, colour=SNP)) + '''
      '''geom_point(size=4) + theme_bw() + '''
      '''geom_text(aes(label=explained_by_homozygotes),'''
      '''colour="black",vjust=0.5, hjust=0.5) + '''
      '''labs(x="Number of Red haired homozygotes", y="MAF") + '''
      '''theme(axis.title=element_text(size=10, colour="black"))''')
    R('''print(pen_p)''')
    R('''dev.off()''')


def findDuplicateVariants(bim_file, take_last=False):
    '''
    identify variants with duplicate position and reference
    alleles
    '''

    # count the number of lines first to get
    # the necessary array sizes
    E.info("getting number of variants")
    lines = 1
    with open(bim_file, "r") as bfile:
        for line in bfile.readlines():
            lines += 1

    E.info("%i variants found" % lines)
    
    # setup index arrays
    var_array = np.empty(lines, dtype=object)
    ref_alleles = np.empty(lines, dtype=object)
    pos_array = np.zeros(lines, dtype=np.int64)
    minor_alleles = np.empty(lines, dtype=object)
    idx = 0

    # find duplicates on position
    with open(bim_file, "r") as bfile:
        for line in bfile.readlines():
            line = line.rstrip("\n")
            varline = line.split("\t")
            var = varline[1]
            pos = int(varline[3])
            ref_allele = varline[-1]
            minor_allele = varline[-2]
            var_array[idx] = var
            ref_alleles[idx] = ref_allele
            minor_alleles[idx] = minor_allele
            pos_array[idx] = pos
            
            idx += 1
   
    # find duplicates using pandas series
    pos_series = pd.Series(pos_array)
    dup_last = pos_series[pos_series.duplicated(take_last=True)]
    dup_first = pos_series[pos_series.duplicated(take_last=False)]
    var_series = pd.Series(var_array)
    ref_series = pd.Series(ref_alleles)
    alt_series = pd.Series(minor_alleles)

    # a few variants have duplicate IDs - count these as duplicates
    # and add to the exclusion list - these won't be identified
    # based on shared position necessarily - force add them
    ref_first = ref_series[ref_series.duplicated(take_last=False)]
    ref_last = ref_series[ref_series.duplicated(take_last=True)]
    ref_dups = set(ref_first.index).union(ref_last.index)

    # union of take first and take last
    dup_all = set(dup_last.index).union(set(dup_first.index))
    dup_complete = dup_all.union(ref_dups)
    dup_idx = np.array([sx for sx in dup_complete])
    dup_idx.sort()

    # make a dataframe to hold all triallelic and duplicate variants
    dup_dict = {"SNP": var_series[dup_idx],
                "BP": pos_series[dup_idx],
                "REF": ref_series[dup_idx],
                "VAR": alt_series[dup_idx]}
    dup_df = pd.DataFrame(dup_dict)

    # some variants may have more than one ID/entry
    # step through using pandas groupby - group on position
    E.info("looking for duplicates and triallelic variants")
    tri_alleles = []
    dups_alleles = []
    overlap_vars = []
    for names, groups in dup_df.groupby(["BP"]):
        # if there is only one reference allele, indicates a
        # triallelic variant, otherwise its probably a duplicate
        # or overlaping INDEL and SNV
        var_lens = groups["VAR"].apply(len)
        if np.mean(var_lens) > 1:
            # probably overlapping variants, exclude, but report
            # separately
            over_vars = groups["SNP"].values.tolist()
            for ovs in over_vars:
                overlap_vars.append(ovs)

        elif len(set(groups["REF"])) == 1:
            tri_vars = groups["SNP"].values.tolist()
            for tri in tri_vars:
                tri_alleles.append(tri)
        else:
            dup_vars = groups["SNP"].values.tolist()
            for dup in dup_vars:
                dups_alleles.append(dup)

    E.info("%i triallelic variants found" % len(tri_alleles))
    E.info("%i duplicate position variants found" % len(dups_alleles))
    E.info("%i overlapping SNVs and INDELs found" % len(overlap_vars))
    
    return dups_alleles, tri_alleles, overlap_vars

def flagExcessHets(hets_file):
    '''
    Take output from Plink 1.9 --het command
    calculate heterozygosity rate and flag individuals
    with heterozygosity > 3 s.d. from the mean
    value.
    This assumes all individuals are from the same
    population, and thus form a homogenous cluster,
    with only outliers at the extremes.
    Visualise the data, if there are multiple apparent
    clusters then filter for ethnicity/ancestry first
    '''
    if hets_file.endswith("gz"):
        compression = "gzip"
    else:
        compression = None

    het_df = pd.read_table(hets_file, header=0, index_col=None,
                           sep="\t", compression=compression)
    nmiss = pd.Series(het_df.loc[:, "N(NM)"], dtype=np.float64)
    nhoms = het_df.loc[:, "O(HOM)"]
    het_df["het_rate"] = (nmiss - nhoms) / nmiss
    # get mean value and std, set upper and lower thresholds
    mean_het = np.mean(het_df.loc[:, "het_rate"].values)
    sd_het = np.std(het_df.loc[:, "het_rate"].values)

    upper = mean_het + (3 * sd_het)
    lower = mean_het - (3 * sd_het)

    hi_hets = het_df[het_df["het_rate"] > upper]
    lo_hets = het_df[het_df["het_rate"] < lower]

    E.info("%i individuals with high heterozygosity" % len(hi_hets))
    E.info("%i individuals with low heterozygosity" % len(lo_hets))

    hi_hets["exclude"] = "high_heterozygosity"
    lo_hets["exclude"] = "low_heterozygosity"
    all_flags = lo_hets.append(hi_hets)

    return all_flags
