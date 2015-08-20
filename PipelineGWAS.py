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
        self.phenotypes = phenotypes
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
                                                               pf)]
            self.map_file = [mf for mf in infiles if re.search(".map",
                                                               mf)]

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
                                                               ff)]
            self.bim_file = [fb for fb in infiles if re.search(".bim",
                                                              fb)]
            self.bed_file = [bf for bf in infiles if re.search(".bed",
                                                               bf)]
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
                                                               gf)]
            self.sample_file = [sf for sf in infiles if re.search(".sample",
                                                                  sf)]
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
                                                                bg)]
            self.sample_file = [sf for sf in infiles if re.search(".sample",
                                                                  sd)]
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
                                                               vf)]
            # check files exist (i.e. are not the default None values)
            try:
                assert self.vcf_file
            except AssertionError:
                raise ValueError(".vcf file is missing, please "
                                 "specify")

        elif self.file_format == "bcf":
            self.bcf_file = [bv for bv in infiles if re.search(".bcf",
                                                               bv)]
            # check files exist (i.e. are not the default None values)
            try:
                assert self.bcf_file
            except AssertionError:
                raise ValueError(".bcf file is missing, please "
                                 "specify")

        elif self.file_format == "GRM_binary":
            self.id_file = [ig for ig in infiles if re.search(".grm.id",
                                                              ig)]
            self.n_file = [gn for gn in infiles if re.search(".grm.N",
                                                             gn)]
            self.bin_file = [gb for gb in infiles if re.search(".grm.bin",
                                                               gb)]
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
                                 "to use GCTA.  Please convert and "
                                 "try again.")

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
        elif file_format == "GRM_binary":
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

    def build_statement(self, infiles, outfile, threads=None):
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

    def genetic_relationship_matrix(self, shape, compression, metric):
        '''
        Calculate genomic pair-wise distance matrix between
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

        if metric in ["cov", "ibc2", "ibc3"]:
            state = self._matrices(matrix_type="grm", shape=shape,
                                   compression=compression, options=metric)
        else:
            E.info("%s metric not recognised.  Running with default Fhat1")
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
        elif filter_type == "exclude_snps":
            self._construct_filters(exclude=filter_value)
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
                      "snp_bp_range": " --from %s --to %s ",
                      "specific_snp": " --snp %s ",
                      "window_size": " --window %s ",
                      "exclude_snp": " --exclude-snp %s ",
                      "range_resolution": "--from-%s %s --to-%s %s ",
                      "covariates_file": " --filter %s ",
                      "covariate_filter": " %s ",
                      "covariate_column": " --mfilter %s ",
                      "missing_phenotype": " --prune "}

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
        Plink2 is capable of much more than just running basic association
        analyses.

        These include file processing, reformating, filtering, data summaries,
        PCA, clustering, GRM calculation (slow and memory intense), etc.

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
                                        "suppress_first": " --list-duplicate-vars suppress-first"}}

        for task, value in kwargs.iteritems():
            if task != "parameter":
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
            flag = " --make-rel "

        if options:
            statement.append(" ".join([flag, shape, compression, options]))
        else:
            statement.append(" ".join([flag, shape, compression]))

        return " ".join(statement)

    def build_statement(self, infiles, outfile, threads=None):
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

        # add output flag
        statement.append(" --out %s " % outfile)

        os.system(" ".join(statement))
