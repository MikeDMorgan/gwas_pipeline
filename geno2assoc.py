'''
geno2assoc.py - perform genome-wide association analyses, filtering and summary
===============================================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Perform task associated with genome-wide analyses of genotyping data.

Usage
-----

.. Example use case

Example::

   python geno2assoc.py

Type::

   python geno2assoc.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import PipelineGWAS as gwas
import re
import random


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--program", dest="program", type="choice",
                      choices=["plink2", "gcta"],
                      help="program to execute genome-wide analysis")

    parser.add_option("--input-file-pattern", dest="infile_pattern", type="string",
                      help="file prefix that identifies a group of files")

    parser.add_option("--input-file-format", dest="file_format", type="choice",
                      choices=["plink", "plink_binary", "oxford",
                               "oxford_binary", "vcf", "GRM_binary",
                               "GRM_gz"],
                      help="format of input files")

    parser.add_option("--phenotypes-file", dest="pheno_file", type="string",
                      help="text file of additional phenotypes")

    parser.add_option("--pheno", dest="pheno", type="string",
                      help="either phenotype file column header or number")

    parser.add_option("--covariates-file", dest="covariate_file", type="string",
                      help="file containing covariates")

    parser.add_option("--covariate-column", dest="covar_col", type="string",
                      help="column number(s) or header(s) to include in "
                      "association model")

    parser.add_option("--association-model", dest="assoc_model", type="choice",
                      choices=["recessive", "dominant", "genotype"],
                      help="model to report from association analysis")

    parser.add_option("--method", dest="method", type="choice",
                      choices=["association", "summary", "format", "matrix",
                               "reml", "pca", "lmm", "simulation"],
                      help="method to apply to genome-wide data")

    parser.add_option("--principal-components", dest="num_pcs", type="int",
                      help="the number of principal components to output")

    parser.add_option("--matrix-shape", dest="matrix_shape", type="choice",
                      choices=["triangle", "square", "square0"],
                      help="output matrix shape.", default="triangle")

    parser.add_option("--matrix-compression", dest="matrix_compress", type="choice",
                      choices=["gz", "bin", "bin4"],
                      help="compression to apply to output matrix file",
                      default="gz")

    parser.add_option("--matrix-form", dest="matrix_form", type="choice",
                      choices=["distance", "grm"],
                      help="type of relationship matrix to calculate")

    parser.add_option("--matrix-metric", dest="matrix_metric", type="choice",
                      choices=["fhat", "cov", "ibc2", "ibc3", "ibs",
                               "genomic", "hamming"],
                      help="value to calculate for diagonal elements of the "
                      "grm. Default is fhat for grm and hamming for distance.")

    parser.add_option("--matrix-options", dest="matrix_options", type="string",
                      help="modifiers of matrix output, see plink documentation "
                      "for details")
    
    parser.add_option("--association-method", dest="assoc_method", type="choice",
                      choices=["linear", "logistic", "assoc", "qassoc"],
                      help="association analysis to run")

    parser.add_option("--permutation", dest="permutation", action="store_true",
                      help="perform association testing by permutation analysis")

    parser.add_option("--repeats", dest="n_perms", type="int",
                      help="number of repetitions for permutation analysis")

    parser.add_option("--association-options", dest="assoc_option", type="string",
                      help="association analysis modifiers")

    parser.add_option("--format-method", dest="format_method", type="choice",
                      choices=["change_format", "change_missing_values",
                               "update_variants", "update_samples", "flip_strands", 
                               "flip_scan", "sort", "merge", "find_duplicates"],
                      help="file formatting to apply to input files")

    parser.add_option("--format-parameter", dest="format_param", type="string",
                      help="formatting parameter, where appropriate")

    parser.add_option("--reformat-type", dest="reformat", type="choice",
                      choices=["plink", "plink_binary", "oxford", "oxford_binary"],
                      help="new format of input files to be reformatted to")

    parser.add_option("--apply-missing", dest="apply_missing", type="choice",
                      choices=["genotype", "phenotype"],
                      help="genotype or phenotype missing values to alter")

    parser.add_option("--update-variant-attribute", dest="variant_update", type="choice",
                      choices=["variant_ids", "missing_id", "chromosome", "centimorgan",
                               "name", "alleles", "map"],
                      help="update variant attributes")

    parser.add_option("--update-sample-attribute", dest="sample_update", type="choice",
                      choices=["sample_ids", "parents", "gender"],
                      help="sample attributes to be updated")

    parser.add_option("--strand-flip-subset", dest="flip_subset", action="store_true",
                      help="apply strand flipping to a subset of samples")

    parser.add_option("--flip-scan-type", dest="scan_param", type="choice",
                      choices=["default", "window", "threshold"],
                      help="strand flipping scan to apply to SNPs")

    parser.add_option("--sort-type", dest="sort_type", type="choice",
                      choices=["none", "natural", "ascii", "file"],
                      help="sort type to input files")

    parser.add_option("--merge-file-format", dest="merge_format", type="choice",
                      choices=["plink", "binary_plink"],
                      help="format of input files to be merged")

    parser.add_option("--merge-mode", dest="merge_mode", type="choice",
                      choices=["default", "original_missing", "new_nonmissing",
                               "no_overwrite", "force", "report_all",
                               "report_nonmissing"],
                      help="merge mode to apply to dealing with merge conflicts")

    parser.add_option("--duplicates-method", dest="dup_method", type="choice",
                      choices=["same_ref", "id_match", "suppress_first"],
                      help="method for identifying and dealing with duplicate "
                      "variants")

    parser.add_option("--summary-method", dest="summary_method", type="choice",
                      choices=["allele_frequency", "missing_data", "hardy_weinberg",
                               "mendel_errors", "inbreeding", "gender_checker",
                               "wrights_fst"],
                      help="summary statistics to calculate")

    parser.add_option("--summary-parameter", dest="sum_param", type="string",
                      help="optional parameters that can be passed to summary "
                      "statistics methods")

    parser.add_option("--genotype-rate", dest="filt_genotype_rate", type="string",
                      help="genotyping rate threshold.  SNPs below this threshold "
                      "will be excluded from analysis")

    parser.add_option("--indiv-missing", dest="filt_missingness", type="string",
                      help="individual missingness rate.  Individuals below "
                      "this threshold will be excluded from analysis")

    parser.add_option("--hardy-weinberg", dest="filt_hwe", type="string",
                      help="hardy-weinberg p-value threshold for SNPs.  SNPs "
                      "with a 2df chisquared p-value below this will be "
                      "filtered out")

    parser.add_option("--min-allele-frequency", dest="filt_min_allele_frequency",
                      type="string",
                      help="only include SNPs with an allele frequency equal to "
                      "or above this threshold")

    parser.add_option("--max-allele-frequency", dest="filt_max_allele_frequency",
                      type="string",
                      help="only include SNPs with an allele frequency equal to "
                      "or below this threshold")

    parser.add_option("--mendelian-error", dest="filt_mendelian_error", type="string",
                      help="exclude individuals/trios with mendelian errors that "
                      "exceed this value")

    parser.add_option("--keep-individuals", dest="filt_keep", type="string",
                      help="a file containing individuals IDs to keep, "
                      "one per row")

    parser.add_option("--remove-invidivals", dest="filt_remove", type="string",
                      help="a file of individual IDs to remove, one per row")

    parser.add_option("--min-quality-score", dest="filt_min_qaul_score", type="string",
                      help="reset the minimum low bound of quality scores for "
                      "variants in a VCF file.  Default is 0")

    parser.add_option("--max-quality-score", dest="filt_max_qual_score", type="string",
                      help="reset the maximum upper bound of quality scores for "
                      "a VCCF file.  Default is Inf")

    parser.add_option("--allow-no-gender", dest="filt_allow_no_sex", type="string",
                      help="allow individuals with gender missing")

    parser.add_option("--enforce-gender", dest="filt_enforce_sex", type="string",
                      help="only include individuals with non-missing gender "
                      "information")

    parser.add_option("--subset-filter", dest="filt_subset_filter", type="choice",
                      choices=["cases", "controls", "males", "females",
                               "founders", "nonfounders"],
                      help="only apply filters to the specific subset of "
                      "individuals supplied")

    parser.add_option("--extract-snps", dest="filt_extract", type="string",
                      help="text file of variant IDs to include in the analysis, "
                      "ignoring all others")

    parser.add_option("--exclude-snps", dest="filt_exclude", type="string",
                      help="a file of variant IDs to exclude from analysis")

    parser.add_option("--restrict-chromosome", dest="filt_chromosome", type="string",
                      help="restict analysis to either a single chromosome, "
                      "or a comma-separated list of chromosomes")

    parser.add_option("--exclude-chromosomes", dest="filt_exclude_chromosome",
                      type="string", help="exclude all variants on these "
                      "chromosome(s)")

    parser.add_option("--autosome-only", dest="filt_autosome", action="store_true",
                      help="if present only autosomal variants will be analysed")

    parser.add_option("--pseudo-autosome", dest="filt_pseudo_autosome", action="store_true",
                      help="include on the pseudo-autosomal region of chromosome X")

    parser.add_option("--ignore-indels", dest="filt_ignore_indels", action="store_true",
                      help="only include bi-allelic single nucleotide "
                      "variants in analysis")

    parser.add_option("--snp-range", dest="filt_snp_bp_range", type="string",
                      help="comma separated list of from, to genome co-ordinates "
                      "within which to include variants for analysis")

    parser.add_option("--snp-id", dest="filt_specific_snp", type="string",
                      help="include a single snp in the analysis given by "
                      "it's variant ID.")

    parser.add_option("--exclude-variant", dest="filt_exclude_snp", type="string",
                      help="exclude a single variant from the analysis, "
                      "given by it's variant ID")

    parser.add_option("--covariate-filter", dest="filt_covariate_filter", type="string",
                      help="covariate column headers or column numbers on which "
                      "to filter on. Requries --covariate-file")

    parser.add_option("--filter-parameter", dest="param", type="string",
                      help="parameter values to be passed to filtering function") 

    parser.add_option("--window-size", dest="window_size", type="string",
                      help="alters the behaviour of the --snp-range and "
                      "--include/exclude snp options.  variants within +/- "
                      "half * window_size (kb) are included")

    parser.add_option("--range-resolution", dest="range_resolution", type="choice",
                      choices=["bp", "kb", "mb"],
                      help="alters the (from, to) range resolution to either bp, "
                      "kb or mb")

    parser.add_option("--output-file-pattern", dest="out_pattern", type="string",
                      help="output file pattern prefix. file suffixes are dependent "
                      "on the task executed")

    parser.add_option("--threads", dest="threads", type="int",
                      help="the number of threads to use for multi-threaded "
                      "processes")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    parser.set_defaults(sum_param=None,
                        dup_method="same_ref",
                        n_perms=None,
                        permutation=False,
                        matrix_shape="triangle",
                        matrix_options=None,
                        matrix_compress="gz",
                        random_seed=random.randint(0, 19999))

    if not options.infile_pattern:
        infiles = (argv[-1]).split(",")
    else:
        infiles = options.infile_pattern

    # create a new filegroup object
    geno_files = gwas.FileGroup(files=infiles,
                                file_format=options.file_format,
                                genotype_format="imputed")

    # add FileGroup object to the gwas program object
    if options.program == "plink2":
        gwas_object = gwas.Plink2(files=geno_files)
        gwas_object.program_call(infiles=geno_files,
                                 outfile=options.out_pattern)
    elif options.program == "gcta":
        gwas_object = gwas.GCTA(files=geno_files)
        gwas_object.program_call(infiles=geno_files,
                                 outfile=options.out_pattern)
    else:
        pass

    # collect filtering options from options
    opt_dict = options.__dict__
    filter_keys = [fx for fx in opt_dict.keys() if re.search("filt", fx)]
    filter_dict = {k: options.__dict__[k] for k in filter_keys if opt_dict[k]}

    # iteratively add filters to GWASProgram object
    for fkey in filter_dict:
        filt_key = fkey.lstrip("filt_")
        filter_value = filter_dict[fkey]
        gwas_object.filter_genotypes(filter_type=filt_key,
                                     filter_value=filter_value)

    # handle summary statistics
    if options.method == "summary":
        if options.summary_method == "allele_frequency":
            gwas_object._output_statistics(allele_frequency=options.sum_param)
        elif options.summary_method == "hardy_weinberg":
            gwas_object._output_statistics(hardy_weinberg=options.sum_param)
        elif options.summary_method == "missing_data":
            gwas_object._output_statistics(missing_data=options.sum_param)
        elif options.summary_method == "mendel_errors":
            gwas_object._output_statistics(mendel_errors=options.sum_param)
        elif options.summary_method == "inbreeding":
            gwas_object._output_statistics(inbreeding=options.sum_param)
        elif options.summary_method == "gender_checker":
            gwas_object._output_statistics(gender_checker=options.sum_param)
        elif options.summary_method == "wrights_fst":
            gwas_object._output_statistics(wrights_fst=options.sum_param)
        else:
            pass
    elif options.method == "pca":
        gwas_object.PCA(n_pcs=options.num_pcs)
    elif options.method == "association":
        gwas_object.run_association(association=options.assoc_method,
                                    permutation=options.permutation,
                                    n_perms=options.n_perms,
                                    random_seed=options.random_seed)
    elif options.method == "format":
        if options.format_method == "change_format":
            gwas_object._run_tasks(change_format=options.reformat,
                                   parameter=options.format_param)
        elif options.format_method == "change_missing_values":
            gwas_object._run_tasks(change_missing_values=options.apply_missing,
                                   parameter=options.format_param)
        elif options.format_method == "update_variants":
            gwas_object._run_tasks(update_variants=options.variant_update,
                                   parameter=options.format_param)
        elif options.format_method == "update_samples":
            gwas_object._run_tasks(update_samples=options.sample_update,
                                   parameter=options.format_param)
        elif options.format_method == "flip_strands":
            if options.flip_subset:
                gwas_object._run_tasks(flip_strands="subset",
                                       parameter=options.format_param)
            else:
                gwas_object._run_tasks(flip_strands="all_samples",
                                       parameter=options.format_param)
        elif options.format_method == "flip_scan":
            gwas_object._run_tasks(flip_scan=options.scan_param,
                                   parameter=options.format_param)
        elif options.format_method == "sort":
            gwas_object._run_tasks(sort=options.sort_type,
                                   parameter=options.format_param)
        elif options.format_method == "merge":
            if options.merge_mode:
                gwas_object._run_tasks(merge_mode=options.merge_mode,
                                   parameter=options.format_param)
            else:
                gwas_object._run_tasks(merge=options.merge_format,
                                   parameter=options.format_param)
        elif options.format_method == "find_duplicates":
            gwas_object._run_tasks(find_duplicates=options.dup_method,
                                   parameter=options.format_param)
        else:
            pass
    elif options.method == "matrix":
        if options.matrix_form == "distance":
            if options.matrix_metric == "hamming":
                gwas_object.hamming_matrix(shape=options.matrix_shape,
                                           compression=options.matrix_compress,
                                           options=options.matrix_options)
            elif options.matrix_metric == "ibs":
                gwas_object.ibs_matrix(shape=options.matrix_shape,
                                       compression=options.matrix_compress,
                                       options=options.matrix_options)
            elif options.matrix_metric == "genomic":
                gwas_object.genome_matrix(shape=options.matrix_shape,
                                          compression=options.matrix_compress,
                                          options=options.matrix_options)
        elif options.matrix_form == "grm":
            gwas_object.genetic_relationship_matrix(shape=options.matrix_shape,
                                                    compression=options.matrix_compress,
                                                    metric=options.matrix_metric,
                                                    options=options.matrix_options)
    else:
        pass

    gwas_object.build_statement(infiles=geno_files,
                                outfile=options.out_pattern,
                                threads=options.threads)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
