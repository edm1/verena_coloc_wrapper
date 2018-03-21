#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy (March 2018)
#
"""
Wrapper for Verena Zuber's colocalisation method

Assumptions:
    - GWAS and cov matrix are referring to the same strand
"""

import sys
import argparse
import pandas as pd
import numpy as np
import Finemap
import pickle

def main():

    # Get args
    args = parse_args()
    chrom, start, end = parse_range(args.range)

    #
    # Load and run 1-D finemapping on left dataset -----------------------------
    #

    # Load left sumstats
    print("Loading left sumstats...")
    left_z, left_efal, left_n = parse_sumstats(args.left_sumstats, chrom, start,
        end, args.left_rsidcol, args.left_chromcol, args.left_poscol,
        args.left_betacol, args.left_secol, args.left_effalcol, args.left_ncol,
        args.sep)
    print("  {0} variants loaded...".format(left_z.shape[0]))

    # Load left cov matrix
    print("Loading left cov matrix...")
    left_cov, left_cov_efal = parse_cov_matrix(args.left_cov, args.left_covmeta)
    print("  {0} variants loaded...".format(left_cov.shape[0]))

    # Harmonise left z-scores and cov matrix
    print("Harmonising left sumstat and cov matrix...")
    left_z, left_cov = harmonise_z_covmat(left_z, left_cov, left_efal,
                                          left_cov_efal)
    print("  {0} overlapping variants...".format(left_z.shape[0]))

    # Run finemap on left
    print("Running finemap on left...")
    left_res = Finemap.finemap(left_z.as_matrix(),
                               left_cov.as_matrix(),
                               left_n,
                               left_z.index.tolist(),
                               kmax=args.left_kmax,
                               kstart=args.left_kstart,
                               max_iter=args.maxiter,
                               prior=args.prior,
                               v_scale=args.v_scale,
                               g=args.g)
    # Write output as a pickled object and pandas df
    save_object(left_res, "{0}.left_configurations.pkl".format(args.outprefix))
    left_res_df = finemap_obj_to_df(left_res)
    left_res_df.to_csv("{0}.left_configurations.tsv".format(args.outprefix),
                       sep="\t", index=None)

    #
    # Load and run 1-D finemapping on 2nd (right) dataset ----------------------
    #

    # Skip if no 2nd dataset is provided
    if args.right_sumstats is None:
        print("No right dataset provided (--right_sumstats). Finishing here!")
        return 0

    # Load right sumstats
    print("Loading right sumstats...")
    right_z, right_efal, right_n = parse_sumstats(args.right_sumstats, chrom,
        start, end, args.right_rsidcol, args.right_chromcol, args.right_poscol,
        args.right_betacol, args.right_secol, args.right_effalcol,
        args.right_ncol, args.sep)
    print("  {0} variants loaded...".format(right_z.shape[0]))

    # Load right cov matrix
    print("Loading right cov matrix...")
    right_cov, right_cov_efal = parse_cov_matrix(args.right_cov,
                                                 args.right_covmeta)
    print("  {0} variants loaded...".format(right_cov.shape[0]))

    # Harmonise right z-scores and cov matrix
    print("Harmonising right sumstat and cov matrix...")
    right_z, right_cov = harmonise_z_covmat(right_z, right_cov, right_efal,
                                            right_cov_efal)
    print("  {0} overlapping variants...".format(right_z.shape[0]))

    # Run finemap on right
    print("Running finemap on right...")
    right_res = Finemap.finemap(right_z.as_matrix(),
                                right_cov.as_matrix(),
                                right_n,
                                right_z.index.tolist(),
                                kmax=args.right_kmax,
                                kstart=args.right_kstart,
                                max_iter=args.maxiter,
                                prior=args.prior,
                                v_scale=args.v_scale,
                                g=args.g)
    # Write output as a pickled object and pandas df
    save_object(right_res, "{0}.right_configurations.pkl".format(args.outprefix))
    right_res_df = finemap_obj_to_df(right_res)
    right_res_df.to_csv("{0}.right_configurations.tsv".format(args.outprefix),
                       sep="\t", index=None)

    #
    # Caluclate 2D joint posteriors (colocalisation)
    #

    print("Running colocalisation between left and right...")

    coloc_evidence, joint_res = left_res.joint_posterior(right_res)

    # Write table of joint posterior configurations
    joint_res_df = joint_posterior_to_df(joint_res)
    joint_res_df.to_csv("{0}.joint_configurations.tsv".format(args.outprefix),
                       sep="\t", index=None)
    save_object(joint_res, "{0}.joint_configurations.pkl".format(args.outprefix))
    # Write file containing single colocalisation evidence score
    outf = "{0}.joint_evidence.tsv".format(args.outprefix)
    with open(outf, "wb") as out_h:
        out_h.write("{0}\n".format(coloc_evidence).encode("utf8"))

    print("Finished!")

    return 0

def joint_posterior_to_df(obj):
    """ Converts POSTGAP joint posterior (colocalisation) in a pandas df
    Args:
        obj (TwoDConfigurationSample)
    Returns:
        Dataframe
    """
    # Add each configuration as a row
    rows = []
    for configuration in obj.configurations:
        snp_labels = ";".join([obj.labels[position] for position in configuration])
        index = obj.configurations[configuration]
        posterior  = obj.posterior[index]
        posterior1 = obj.posterior1[index]
        posterior2 = obj.posterior2[index]
        rows.append([snp_labels, posterior1, posterior2, posterior])
    df = pd.DataFrame(rows, columns=["snps", "left_posterior", "right_posterior", "joint_posterior"])
    df = df.sort_values("joint_posterior", ascending=False)
    return df

def finemap_obj_to_df(obj):
    """ Converts POSTGAP finemap object to a pandas df
    Args:
        obj (OneDConfigurationSample)
    Returns:
        Dataframe
    """
    # Add each configuration as a row
    rows = []
    for configuration in obj.configurations:
        index_of_configuration = obj.configurations[configuration]
        posterior = obj.posterior[index_of_configuration]
        prior = np.exp(obj.log_prior[index_of_configuration])
        logBF = obj.log_BF[index_of_configuration]
        snps = ';'.join([obj.labels[position] for position in configuration])
        rows.append([snps, prior, logBF, posterior])
    df = pd.DataFrame(rows, columns=["snps", "prior", "logBF", "posterior"])
    df = df.sort_values("posterior", ascending=False)

    return df

def save_object(obj, filename):
    """ Write object as a pickle
    Args:
        obj: any object
        filename: output file
    """
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
    return 0

def harmonise_z_covmat(z, cov, z_efal, cov_efal):
    """ Harmonise z-scores and covariance matrix. This includes:
        1. Take intersection of rsids between the two
        2. sort into the same orders
        3. if effect alleles are not the same, flip in the cov matrix
    Args:
        z (Series)
        cov (Dataframe)
        z_efal (Series)
        cov_efal (Series)
    Returns:
        z (Series), cov (Dataframe)
    """
    intersect = z.index[z.index.isin(cov.index)]
    z = z.loc[intersect]
    cov = cov.loc[intersect, intersect]

    # Make harmonisation vector to flip cov if effect alleles don't match
    vec = [1 if z_efal[snp] == cov_efal[snp] else -1 for snp in intersect]
    vec = pd.Series(vec, index=intersect)

    # Multiply the cov matrix by the harmonisation vector across both axes
    cov = cov.multiply(vec, axis=0)
    cov = cov.multiply(vec, axis=1)

    return z, cov

def parse_cov_matrix(in_cov, in_meta):
    """ Parses ldstore cov matrix and associated meta data file
    Args:
        in_cov (str): file containing matrix of correlations
        in_meta (str): file containing variant meta-data
    Returns:
        cov matrix (dataframe), effect alleles (dict)
    """
    cov = pd.read_csv(in_cov, sep=" ", header=None)
    meta = pd.read_csv(in_meta, sep=" ", header=0)
    cov.index = meta.RSID
    meta.index = meta.RSID
    cov.columns = meta.RSID

    # Remove duplicates
    isdupe = meta.index.duplicated()
    cov = cov.loc[~isdupe, ~isdupe]
    meta = meta[~isdupe]

    return cov, meta["A_allele"].to_dict()

def parse_sumstats(inf, chrom, start, end, rsid_col, chrom_col, pos_col,
                   beta_col, se_col, effal_col, n_col, sep):
    """ Parse a summary stat file and return zscores
    Args:
        inf (str):   input summary stat file
        chrom (str): chrom to extract
        start (int): position to extract from
        end (int):   position to extract to
        ...column names
    Returns:
        zscores (Series), effect alleles (dict), n (int)
    """
    # Load iteratively
    iter_csv = pd.read_csv(inf, sep=sep, header=0, iterator=True,
               chunksize=100000)
    iter_chunks = []
    for chunk in iter_csv:
        chunk[chrom_col] = chunk[chrom_col].astype(str)
        iter_chunks.append(chunk[((chunk[chrom_col] == chrom) &
                               (chunk[pos_col] >= start) &
                               (chunk[pos_col] <= end))])
    df = pd.concat(iter_chunks)
    # df.to_csv("temp/test_sumstat_slice.tsv", index=None, sep="\t") # DEBUG

    # Only keep rows with a valid rsid
    df = df[df[rsid_col].str.startswith("rs")]
    df.index = df[rsid_col]

    # Remove duplicates
    df = df.loc[~df.index.duplicated(), :]

    # Calculate z-scores
    zscores = df[beta_col] / df[se_col]

    # Get mean N
    n_mean = int(np.mean(df[n_col]))

    return zscores, df[effal_col].to_dict(), n_mean

def parse_range(range_str):
    """ Parse the chrom, start and end from the range string
    Args:
        range_str (str): In format chrom:start-end
    Returns:
        chrom (str), start (int), end (int)
    """
    chrom, start_end = range_str.split(":")
    start, end = start_end.split("-")
    return str(chrom), int(start), int(end)

def parse_args():
    """ Load command line args.
    """
    parser = argparse.ArgumentParser()
    # Input args
    parser.add_argument('--left_sumstats', metavar="<file>", help=('Summary statistics file'), type=str, required=True)
    parser.add_argument('--left_cov', metavar="<file>", help=("Covariance matrix - correlation structure between variants"), type=str, required=True)
    parser.add_argument('--left_covmeta', metavar="<file>", help=("Covariance matrix SNP info (from LDstore)"), type=str, required=True)
    parser.add_argument('--range', metavar="<str>", help=("Genomic range in format chrom:start-end"), type=str, required=True)
    parser.add_argument('--right_sumstats', metavar="<file>", help=('(Optional) Summary statistics file'), type=str, required=False)
    parser.add_argument('--right_cov', metavar="<file>", help=("Covariance matrix - leave blank if same as left_cov"), type=str, required=False)
    parser.add_argument('--right_covmeta', metavar="<file>", help=("Covariance matrix SNP info (from LDstore)"), type=str, required=False)
    # Output args
    parser.add_argument('--outprefix', metavar="<str>", help=("Output prefix"), type=str, required=True)
    # Finemapping args
    parser.add_argument('--left_kmax', metavar="<int>", help=('Maximum number of causal variants (default: 5)'), default=5, type=int)
    parser.add_argument('--left_kstart', metavar="<int>", help=('Full exploration of sets with #kstart causal variants (default: 1)'), default=1, type=int)
    parser.add_argument('--right_kmax', metavar="<int>", help=('Maximum number of causal variants (default: left_kmax)'), type=int)
    parser.add_argument('--right_kstart', metavar="<int>", help=('Full exploration of sets with #kstart causal variants (default: left_kstart)'), type=int)
    parser.add_argument('--v_scale', metavar="<float>", help=('Prior variance of the independence prior (default: 0.0025)'), default=0.0025, type=float)
    parser.add_argument('--maxiter', metavar="<int>", help=('Max iterations of stochastic search (default: 100000)'), default=100000, type=int)
    parser.add_argument('--prior', metavar="<str>", help=('Choice of "independence" or "gprior" (default: independence)'), default="independence", type=str, choices=["independence", "gprior"])
    parser.add_argument('--g', metavar="<str>", help=('g-parameter of the g-prior (default: BRIC). see Mixtures of g Priors for Bayesian Variable Selection Liang et al 2008'), default="BRIC", type=str, choices=["BRIC", "BIC", "RIC"])
    # File parsing args
    parser.add_argument('--sep', metavar="<str>", help=('Column sep (default: tab)'), type=str, default="\t")
    parser.add_argument('--left_rsidcol', metavar="<str>", help=('RSID column (default: rsid)'), default="rsid", type=str)
    parser.add_argument('--left_chromcol', metavar="<str>", help=('Chromosome column (default: chrom)'), default="chrom", type=str)
    parser.add_argument('--left_poscol', metavar="<str>", help=('Position column (default: pos)'), default="pos", type=str)
    parser.add_argument('--left_betacol', metavar="<str>", help=('Beta column (default: beta)'), default="beta", type=str)
    parser.add_argument('--left_secol', metavar="<str>", help=('Standard error column (default: se)'), default="se", type=str)
    parser.add_argument('--left_effalcol', metavar="<str>", help=('Effect allele column (default: effect_allele)'), default="effect_allele", type=str)
    parser.add_argument('--left_ncol', metavar="<str>", help=('Sample size column (default: n)'), default="n", type=str)
    parser.add_argument('--right_rsidcol', metavar="<str>", help=('RSID column (default: left_rsidcol)'), type=str)
    parser.add_argument('--right_chromcol', metavar="<str>", help=('Chromosome column (default: left_chromcol)'), type=str)
    parser.add_argument('--right_poscol', metavar="<str>", help=('Position column (default: left_poscol)'), type=str)
    parser.add_argument('--right_betacol', metavar="<str>", help=('Beta column (default: left_betacol)'), type=str)
    parser.add_argument('--right_secol', metavar="<str>", help=('Standard error column (default: left_secol)'), type=str)
    parser.add_argument('--right_effalcol', metavar="<str>", help=('Effect allele columns (default: left_effalcol)'), type=str)
    parser.add_argument('--right_ncol', metavar="<str>", help=('Sample size column (default: n)'), default="n", type=str)

    args = parser.parse_args()

    # If right argument is not set, set to same as left argument
    for rarg in ["right_cov", "right_covmeta", "right_kmax", "right_kstart",
                 "right_rsidcol", "right_chromcol", "right_poscol",
                 "right_betacol", "right_secol", "right_effalcol",
                 "right_ncol"]:
        if getattr(args, rarg) is None:
            larg_val = getattr(args, rarg.replace("right_", "left_"))
            setattr(args, rarg, larg_val)

    return args

if __name__ == '__main__':

    main()
