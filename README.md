# verena_coloc_wrapper
Wrapper scripts for Verena's (yet unnamed) GWAS co-localisation methods

### Set up environment
```
# Install dependencies into isolated environment
conda env create -n vloc --file environment.yaml

# Activate environment
source activate vloc
```

### Usage
The wrapper will run 1D-finemapping if only `--left_sumstats` dataset is specified. If a second dataset `--right_sumstats` is specified, it will run 1D-finemapping, followed by co-localisation between the left and right datasets.

```
$ python scripts/vloc.py --help
usage: vloc.py [-h] --left_sumstats <file> --left_cov <file> --left_covmeta
               <file> --range <str> [--right_sumstats <file>]
               [--right_cov <file>] [--right_covmeta <file>] --outprefix <str>
               [--left_kmax <int>] [--left_kstart <int>] [--right_kmax <int>]
               [--right_kstart <int>] [--v_scale <float>] [--maxiter <int>]
               [--prior <str>] [--g <str>] [--sep <str>]
               [--left_rsidcol <str>] [--left_chromcol <str>]
               [--left_poscol <str>] [--left_betacol <str>]
               [--left_secol <str>] [--left_effalcol <str>]
               [--left_ncol <str>] [--right_rsidcol <str>]
               [--right_chromcol <str>] [--right_poscol <str>]
               [--right_betacol <str>] [--right_secol <str>]
               [--right_effalcol <str>] [--right_ncol <str>]

optional arguments:
  -h, --help            show this help message and exit
  --left_sumstats <file>
                        Summary statistics file
  --left_cov <file>     Covariance matrix - correlation structure between
                        variants
  --left_covmeta <file>
                        Covariance matrix SNP info (from LDstore)
  --range <str>         Genomic range in format chrom:start-end
  --right_sumstats <file>
                        (Optional) Summary statistics file
  --right_cov <file>    Covariance matrix - leave blank if same as left_cov
  --right_covmeta <file>
                        Covariance matrix SNP info (from LDstore)
  --outprefix <str>     Output prefix
  --left_kmax <int>     Maximum number of causal variants (default: 5)
  --left_kstart <int>   Full exploration of sets with #kstart causal variants
                        (default: 1)
  --right_kmax <int>    Maximum number of causal variants (default: left_kmax)
  --right_kstart <int>  Full exploration of sets with #kstart causal variants
                        (default: left_kstart)
  --v_scale <float>     Prior variance of the independence prior (default:
                        0.0025)
  --maxiter <int>       Max iterations of stochastic search (default: 100000)
  --prior <str>         Choice of "independence" or "gprior" (default:
                        independence)
  --g <str>             g-parameter of the g-prior (default: BRIC). see
                        Mixtures of g Priors for Bayesian Variable Selection
                        Liang et al 2008
  --sep <str>           Column sep (default: tab)
  --left_rsidcol <str>  RSID column (default: rsid)
  --left_chromcol <str>
                        Chromosome column (default: chrom)
  --left_poscol <str>   Position column (default: pos)
  --left_betacol <str>  Beta column (default: beta)
  --left_secol <str>    Standard error column (default: se)
  --left_effalcol <str>
                        Effect allele column (default: effect_allele)
  --left_ncol <str>     Sample size column (default: n)
  --right_rsidcol <str>
                        RSID column (default: left_rsidcol)
  --right_chromcol <str>
                        Chromosome column (default: left_chromcol)
  --right_poscol <str>  Position column (default: left_poscol)
  --right_betacol <str>
                        Beta column (default: left_betacol)
  --right_secol <str>   Standard error column (default: left_secol)
  --right_effalcol <str>
                        Effect allele columns (default: left_effalcol)
  --right_ncol <str>    Sample size column (default: n)
```
