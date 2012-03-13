\name{run.detsel}
\Rdversion{1.0}
\alias{run.detsel}
\title{Create Simulated Data}
\description{This runs the simulated data, using a coalescent-based algorithm.
}
\usage{
  run.detsel(example)
}
\arguments{
  \item{example}{a logical variable, which is TRUE if the user wants to run a toy example with the example file.}
}
\details{
Once the \code{\link{read.data}} command line has been executed, \code{\link{run.detsel}} executes the simulations. The user is first asked to provide the total number of simulations for the entire set of parameter values (default: 500000). For bi-allelic data, I recommend to run no less than 1000000 simulations to estimate correctly the \emph{P}-values for each empirical locus. With less than a million simulations, indeed, simulation tests have shown that the \emph{P}-values may be biased. The user is then asked to provide the average mutation rate, and the mutation model: type \sQuote{0} for the infinite allele model, where each mutation creates a new allelic state in the population; type \sQuote{1} for the stepwise mutation model, where each mutation consists in increasing or decreasing by one step, the size of the current allele; and type any integer \eqn{k} (with \eqn{k > 1}) for a k-allele model, where each mutation consists in drawing randomly one allele among k possible states, provided it is different from the current state. For example, for SNP data, type \sQuote{2}. Finally, the user is asked to provide the number of distinct sets of nuisance parameters (the default is a single set of parameters). Because of the uncertainty in the nuisance parameter values, it is recommended to perform simulations using different combinations of values for the ancestral population size, divergence time and bottleneck parameters. Then, the user is asked to provide as many sets of parameters as he/she indicated. Each set comprises four parameters that should be given in the following order: \eqn{t, N_0,t_0} and \eqn{N_e}. Here, the user must chose parameter values, including the mutation model, that correspond to his/her knowledge of the biological model.

The command line \code{\link{run.detsel}} creates a list of files named \sQuote{Pair_i_j_ni_nj.dat}, where i and j are the indices of populations pairs, and ni and nj are the sample sizes of populations i and j, respectively. Because some marker loci may have missing data, several \sQuote{Pair_i_j_ni_nj.dat} files may be created for a given pair of populations. Simulating the exact sample size for each locus is required to precisely calculate the empirical \emph{P}-values, especially for bi-allelic markers. Note that if negative multi-locus Fi estimates are observed for a pairwise comparison, then the simulations will not be run for that pair. Each line of the \sQuote{Pair_i_j_ni_nj.dat} files contains the locus-specific estimates of \eqn{F_1} and \eqn{F_2}, Weir and Cockerham's estimate of \eqn{F_{ST}} (Weir and Cockerham 1984), Nei's heterozygosity (\eqn{H_e}), and the number of alleles at that locus in the pooled sample. The command line run.detsel() also creates a file named out.dat that contains the estimates of the above statistics averaged over all the simulated data. In the file out.dat, each line corresponds to the pairwise analysis of populations i and j with sample sizes and ni and nj. Each line contains (in that order): the name of the output simulation file (\sQuote{Pair_i_j_ni_nj.dat}), the multi-locus estimates of \eqn{F_1} and \eqn{F_2}, and Weir and Cockerham's estimate of \eqn{F_{ST}} (Weir and Cockerham 1984). An important point to consider is to make sure that for each pairwise comparison, the   and   estimates averaged over the simulated data (in the file out.dat) closely match to the observed values in the real dataset (in the file infile.dat). If not, this suggests that the simulated datasets do not fit to the observed data, which urges to choose other parameter values for the nuisance parameters. 
}
\value{
The output files are saved in the current directory.
}
\examples{
## This is to generate an example file in the working directory.
make.example.files()

## This will read an input file named 'data.dat' that contains co-dominant markers,
## and a maximum allele frequency of 0.99 will be applied (i.e., by removing 
## marker loci in the observed and simulated datasets that have an allele with
## frequency larger than 0.99).
read.data(infile = 'data.dat',dominance = FALSE,maf = 0.99)

## The following command line executes the simulations:
run.detsel(example = TRUE)
}
\references{
Weir, B. S., and Cockerham, C. C. (1984) Estimating F-statistics for the analysis of population structure, \emph{Evolution} \bold{38}: 1358-1370.
}

