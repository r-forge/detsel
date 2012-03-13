\name{read.data}
\Rdversion{1.0}
\alias{read.data}
\title{Read Data}
\description{
Read the data file in DetSel format.
}
\usage{
  read.data(infile,dominance,maf,a,b)
}
\arguments{
  \item{infile}{An input file in DetSel format.}
  \item{dominance}{A a logical variable, which is FALSE if co-dominant data are considered (e.g., microsatellite markers, SNPs, etc.), or TRUE, if bi-allelic dominant data are considered (e.g., AFLPs).}
  \item{maf}{The maximum allele frequency (the frequency of the most frequent allele over the full sample) to be considered in both the input file and the simulated data.}
  \item{a,b}{The parameters for the beta prior distribution, used in Zhivotovsky's (1999) Bayesian method to compute the underlying allele frequencies. The default values are a = b = 0.25, as suggested by Mark A. Beaumont in the DFdist manual, yet the user may alternatively chose to use Zhivotovsky's equation (13) to compute estimates of a and b from the data. Note that neither the parameter a nor the parameter b are not needed if dominance = FALSE.}
}
\details{
The input file should be a space- or tab-delimited ASCII text file. The first line is a 0 / 1 indicator. \sQuote{0} indicates that the data matrix for each locus is a populations x alleles matrix; \sQuote{1} indicates that the data matrix for each locus is an alleles x populations matrix. The second line contains the number of populations. The third line contains the number of loci. Then, the data for each locus consists in the number of alleles at that locus, followed by the data matrix at that locus, with each row corresponding to the same allele (if the indicator variable is 1) or to the same population (if the indicator variable is 0). For dominant data, the data consists in the number of genotypes, not the number of alleles. It is important to note that the frequency of the homozygote individuals for the recessive allele appear first in either the rows or columns of the data matrix. In the following example, the data consists in 2 populations and 2 loci, with 5 alleles at the first locus and 8 alleles at the second locus.

	0\cr
	2\cr
	2\cr

	5\cr
	1	0	4	10	5\cr
	0	1	13	0	6\cr

	8\cr
	3	1	1	0	0	0	1	14\cr
	6	0	2	1	2	5	2	2\cr
	
Spaces and blank lines can be included as desired.

For dominant data, it is important to note that the frequency of the homozygote individuals for the recessive allele appears first in either the rows or columns of the data matrix.

The command line \code{\link{read.data}} creates a file named \sQuote{infile.dat}, a file named \sQuote{sample_sizes.dat} and a set of files named \sQuote{plot_i_j.dat} where \eqn{i} and \eqn{j} correspond to population numbers, so that each file \sQuote{plot_i_j.dat} corresponds to the pairwise analysis of populations \eqn{i} and j. In the file infile.dat, each line corresponds to the pairwise analysis of populations \eqn{i} and \eqn{j}. Each line contains (in that order): the name of the output simulation file, the numbers \eqn{i} and \eqn{j}, the multi-locus estimates of \eqn{F_1} and \eqn{F_2}, and Weir and Cockerham's (1984) estimate of \eqn{F_{ST}}. The file sample_sizes.dat contains sample sizes information, for internal use only. In the files \sQuote{plot_i_j.dat}, each line corresponds to one locus observed in the data set. Each line contains (in that order): the locus-specific estimates of \eqn{F_1} and \eqn{F_2}, Weir and Cockerham's (1984) estimate of \eqn{F_{ST}}, Nei's heterozygosity (\eqn{H_e}), the number of alleles at that locus in the pooled sample, and the rank of the locus in the data set.}
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
}
\references{
Weir, B. S., and Cockerham, C. C. (1984) Estimating F-statistics for the analysis of population structure, \emph{Evolution} \bold{38}: 1358-1370.

Zhivotovsky, L. A. (1999) Estimating population structure in diploids with multilocus dominant DNA markers, \emph{Molecular Ecology} \bold{8}, 907-913
}

