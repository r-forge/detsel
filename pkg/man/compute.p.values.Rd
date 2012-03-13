\name{compute.p.values}
\Rdversion{1.0}
\alias{compute.p.values}
\title{Compute Empirical p-values}
\description{This command compute empirical p-values.
}
\usage{
  compute.p.values(x.range,y.range,n.bins,m)
}
\arguments{
  \item{x.range}{the range of values in the x-axis, respectively, which takes the default values x.range = c(-1,1)}
  \item{y.range}{the range of values in the y-axis, respectively, which takes the default values y.range = c(-1,1)}
  \item{n.bins}{the size of the 2-dimensional array of n x n square cells used to bin the F_1 and F_2 estimates, which takes the default value n.bins = c(100,100)}
  \item{m}{the smoothing parameters of the ASH algorithm, which takes the default value m = c(2,2)}
}
\details{
compute.p.values(x.range,y.range,n.bins,m) produces an output file, named \sQuote{P-values_i_j.dat}, with the \emph{P}-value associated with each observation. To that end, the cumulative distribution function (CDF) is evaluated empirically from the joint distribution of all the pairwise observations (\eqn{F_1},\eqn{F_2}) within the simulated dataset. Then, the empirical \emph{P}-value for a given marker locus i is calculated as one minus the CDF evaluated at locus i. For multi-allelic markers, the joint distribution of all the pairwise observations (\eqn{F_1},\eqn{F_2}) within the simulated dataset is computed from a 2-dimensional array, where the (\eqn{F_1},\eqn{F_2}) pairs are binned, and then smoothed using the Average Shifted Histogram (ASH) algorithm (Scott 1992) as implemented in the "ash" R package. Because the distribution of (\eqn{F_1},\eqn{F_2}) estimates for bi-allelic markers is discontinuous with many ties, the CDF is computed instead by enumerating all (\eqn{F_1},\eqn{F_2}) pairs in the simulated data.
}
\value{
The output files are saved in the working directory.
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

## This compute empirical \emph{P}-values, assuming a range of values from -1 to 1
## in both dimensions, a grid of 50 x 50 bins, and a smoothing parameter m = 3
## in both dimensions.
compute.p.values(x.range = c(-1,1),y.range = c(-1,1),n.bins = c(50,50),m = c(3,3))
}
\references{
	Scott, D. W. (1992) Multivariate density estimation: theory, practice, and visualization, John Wiley, New York.
}

