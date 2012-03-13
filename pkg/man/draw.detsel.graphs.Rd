\name{draw.detsel.graphs}
\Rdversion{1.0}
\alias{draw.detsel.graphs}
\title{Plot Graphical Outputs}
\description{This command plots graphical outputs for DetSel analyses.
}
\usage{
	draw.detsel.graphs(i,j,x.range,y.range,n.bins,m,alpha,pdf,outliers)
}
\arguments{
  \item{i}{population index}
  \item{j}{population index}
  \item{x.range}{the range of values in the x-axis, respectively, which takes the default values x.range = c(-1,1)}
  \item{y.range}{the range of values in the y-axis, respectively, which takes the default values y.range = c(-1,1)}
  \item{n.bins}{the size of the 2-dimensional array of n x n square cells used to bin the F_1 and F_2 estimates, which takes the default value n.bins = c(100,100)}
  \item{m}{the smoothing parameters of the ASH algorithm, which takes the default value m = c(2,2)}
  \item{alpha}{the alpha-level (hence 1 - alpha is the proportion of the distribution within the plotted envelope), which takes the default value alpha = 0.05}
  \item{pdf}{a logical variable, which is TRUE if the user wants graphics to be plotted in a pdf file}
  \item{outliers}{an optional vector that represents a list of candidate outliers, defined by the user}
}
\details{
Once the \code{\link{run.detsel}} and \code{\link{compute.p.values}} command lines have been executed, the function \code{\link{draw.detsel.graphs}} can be used to plot graphs with an estimation of the density of \eqn{F_1} and \eqn{F_2} estimates, as detailed in the appendix in Vitalis et al (2001). Note that if the arguments i and j are missing, then all the population pairs are plotted. It is noteworthy that our estimation of the density of the \eqn{F_1} and \eqn{F_2} estimates might be discontinuous, because of the discrete nature of the data (the allele counts). This is particularly true when the number of alleles upon which the distribution is conditioned is small. The command line \code{\link{draw.detsel.graphs}} produces as many conditional distributions per population pair as there are different allele numbers in the pooled sample. All the observed data points are plotted in each graph. The outlier loci are plotted with a star symbol. For the latter, the locus number (i.e., its rank in the data file) is provided on the graph. If the user choses not to provide a pre-defined list of outliers, then the outlier represent all the markers for which the empirical \emph{P}-value is below the threshold alpha-level,
}
\value{
The pdf files are created in the current directory.
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

## This compute empirical P-values, assuming a range of values from -1 to 1
## in both dimensions, a grid of 50 x 50 bins, and a smoothing parameter m = 3
## in both dimensions.
compute.p.values(x.range = c(-1,1),y.range = c(-1,1),n.bins = c(50,50),m = c(3,3))

## This plots (on the screen) the 99% confidence regions corresponding to the
## pair of populations 1 and 2, using a 50 x 50 2-dimensions array.
draw.detsel.graphs(i = 1,j = 2,n.bins = c(50,50),alpha = 0.01,pdf = FALSE) 
}
\references{
	Vitalis, R., Dawson, K., and Boursot, P. (2001) Interpretation of variation across marker loci as evidence of selection, \emph{Genetics} \bold{158}: 1811-1823.
}

