\name{DetSel}
\Rdversion{1.0.2}
\alias{DetSel}
\title{A R-package to Detect Marker Loci Responding to Selection}
\description{
In the new era of population genomics, surveys of genetic polymorphism (\sQuote{genome scans}) offer the opportunity to distinguish locus-specific from genome wide effects at many loci. Identifying presumably neutral regions of the genome that are assumed to be influenced by genome-wide effects only, and excluding presumably selected regions, is therefore critical to infer population demography and phylogenetic history reliably. Conversely, detecting locus-specific effects may help identify those genes that have been, or still are, targeted by natural selection. The software package DetSel has been developed to identify markers that show deviation from neutral expectation in pairwise comparisons of diverging populations. Recently, two major improvements have been made: the analysis of dominant markers is now supported, and the estimation of empirical \emph{P}-values has been implemented. These features, which are described below, have been incorporated into an R package, which replaces the stand-alone DetSel software package.
}
\author{Renaud Vitalis \email{vitalis@supagro.inra.fr}}
\references{
Vitalis, R., Dawson, K., and Boursot, P. (2001) Interpretation of variation across marker loci as evidence of selection, \emph{Genetics} \bold{158}, 1811--1823.

Vitalis, R., Dawson, K., Boursot, P., and Belkhir, K. (2003) DetSel 1.0: a computer program to detect markers responding to selection, \emph{Journal of Heredity} \bold{94}, 429--431.

Vitalis R. (2012) DETSEL: a R-package to detect marker loci responding to selection, \emph{in} Data Production and Analysis in Population Genomics (Pompanon F. and Bonin A., eds), pp. 277--293 \emph{Methods in Molecular Biology, vol. 888}, Humana Press, USA.
}