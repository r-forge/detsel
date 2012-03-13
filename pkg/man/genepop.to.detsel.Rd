\name{genepop.to.detsel}
\Rdversion{1.0}
\alias{genepop.to.detsel}
\title{Convert Input File}
\description{This command converts a data file in Genepop format (see Raymond and Rousset 1995; Rousset 2007), which name can be specified using the infile argument into a data file in DetSel format, which name can be specified using the outfile argument.
}
\usage{
  genepop.to.detsel(infile,outfile = 'data.dat')
}
\arguments{
  \item{infile}{An input file in Genepop format.}
  \item{outfile}{An output file in DetSel format.}
}
\details{
	This command is only available for co-dominant data.
    The output file is a space-delimited ASCII text file. The first line is a 0 / 1 indicator. \sQuote{0} indicates that the data matrix for each locus is a populations x alleles matrix; \sQuote{1} indicates that the data matrix for each locus is an alleles x populations matrix. The second line contains the number of populations. The third line contains the number of loci.
    Then, the data for each locus consists in the number of alleles at that locus, followed by the data matrix at that locus, with each row corresponding to the same allele (if the indicator variable is 1) or to the same population (if the indicator variable is 0).
    For dominant data, the data consists in the number of genotypes, not the number of alleles. It is important to note that the frequency of the homozygote individuals for the recessive allele appear first in either the rows or columns of the data matrix. In the following example, the data consists in 2 populations and 2 loci, with 5 alleles at the first locus and 8 alleles at the second locus.

    0\cr
    2\cr
    2\cr

    5\cr
    1	0	4	10	5\cr
    0	1	13	0	6\cr

    8\cr
	3	1	1	0	0	0	1	14\cr
	6	0	2	1	2	5	2	2\cr
}
\value{
The output file is saved in the working directory.
}
\examples{
## This is to generate an example file in the working directory.
make.example.files()
## This is to convert the example file in genepop format named 'genepop.dat',
## into a file in DetSel format named 'converted_data.dat'
genepop.to.detsel(infile = 'data.gen',outfile = 'data-converted.dat')
}
\references{
Raymond, M., and Rousset, F. (1995) Genepop (version 1.2): population genetics software for exact tests and ecumenicism, \emph{Journal of Heredity} \bold{86}: 248-249.

Rousset, F. (2008) genepop'007: a complete re-implementation of the genepop software for Windows and Linux, \emph{Molecular Ecology Notes} \bold{8}: 103-106.
}

