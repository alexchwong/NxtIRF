% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{example-NxtIRF-data}
\alias{example-NxtIRF-data}
\alias{mock_genome}
\alias{mock_gtf}
\alias{example_bams}
\alias{NxtIRF_example_bams}
\alias{NxtIRF_example_NxtSE}
\alias{get_mappability_exclusion}
\title{NxtIRF Examples}
\usage{
mock_genome(destination_path = tempdir())

mock_gtf(destination_path = tempdir())

get_mappability_exclusion(
  genome_type = c("hg38", "hg19", "mm10", "mm9"),
  destination_path = tempdir()
)

example_bams(destination_path = tempdir())

NxtIRF_example_bams()

NxtIRF_example_NxtSE()
}
\arguments{
\item{destination_path}{The directory to place the downloaded example files.
If left blank, gives the direct BiocFileCache path of the resource.}

\item{genome_type}{Either one of \code{hg38}, \code{hg19}, \code{mm10} or \code{mm9}}
}
\value{
See Functions section below.
}
\description{
Contains files that provides a workable example for the
NxtIRF package.\cr\cr
A mock reference, with genome sequence (FASTA) and gene annotation (GTF)
files are provided, based on the genes SRSF1, SRSF2, SRSF3, TRA2A, TRA2B,
TP53 and NSUN5, of which sequences are used to construct an artificial
chromosome Z. This was generated based on release-94 of Ensembl GRCh38 (hg38)
reference.\cr\cr
Also, there are 6 example bam files based on samples from the
Leucegene dataset (GSE67039). Bam files are constructed
based on the complete bam files of 6 samples from Leucegene,
subsetted by regions containing the 7 above genes. Then, the reads of these
subsetted BAMs were realigned to the mock reference using STAR.\cr\cr
Additionally, there are files for Mappability exclusion regions generated
using NxtIRF, suitable for use in generating references based on hg38,
hg19, mm10 and mm9 genomes.
}
\section{Functions}{
\itemize{
\item \code{mock_genome}: Makes a copy of the NxtIRF mock genome FASTA
file, in the given destination path

\item \code{mock_gtf}: Makes a copy of the NxtIRF mock gene 
annotation GTF file, in the given destination path

\item \code{get_mappability_exclusion}: Returns a copy of Mappability Exclusion 
BED file for the specified genome, in the given destination path

\item \code{example_bams}: Makes a copy of the Leucegene example 
BAM files, aligned to the NxtIRF mock genome, in the given destination path

\item \code{NxtIRF_example_bams}: Returns a 2-column data frame, containing 
sample names and sample paths (in tempdir()) of example BAM files

\item \code{NxtIRF_example_NxtSE}: Returns a (in-memory / realized) NxtSE object 
generated using the NxtIRF mock reference and example BAM files
}}

\examples{
# returns the location of the genome.fa file of the mock reference

mock_genome() 

# returns the location of the transcripts.gtf file of the mock reference

mock_gtf() 

# returns the location of the Mappability exclusion BED for hg38

get_mappability_exclusion("hg38") 

# returns the locations of the 6 example bam files

example_bams() 

# returns a data frame with the first column as sample names, and the 
# second column as BAM paths

NxtIRF_example_bams() 

# Returns a NxtSE object created by the example bams aligned to the

# mock NxtSE reference

NxtIRF_example_NxtSE() 


}
\references{
Generation of the mappability files was performed using NxtIRF using
a method analogous to that described in:

Middleton R, Gao D, Thomas A, Singh B, Au A, Wong JJ, Bomane A, Cosson B,
Eyras E, Rasko JE, Ritchie W.
IRFinder: assessing the impact of intron retention on mammalian gene
expression.
Genome Biol. 2017 Mar 15;18(1):51.
\doi{10.1186/s13059-017-1184-4}
}
\seealso{
\code{\link[=BuildReference]{BuildReference()}}, \code{\link[=IRFinder]{IRFinder()}}, \code{\link[=CollateData]{CollateData()}}, \code{\link[=MakeSE]{MakeSE()}}
}
\keyword{package}
