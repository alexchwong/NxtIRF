% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CollateData.R
\name{CollateData}
\alias{CollateData}
\title{Processes data from IRFinder output}
\usage{
CollateData(
  Experiment,
  reference_path,
  output_path,
  IRMode = c("SpliceOverMax", "SpliceMax"),
  samples_per_block = 16,
  n_threads = 1
)
}
\arguments{
\item{Experiment}{A 2-column data frame (generated by \code{\link[=Find_Samples]{Find_Samples()}}),
with the first column designating the sample names, and the 2nd column
containing the primary IRFinder output file (of type \code{.txt.gz}). A
third optional column can contain the coverage files of the corresponding
samples. NB: all other columns are ignored.\cr\cr
We recommend using the output of \code{\link[=Find_IRFinder_Output]{Find_IRFinder_Output()}} here.}

\item{reference_path}{The path to the reference generated by
\code{\link[=BuildReference]{BuildReference()}}}

\item{output_path}{The path for the output files to be generated by this
function.}

\item{IRMode}{The algorithm to calculate 'splice abundance' in IR
quantification. The original algorithm by Middleton et al (2017) proposes
\code{SpliceMax}, which calculates the number of mapped splice events
that share the boundary coordinate of either the left or right flanking
exon (SpliceLeft, SpliceRight) and defines splice abundance as the larger
of the two values. NxtIRF proposes a new algorithm, \code{SpliceOverMax},
to account for the possibility that the major isoform shares neither
boundary, but arises from either of the flanking "exon islands". Exon
islands are contiguous regions covered by exons from any transcript
(except those designated as \code{retained_intron} or
\code{sense_intronic}), and are separated by
obligate intronic regions (genomic regions that are introns for all
transcripts). Note for introns that are internal to a single exon island
(i.e. akin to "known-exon" introns from IRFinder), \code{SpliceOverMax}
uses \link[GenomicRanges]{findOverlaps} to summate competing mapped
splice reads.}

\item{samples_per_block}{How many samples to process per thread.}

\item{n_threads}{The number of threads to use. On low
memory systems, reduce the number of \code{n_threads} and \code{samples_per_block}}
}
\value{
None. \code{CollateData()} writes to the directory given by
\code{output_path}
}
\description{
CollateData unifies a list of IRFinder output files belonging to an
experiment. It is assumed every sample is analysed using the same IRFinder
reference. The combination of junction counts and IR quantification from
IRFinder is used to calculate percentage spliced in (PSI) of alternative
splice events, and percent intron retention (PIR) of retained introns. Also,
QC information is extracted, and data is organised for quick
downstream access such as \code{\link[=MakeSE]{MakeSE()}}.
}
\examples{
bams = NxtIRF_example_bams()
IRFinder(bams$BAM, bams$sample,
  reference_path = file.path(tempdir(), "Reference"),
  output_path = file.path(tempdir(), "IRFinder_output")
)
expr = Find_IRFinder_Output(file.path(tempdir(), "IRFinder_output"))
CollateData(expr, 
  reference_path = file.path(tempdir(), "Reference"),
  output_path = file.path(tempdir(), "NxtIRF_output")
)
}
\seealso{
\code{\link[=IRFinder]{IRFinder()}}, \code{\link[=MakeSE]{MakeSE()}}
}
