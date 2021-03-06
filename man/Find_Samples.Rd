% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CollateData.R
\name{Find_Samples}
\alias{Find_Samples}
\alias{Find_FASTQ}
\alias{Find_Bams}
\alias{Find_IRFinder_Output}
\title{Convenience Function to (recursively) find all files in a folder.}
\usage{
Find_Samples(
  sample_path,
  suffix = ".txt.gz",
  suffix_type = "path",
  use_subdir = FALSE
)

Find_FASTQ(
  sample_path,
  paired = TRUE,
  fastq_suffix_type = c("fastq", "fq", "fastq.gz", "fq.gz"),
  ...
)

Find_Bams(sample_path, ...)

Find_IRFinder_Output(sample_path, ...)
}
\arguments{
\item{sample_path}{The path in which to recursively search for files
that match the given \code{suffix}}

\item{suffix}{A vector of or or more strings that specifies the file suffix
(e.g. '.bam' denotes BAM files, whereas ".txt.gz" denotes gzipped txt
files).}

\item{suffix_type}{A string vector that determines the column
names of the files retrieved by \code{suffix}. Must be the same length as
\code{suffix}}

\item{use_subdir}{Whether to assume the directory name containing
the found files denote the sample name. If \code{FALSE}, the base name
of the file is assumed to be the sample name. See below example.}

\item{paired}{Whether to expect single FASTQ files (of the format
"sample.fastq"), or
paired files (of the format "sample_1.fastq", "sample_2.fastq")}

\item{fastq_suffix_type}{The name of the FASTQ suffix. Options are:
"fastq", "fastq.gz", "fq", or "fq.gz"}

\item{...}{Additional parameters to be passed to Find_Samples}
}
\value{
A multi-column data frame with the first column containing
the sample name, and subsequent columns being the file paths with suffix
as determined by \code{suffix}.
}
\description{
Often, output files (whether it be raw sequencing files, aligned sequences)
in BAM files, or IRFinder output files, are stored in a single folder.
Sometimes, these are named according to the samples they represent. Other
times, they have generic names but are partitioned in sub-folders named
by sample names. This function (recursively) finds all files and extracts
sample names assuming the files are named by sample names
(\code{use_subdir = FALSE}). Alternately, the names can be derived from the
subfolders one level higher (\code{use_subdir = TRUE})
}
\section{Functions}{
\itemize{
\item \code{Find_FASTQ}: Returns all FASTQ files in a given folder

\item \code{Find_Bams}: Returns all BAM files in a given folder

\item \code{Find_IRFinder_Output}: Returns all IRFinder output files in a given folder
}}

\examples{
# Retrieve all BAM files in a given folder, named by sample names
bam_path = dirname(example_bams())[1]
df.bams = Find_Samples(sample_path = bam_path, 
  suffix = ".bam", suffix_type = "BAM", use_subdir = FALSE)
# equivalent to:
df.bams = Find_Bams(bam_path)

# Retrieve all IRFinder output files in a given folder, 
# named by sample names

expr = Find_IRFinder_Output(file.path(tempdir(), "IRFinder_output"))
# equivalent to:
expr = Find_Samples(file.path(tempdir(), "IRFinder_output"),
    c(".txt.gz", ".cov"), 
    c("irf_file", "cov_file")
)
}
