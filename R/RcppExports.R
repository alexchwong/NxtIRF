# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Has_OpenMP <- function() {
    .Call(`_NxtIRF_Has_OpenMP`)
}

IRF_Check_Cov <- function(s_in) {
    .Call(`_NxtIRF_IRF_Check_Cov`, s_in)
}

IRF_RLE_From_Cov <- function(s_in, seqname, start, end, strand) {
    .Call(`_NxtIRF_IRF_RLE_From_Cov`, s_in, seqname, start, end, strand)
}

IRF_RLEList_From_Cov <- function(s_in, strand) {
    .Call(`_NxtIRF_IRF_RLEList_From_Cov`, s_in, strand)
}

IRF_gunzip <- function(s_in, s_out) {
    .Call(`_NxtIRF_IRF_gunzip`, s_in, s_out)
}

IRF_gunzip_DF <- function(s_in, s_header_begin) {
    .Call(`_NxtIRF_IRF_gunzip_DF`, s_in, s_header_begin)
}

IRF_main <- function(bam_file, reference_file, output_file, verbose = TRUE, n_threads = 1L) {
    .Call(`_NxtIRF_IRF_main`, bam_file, reference_file, output_file, verbose, n_threads)
}

IRF_main_multi <- function(reference_file, bam_files, output_files, max_threads, verbose = TRUE) {
    .Call(`_NxtIRF_IRF_main_multi`, reference_file, bam_files, output_files, max_threads, verbose)
}

IRF_GenerateMappabilityReads <- function(genome_file, out_fa, read_len, read_stride, error_pos) {
    .Call(`_NxtIRF_IRF_GenerateMappabilityReads`, genome_file, out_fa, read_len, read_stride, error_pos)
}

IRF_GenerateMappabilityRegions <- function(bam_file, output_file, threshold, includeCov, verbose, n_threads = 1L) {
    .Call(`_NxtIRF_IRF_GenerateMappabilityRegions`, bam_file, output_file, threshold, includeCov, verbose, n_threads)
}

