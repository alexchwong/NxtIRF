# wrappers to R/C++

.run_IRFinder = function(
        reference_path = "./Reference", 
        bamfiles = "Unsorted.bam", 
        output_files = "./Sample",
        max_threads = max(parallel::detectCores() - 2, 1),
        Use_OpenMP = TRUE,
        run_featureCounts = FALSE,
        overwrite_IRFinder_output = FALSE,
        verbose = TRUE
    ) {
    .validate_reference(reference_path)
    s_bam = normalizePath(bamfiles)
    s_ref = normalizePath(reference_path)
    
    .irfinder_validate_args(s_bam, s_ref, max_threads, output_files)
    
    ref_file = normalizePath(file.path(s_ref, "IRFinder.ref.gz"))

    message("Running IRFinder ", appendLF = FALSE)
    n_threads = floor(max_threads)
    
    # OpenMP version currently causes C stack usage errors. Disable for now
    if(Has_OpenMP() > 0 & Use_OpenMP) {
        # n_threads = min(n_threads, length(s_bam))
        IRF_main_multi(ref_file, s_bam, output_files, n_threads, verbose)
    } else {
        # Use BiocParallel
        n_rounds = ceiling(length(s_bam) / floor(max_threads))
        n_threads = ceiling(length(s_bam) / n_rounds)

        BPPARAM = BiocParallel::bpparam()
        if(Sys.info()["sysname"] == "Windows") {
            BPPARAM_mod = BiocParallel::SnowParam(n_threads)
            message(paste("Using SnowParam", BPPARAM_mod$workers, "threads"))
        } else {
            BPPARAM_mod = BiocParallel::MulticoreParam(n_threads)
            message(paste("Using MulticoreParam", BPPARAM_mod$workers, 
                "threads"))
        }

        row_starts = seq(1, by = n_threads, length.out = n_rounds)
        for(i in seq_len(n_rounds)) {
            selected_rows_subset = seq(row_starts[i], 
                min(length(s_bam), row_starts[i] + n_threads - 1)
            )
            BiocParallel::bplapply(selected_rows_subset,
                function(i, s_bam, reference_file, output_files, verbose, overwrite) {
                    .irfinder_run_single(s_bam[i], reference_file, output_files[i], 
                        verbose, overwrite)
                }, 
                s_bam = s_bam,
                reference_file = ref_file,
                output_files = output_files,
                verbose = verbose,
                overwrite = overwrite_IRFinder_output,
                BPPARAM = BPPARAM_mod
            )
        }
    }
    if(run_featureCounts == TRUE) {
        .irfinder_run_featureCounts(reference_path, output_files, 
            s_bam, n_threads)
    }
}

.irfinder_run_single <- function(bam, ref, out, verbose, overwrite, n_threads = 1) {
    file_gz = paste0(out, ".txt.gz")
    file_cov = paste0(out, ".cov")
    bam_short = file.path(basename(dirname(bam)), basename(bam))
    if(overwrite ||
        !(file.exists(file_gz) | file.exists(file_cov))) {
        ret = IRF_main(bam, ref, out, verbose, n_threads)
        # Check IRFinder returns all files successfully
        if(ret != 0) {
            .log(paste(
                "IRFinder exited with errors, see error messages above"))
        } else if(!file.exists(file_gz)) {
            .log(paste(
                "IRFinder failed to produce", file_gz))
        } else if(!file.exists(file_cov)) {
            .log(paste(
                "IRFinder failed to produce", file_cov))
        } else {
            message(paste("IRFinder processed", bam_short))
        }
    } else {
        message(paste("IRFinder output for", 
            bam_short, "already exists, skipping..."))
    }
}

.irfinder_run_featureCounts <- function(reference_path, output_files, 
        s_bam, n_threads) {
    NxtIRF.CheckPackageInstalled("Rsubread", "2.4.0")
    gtf_file <- Get_GTF_file(reference_path)
    
    # determine paired_ness, strandedness, assume all BAMS are the same
    data.list = get_multi_DT_from_gz(
        normalizePath(paste0(output_files[1], ".txt.gz")), 
        c("BAM", "Directionality")
    )
    stats = data.list$BAM
    direct = data.list$Directionality

    paired = (stats$Value[3] == 0 & stats$Value[4] > 0) || 
        (stats$Value[3] > 0 && stats$Value[4] / stats$Value[3] / 1000)
    strand = direct$Value[9]
    if(strand == -1) strand = 2
    
    res = Rsubread::featureCounts(
        s_bam,
        annot.ext = gtf_file,
        isGTFAnnotationFile = TRUE,
        strandSpecific = strand,
        isPairedEnd = paired,
        requireBothEndsMapped = paired,
        nthreads = n_threads
    )
    # Append to existing main.FC.Rds if exists:
    if(file.exists(file.path(dirname(output_files[1]), "main.FC.Rds"))) {
        res.old = readRDS(
            file.path(dirname(output_files[1]), "main.FC.Rds"))

        # Check md5 of annotation to show same reference was used
        anno.old = res.old$annotation[, 
            c("GeneID", "Chr", "Start", "End", "Strand")]
        anno.new = res$annotation[, 
            c("GeneID", "Chr", "Start", "End", "Strand")]
        # md5.old = with(res.old$annotation, openssl::md5(paste(
            # GeneID, Chr, Start, End, Strand, collapse=" ")))
        # md5 = with(res$annotation, openssl::md5(paste(
            # GeneID, Chr, Start, End, Strand, collapse=" ")))
        md5.old.stat = openssl::md5(paste(res.old$stat$Status, collapse=" "))
        md5.stat = openssl::md5(paste(res$stat$Status, collapse=" "))
        if(identical(anno.old, anno.new) & md5.stat == md5.old.stat) {
            new_samples = res$targets[!(res$targets %in% res.old$targets)]
            res$targets = c(res.old$targets, new_samples)
            res$stat = cbind(res.old$stat, res$stat[,new_samples])        
            res$counts = cbind(res.old$counts, res$counts[,new_samples])
        }
    }   
    if(all(c("counts", "annotation", "targets", "stat") %in% names(res))) {
        saveRDS(res, file.path(dirname(output_files[1]), "main.FC.Rds"))
    }
    message(paste("featureCounts ran succesfully; saved to",
        file.path(dirname(output_files[1]), "main.FC.Rds")))
}


.irfinder_validate_args <- function(s_bam, s_ref, max_threads, output_files) {
    if(max_threads != 1 && max_threads > parallel::detectCores()) {
        .log(paste("In .run_IRFinder(), ",
            max_threads, " threads is not allowed for this system"))
    }
    if(!all(file.exists(s_bam))) {
        .log(paste("In .run_IRFinder(), ",
            paste(unique(s_bam[!file.exists(s_bam)]),
                collapse = ""),
            " - files not found"))
    }    

    if(!all(dir.exists(dirname(output_files)))) {
        .log(paste("In .run_IRFinder(), ",
            paste(unique(dirname(
                    output_files[!dir.exists(dirname(output_files))])),
                collapse = ""),
            " - directories not found"))
    }

    if(!(length(s_bam) == length(output_files))) {
        .log(paste("In .run_IRFinder(), ",
            "Number of output files and bam files must be the same"))
    }
    return(TRUE)
}

run_IRFinder_GenerateMapReads = function(genome.fa = "", out.fa, 
    read_len = 70, read_stride = 10, error_pos = 35) {
    return(
        IRF_GenerateMappabilityReads(normalizePath(genome.fa), 
            file.path(normalizePath(dirname(out.fa)), basename(out.fa)),
            read_len = read_len, 
            read_stride = read_stride, 
            error_pos = error_pos)
    )
}

run_IRFinder_MapExclusionRegions = function(bamfile = "", output_file, 
        threshold = 4, includeCov = FALSE, n_threads = 1) {
    s_bam = normalizePath(bamfile)
    if(!file.exists(s_bam)) {
        .log(paste("In run_IRFinder_MapExclusionRegions(),",
            s_bam, "does not exist"))
    }
    return(
        IRF_GenerateMappabilityRegions(s_bam, 
            output_file,
            threshold = threshold,
            includeCov = includeCov,
            verbose = TRUE, n_threads
        )
    )
}

run_Gunzip = function(infile = "", outfile) {
    file_to_read = normalizePath(infile)
    if(!file.exists(file_to_read)) {
        .log(paste("In run_Gunzip(),",
            file_to_read, "does not exist"))
    }
    if(!dir.exists(dirname(outfile))) {
        .log(paste("In run_Gunzip(),",
            dirname(outfile), "does not exist"))
    }
    IRF_gunzip(file_to_read, outfile)
}

get_multi_DT_from_gz = function(infile = "", 
        block_headers = c("Header1", "Header2")) {
    file_to_read = normalizePath(infile)
    if(!file.exists(file_to_read)) {
        .log(paste("In get_multi_DT_from_gz(),",
            file_to_read, "does not exist"))
    }
    df.list = IRF_gunzip_DF(file_to_read, block_headers)
    for(i in seq_len(length(df.list))) {
        for(j in seq_len(length(df.list[[i]]))) {
            suppressWarnings({
                if(all(df.list[[i]][[j]] == "NA" | 
                        !is.na(as.numeric(df.list[[i]][[j]])))) {
                    df.list[[i]][[j]] = as.numeric(df.list[[i]][[j]])
                }
            })
        }
        df.list[[i]] = as.data.table(df.list[[i]])
    }
    return(df.list)
}

#' Calls NxtIRF's C++ function to retrieve coverage
#'
#' This function returns an RLE or RLEList containing coverage data from the
#' given COV file
#' @param file The file name of the COV file
#' @param seqname Either blank, or a character string denoting the chromosome 
#'  name
#' @param start The 0-based start coordinate 
#' @param end The 0-based end coordinate
#' @param strand Either "*", "+", or "-"
#' @return If seqname is left as "", returns an RLEList of the whole BAM file.
#'   If seqname and coordinates are given, returns an RLE containing the
#'   chromosome coordinate. Coordinates outside the given range will be set to 0
#' @examples
#' se <- NxtIRF_example_NxtSE()
#' 
#' cov_file <- covfile(se)[1]
#'
#' cov <- GetCoverage(cov_file, seqname = "chrZ", 
#'   start = 10000, end = 20000,
#'   strand = "*"
#' )
#' @export
GetCoverage <- function(file, seqname = "", start = 0, end = 0, 
        strand = c("*", "+", "-")) {
    strand = match.arg(strand)
    if(!(strand %in% c("*", "+", "-"))) {
        .log(paste("In GetCoverage(),",
            "Invalid strand. '*', '+' or '-'"))
    }
    strand_int = ifelse(strand == "*", 2, 
        ifelse(strand == "+", 1, 0))
        
    if(!is.numeric(start) || !is.numeric(end) || 
            (as.numeric(start) > as.numeric(end))) {
        .log(paste("In GetCoverage(),",
            "Null or negative regions not allowed"))
    }
    if(seqname == "") {
        raw_list = IRF_RLEList_From_Cov(normalizePath(file), strand_int)
        final_list = list()
        if(length(raw_list) > 0) {
            for(i in seq_len(length(raw_list))) {
                final_list[[i]] = S4Vectors::Rle(
                    raw_list[[i]]$values, raw_list[[i]]$lengths
                )
            }
        } else {
            return(NULL)
        }
        final_RLE = as(final_list, "RleList")
        names(final_RLE) = names(raw_list)
        return(final_RLE)
    } else if(end == 0) {
        raw_RLE = IRF_RLE_From_Cov(
            normalizePath(file), as.character(seqname), 
            0,0, strand_int
        )
        final_RLE = S4Vectors::Rle(raw_RLE$values, raw_RLE$lengths)
    } else {
        raw_RLE = IRF_RLE_From_Cov(
            normalizePath(file), as.character(seqname), 
            round(as.numeric(start)), round(as.numeric(end)), 
            strand_int
        )
        final_RLE = S4Vectors::Rle(raw_RLE$values, raw_RLE$lengths)
    }
}

#' Validates the given file as a valid COV file
#' @param coverage_files A vector containing the file names of files to be
#'   checked
#' @return `TRUE` if all files are valid COV files. `FALSE` otherwise
#' @examples
#' se = NxtIRF_example_NxtSE()
#'
#' cov_files = covfile(se)
#'
#' IsCOV(cov_files) # returns true if these are true COV files
#' @md
#' @export
IsCOV = function(coverage_files) {
    for(i in coverage_files) {
        if(file.exists(i) && IRF_Check_Cov(normalizePath(i))) {
            # do nothing
        } else {
            return(FALSE)
        }
    }    
    return(TRUE)
}