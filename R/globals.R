globalVariables(c(":=","."))

buildref_version <- "0.99.0"

is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

is_valid <- function(x) {
    !missing(x) &&
    !is.null(x) && length(x) > 0 && 
        (isS4(x) || !is.na(x)) && 
        (!is.character(x) || (x != "" && x != "(none)"))
}

.log <- function(msg = "", 
        type = c("error", "warning", "silent", "message"),
        use_system_time = TRUE,
        ...
) {
    type = match.arg(type)
    if(use_system_time) {
        msg = paste(format(Sys.time(), "%b %d %X"), msg)
    }
    if(type == "error") {
        stop(msg, call. = FALSE)
    } else if(type == "warning") {
        warning(msg)
    } else {
        message(msg)
    }
}

.check_package_installed <- function(
        package = "DESeq2", version = "0.0.0",
        returntype = c("error", "warning", "silent")
) {
    res <- tryCatch(
        ifelse(packageVersion(package) >= version, TRUE, FALSE),
        error = function(e) FALSE)
    if (!res) {
        returntype <- match.arg(returntype)
        .log(paste(package, "version", version, "is not installed;",
            "and is required for this function"), type = returntype)
    }
    return(res)
}

.validate_threads <- function(n_threads, as_BPPARAM = TRUE, ...) {
    n_threads_to_use = as.numeric(n_threads)
    if(is.na(n_threads_to_use)) {
        .log("n_threads must be a numeric value")
    }
    if(n_threads_to_use > (parallel::detectCores() - 2) ) {
        n_threads_to_use = max(1, parallel::detectCores() - 2)
    }
    if(as_BPPARAM) {
        if(Sys.info()["sysname"] == "Windows") {
            BPPARAM_mod = BiocParallel::SnowParam(n_threads_to_use, ...)
            message(paste("Using SnowParam", BPPARAM_mod$workers, "threads"))
        } else {
            BPPARAM_mod = BiocParallel::MulticoreParam(n_threads_to_use, ...)
            message(paste("Using MulticoreParam", BPPARAM_mod$workers, "threads"))
        }
        return(BPPARAM_mod)
    } else {
        return(n_threads_to_use)
    }

}

# Semi-join a data.table. Equivalent to dplyr::semi_join(A, B, by = by)
.semi_join_DT <- function(A, B, by, nomatch = 0) {
    A[A[B, on = by, which = TRUE, nomatch = nomatch]]
}

# Converts data table to GRanges object, preserving info
.grDT <- function(DT, ...) {
    if (nrow(DT) == 0) return(GenomicRanges::GRanges(NULL))
    makeGRangesFromDataFrame(as.data.frame(DT), ...)
}

.grlGaps <- function(grl) {
    psetdiff(unlist(range(grl), use.names = TRUE), grl)
}

.make_path_relative <- function(files, relative_to) {
    if (!all(file.exists(files))) .log("Some files do not exist")
    if (!all(file.exists(relative_to))) .log("Some directories do not exist")
    if (length(relative_to) == 1) relative_to <- rep(relative_to, length(files))

    if (Sys.info()["sysname"] == "Windows") {
        files <- normalizePath(files, winslash = "/")
        relative_to <- normalizePath(relative_to, winslash = "/")
    } else {
        files <- normalizePath(files)
        relative_to <- normalizePath(relative_to)
    }
    out <- c()
    for (i in seq_len(length(files))) {
        f <- files[i]
        base <- relative_to[i]

        common <- sub("^([^|]*)[^|]*(?:\\|\\1[^|]*)$", "^\\1/?",
            paste0(base, "|", f))

        out <- c(out, paste0(gsub("[^/]+/?", "../", sub(common, "", base)),
            sub(common, "", f)))
    }
    return(out)
}

# Compatibility for running inside a shiny withProgress block
dash_progress <- function(message = "", total_items = 1, add_msg = FALSE) {
    if (total_items != round(total_items) | total_items < 1) {
        .log("dash_progress needs at least 1 item")
    }
    if (add_msg) {
        .log(message, "message")
    }
    has_shiny <- .check_package_installed(
        package = "shiny", returntype = "silent")
    if (has_shiny) {
        session <- shiny::getDefaultReactiveDomain()
        if (!is.null(session)) {
            shiny::incProgress(1 / total_items, message = message)
        }
    }
}

# Equivalent to shiny::withProgress with compatibility if shiny is missing.
dash_withProgress <- function(expr, min = 0, max = 1,
    value = min + (max - min) * 0.1,
    message = NULL, detail = NULL,
    # style = getShinyOption("progress.style", default = "notification"),
    # session = getDefaultReactiveDomain(),
    env = parent.frame(), quoted = FALSE
) {

    has_shiny <- .check_package_installed(
        package = "shiny", returntype = "silent")
    if (has_shiny) {
        session <- shiny::getDefaultReactiveDomain()
        if (!is.null(session)) {
            shiny::withProgress(expr = expr, min = min, max = max,
                value = value, message = message, detail = detail,
                env = env, quoted = quoted
            )
        } else {
            if (!quoted) expr <- substitute(expr)
            eval(expr, env)
        }
    } else {
        if (!quoted) expr <- substitute(expr)
        eval(expr, env)
    }
}
