#' nmr_dataset (S3 class)
#'
#' An `nmr_dataset` represents a set of NMR samples.
#' It is defined as an S3 class, and it can be treated as a regular list.
#'
#' It currently has the following elements:
#'
#'
#' - `metadata`: A list of data frames. Each data frame contains metadata of
#' a given area (acquisition parameters, preprocessing parameters, general sample information...)
#'
#' - `axis`: A list with length equal to the dimensionality of the data.
#' For 1D spectra it is a list with a numeric vector
#'
#' - `data_*`: Data arrays with the actual spectra. The first index represents
#' the sample, the rest of the indices match the length of each `axis`.
#' Typically `data_1r` is a matrix with one sample on each row and the chemical
#' shifts in the columns.
#'
#' - `num_samples`: The number of samples in the dataset
#'
#' @name nmr_dataset
#' @family AlpsNMR dataset objects
#' @seealso [Functions to save and load these objects][load_and_save_functions]
#' @examples
#' metadata_1D <- list(external = data.frame(NMRExperiment = c("10", "20")))
#' # Sample 10 and Sample 20 can have different lengths (due to different setups)
#' data_fields_1D <- list(data_1r = list(runif(16), runif(32)))
#' # Each sample has its own axis list, with one element (because this example is 1D)
#' axis_1D <- list(list(1:16), list(1:32))
#' my_1D_data <- new_nmr_dataset(metadata_1D, data_fields_1D, axis_1D)
NULL


#' Read NMR samples
#'
#' These functions load samples from files and return a [nmr_dataset].
#'
#' @name nmr_read_samples
#' @param sample_names A character vector with file or directory names.
#' @param samples_dir A directory or directories that contain multiple samples
#' @param format Either "bruker" or "jdx"
#' @param metadata_only A logical, to load only metadata (default: `FALSE`)
#' @param pulse_sequence If it is set to a pulse sequence
#'                                             ("NOESY", "JRES", "CPMG"...) it will only load
#'                                             the samples that match that pulse sequence.
#' @inheritDotParams read_bruker_pdata
#' @seealso [read_bruker_pdata()]
#' @return a [nmr_dataset] object
NULL

#' @rdname nmr_read_samples
#' @family import/export functions
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#'
nmr_read_samples_dir <- function(samples_dir,
    format = "bruker",
    pulse_sequence = NULL,
    metadata_only = FALSE,
    ...) {
    samples_dir <- as.character(samples_dir)
    dirs_that_dont_exist <- !dir.exists(samples_dir)
    if (any(dirs_that_dont_exist)) {
        rlang::abort(c("These directories do not exist:", samples_dir[dirs_that_dont_exist]))
    }
    if (format == "bruker") {
        all_samples <-
            c(
                list.dirs(
                    path = samples_dir,
                    full.names = TRUE,
                    recursive = FALSE
                ),
                list.files(
                    path = samples_dir,
                    full.names = TRUE,
                    pattern = ".*zip$"
                )
            )
    } else if (format == "jdx") {
        all_samples <-
            list.files(
                path = samples_dir,
                full.names = TRUE,
                pattern = ".*jdx$"
            )
    } else {
        stop("Unsupported sample format: ", format)
    }
    all_samples <- stringr::str_sort(all_samples, numeric = TRUE)
    dataset <- nmr_read_samples(
        sample_names = all_samples,
        format = format,
        pulse_sequence = pulse_sequence,
        metadata_only = metadata_only,
        ...
    )
    return(dataset)
}


#' @rdname nmr_read_samples
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' zip_files <- fs::dir_ls(dir_to_demo_dataset, glob = "*.zip")
#' dataset <- nmr_read_samples(sample_names = zip_files)
#'
nmr_read_samples <- function(sample_names,
    format = "bruker",
    pulse_sequence = NULL,
    metadata_only = FALSE,
    ...) {
    nn <- names(sample_names)
    if (all(sample_names == nn)) {
        nn <- NULL
    }
    sample_names <- as.character(sample_names)

    if (is.null(nn)) {
        names(sample_names) <- create_sample_names(sample_names)
    } else {
        names(sample_names) <- nn
    }

    if (format == "bruker") {
        samples <- nmr_read_samples_bruker(
            sample_names = sample_names,
            metadata_only = metadata_only,
            pulse_sequence = pulse_sequence,
            ...
        )
    } else if (format == "jdx") {
        # otherwise the jdx format
        samples <- nmr_read_samples_jdx(
            sample_names = sample_names,
            metadata_only = metadata_only
        )
    } else {
        stop("Unsupported format")
    }
    return(samples)
}

create_sample_names <- function(x) {
    # 1. Separate the zip path: "/a/b/c.zip!/d/e" -> "/a/b/c.zip"
    # 1. Use the basename without the extension: "c"
    # 2. If there are repeated names, prepend the dirname: "b/c"
    # 3. If there are repeated names, use vctrs::name_repair()
    has_zip_path <- grepl("\\.zip!.*$", x)
    x_without_zip_path <- x
    x_without_zip_path[has_zip_path] <- gsub(
        pattern = "\\.zip!.*$",
        replacement = ".zip",
        x_without_zip_path[has_zip_path]
    )
    remove_extensions <- tools::file_path_sans_ext(x_without_zip_path)
    xnames <- basename(remove_extensions)
    if (anyDuplicated(xnames) == 0) {
        return(xnames)
    }
    prepended_dirnames <- basename(dirname(remove_extensions))
    xnames2 <- paste(prepended_dirnames, xnames, sep = "/")
    if (anyDuplicated(xnames2) == 0) {
        return(xnames2)
    }
    vctrs::vec_as_names(xnames2, repair = "unique")
}

nmr_read_sample_bruker <- function(sample_path, pulse_sequence = NULL, metadata_only = FALSE, ...) {
    is_zip <- grepl("\\.zip$", sample_path) || grepl("\\.zip!.*$", sample_path)
    if (!is_zip) {
        sampl_dir <- normalizePath(sample_path)
    } else {
        zip_fn <- gsub(
            pattern = "(.*\\.zip)!(.*)$",
            replacement = "\\1",
            sample_path
        )
        zip_fn <- normalizePath(zip_fn)
        zip_basename <- basename(tools::file_path_sans_ext(zip_fn))
        if (grepl("(.*\\.zip)!(.*)$", sample_path)) {
            zip_subdir <- gsub(
                pattern = "(.*\\.zip)!(.*)$",
                replacement = "\\2",
                sample_path
            )
        } else {
            # FIXME:
            # By default, we assume the zip file has one folder with the same
            # name. This may not always be true, we could do some smart detection
            # here
            zip_subdir <- zip_basename
        }

        sampl_temp_dir <- tempfile(pattern = paste0("nmr_sample_", zip_basename, "_"))
        utils::unzip(zip_fn, exdir = sampl_temp_dir)
        on.exit({unlink(sampl_temp_dir, recursive = TRUE)})
        sampl_dir <- file.path(sampl_temp_dir, zip_subdir)
    }
    # Ignore internal TopSpin directory used for sample processing
    if (basename(sampl_dir) == "98888") {
        return(NULL)
    }
    meta <- read_bruker_metadata(sampl_dir)
    if (is_zip) {
        meta$info$file_format <- "Zipped Bruker NMR directory"
    }
    meta$info$sample_path <- sample_path
    if (!is.null(pulse_sequence) &&
        toupper(meta$info$pulse_sequence) != toupper(pulse_sequence)) {
        return(NULL)
    }
    if (metadata_only) {
        pdata <- NULL
    } else {
        pdata <- read_bruker_pdata(sample_path = sampl_dir, ...)
    }
    output <- bruker_merge_meta_pdata(meta, pdata)
    output
}

nmr_read_samples_bruker <-
    function(sample_names,
    pulse_sequence = NULL,
    metadata_only = FALSE,
    ...) {
        if (length(sample_names) == 0) {
            stop("No samples to load")
        }
        list_of_samples <- BiocParallel::bplapply(
            X = sample_names,
            FUN = function(sampl, pulse_sequence, metadata_only, ...) {
                loaded_sample <- 
                    tryCatch(
                        {
                            nmr_read_sample_bruker(sampl, pulse_sequence = pulse_sequence, metadata_only = metadata_only, ...)
                        },
                        error = function(err) {
                            msg <- conditionMessage(err)
                            rlang::warn(
                                message = c(
                                    "Error loading a sample",
                                    "i" = glue::glue("The sample '{sampl}' failed to load"),
                                    "i" = glue::glue("The underlying error message is: {msg}")
                                ),
                                underlying_error = err
                            )
                            return(err)
                        }
                    )
                return(loaded_sample)
            },
            pulse_sequence = pulse_sequence,
            metadata_only = metadata_only,
            ...
        )
        
        # Remove samples that could not be loaded:
        any_error <- purrr::map_lgl(list_of_samples, function(s) inherits(s, "error"))
        list_of_errors <- list_of_samples[any_error]
        list_of_samples <- list_of_samples[!any_error]
        sample_names <- sample_names[!any_error]


        if (length(list_of_samples) == 0) {
            rlang::abort(
                message = c(
                    "No samples could be loaded",
                    "i" = "You can check the underlying error messages with rlang::last_error()$error_list"
                ),
                error_list = list_of_errors
            )
        }

        # merge the sample information:
        all_fields <-
            unique(do.call(c, lapply(list_of_samples, function(x) {
                names(x)
            })))

        axis_fields <- "axis"
        data_fields <-
            all_fields[grepl(pattern = "^data_.*", x = all_fields)]
        metadata_fields <-
            setdiff(all_fields, c(axis_fields, data_fields))

        sample_meta <- list()
        for (meta_field in metadata_fields) {
            sample_meta[[meta_field]] <-
                list_of_lists_to_tibble(purrr::map(list_of_samples, meta_field))
            if (ncol(sample_meta[[meta_field]]) > 0) {
                colnames(sample_meta[[meta_field]]) <-
                    paste(meta_field, colnames(sample_meta[[meta_field]]), sep = "_")
            }
        }

        if (!is.null(names(sample_names))) {
            nmr_experiment_col <- names(sample_names)
        } else {
            nmr_experiment_col <- sample_meta[["info"]][["info_NMRExperiment"]]
        }
        nmr_experiment_col <- vctrs::vec_as_names(nmr_experiment_col, repair = "unique")
        sample_meta <- purrr::map(
            sample_meta,
            function(x) {
                x %>%
                    dplyr::mutate(NMRExperiment = nmr_experiment_col) %>%
                    dplyr::select("NMRExperiment", dplyr::everything())
            }
        )
        sample_meta[["external"]] <- tibble::tibble(NMRExperiment = nmr_experiment_col)
        data_fields_full <- list()
        axis <- NULL
        if (!metadata_only) {
            for (data_field in data_fields) {
                data_fields_full[[data_field]] <- purrr::map(list_of_samples, data_field)
            }
            axis <- purrr::map(list_of_samples, "axis")
        }
        samples <- new_nmr_dataset(
            metadata = sample_meta,
            data_fields = data_fields_full,
            axis = axis
        )
        return(samples)
    }

# @rdname nmr_read_samples
nmr_read_samples_jdx <-
    function(sample_names, metadata_only = FALSE) {
        nn <- names(sample_names)
        sample_names <- normalizePath(sample_names, mustWork = FALSE)
        names(sample_names) <- nn
        raw_samples <- read_jdx(sample_names, metadata_only = metadata_only)
        # Assume 1-D
        if (!metadata_only) {
            block_with_data_per_sample <-
                vapply(
                    raw_samples,
                    FUN = function(sample) {
                        block_with_xydata <-
                            vapply(
                                sample$block,
                                function(block) {
                                    "XYDATA" %in% names(block)
                                }, logical(1)
                            )
                        if (sum(block_with_xydata) == 1) {
                            return(which(block_with_xydata))
                        } else {
                            return(-1)
                        }
                    },
                    numeric(1)
                )
            if (any(block_with_data_per_sample == -1)) {
                stop("Loading samples: ", sample_names[block_with_data_per_sample == -1], "failed")
            }
        }
        num_samples <- length(raw_samples)

        # Metadata:
        metadata <-
            dplyr::bind_rows(lapply(raw_samples, create_df_from_jdx_sample))
        metadata$file_name <- sample_names
        # Make a reasonable NMRExperiment:
        # 1. If it is given as the names of sample_names
        # 2. If it is found in the metadata
        # 3. Based on the filename
        if (!is.null(nn)) {
            if (anyDuplicated(nn) > 0) {
                rlang::abort("names of samples must be unique")
            }
            NMRExperiments <- nn
        } else if (!"NMRExperiment" %in% colnames(metadata)) {
            NMRExperiments <- basename(sample_names)
            if (any(duplicated(NMRExperiments))) {
                NMRExperiments <- sample_names
            }
            NMRExperiments <- vctrs::vec_as_names(NMRExperiments, repair = "unique")
        } else {
            NMRExperiments <- metadata$NMRExperiment
        }
        metadata$NMRExperiment <- NMRExperiments
        metadata <-
            dplyr::select(metadata, "NMRExperiment", dplyr::everything())
        metadata_external <- tibble::tibble(NMRExperiment = metadata$NMRExperiment)


        axis <- NULL
        data_fields <- list()
        if (!metadata_only) {
            data_fields[["data_1r"]] <-
                vector(mode = "list", length = num_samples)
            axis <- vector(mode = "list", length(num_samples))
            for (sample_idx in seq_along(raw_samples)) {
                xydata <-
                    raw_samples[[sample_idx]]$blocks[[block_with_data_per_sample[sample_idx]]][["XYDATA"]]
                data_fields[["data_1r"]][[sample_idx]] <- xydata$y
                axis[[sample_idx]] <- list(x = xydata$x)
            }
        }
        samples <-
            new_nmr_dataset(
                metadata = list(
                    external = metadata_external,
                    metadata = metadata
                ),
                data_fields = data_fields,
                axis = axis
            )
        return(samples)
    }


#' Object is of [nmr_dataset] class
#' @param x An object
#' @return `TRUE` if the object is an [nmr_dataset], `FALSE` otherwise
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' is(dataset)
#'
is.nmr_dataset <- function(x) {
    inherits(x, "nmr_dataset")
}


#' Extract parts of an nmr_dataset
#' @param x an [nmr_dataset] object
#' @param i indices of the samples to keep
#' @return an nmr_dataset with the extracted samples
#' @family subsetting functions
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset2 <- dataset[1:3] # get the first 3 samples
#'
`[.nmr_dataset` <- function(x, i) {
    output <- x
    output$metadata <- purrr::map(output$metadata, function(metad) {
        metad[i, , drop = FALSE]
    })
    data_fields <-
        names(unclass(output))[grepl(pattern = "^data_.*", x = names(unclass(output)))]

    output[["axis"]] <- output[["axis"]][i]
    for (data_field in data_fields) {
        output[[data_field]] <- output[[data_field]][i]
    }
    output$num_samples <- nrow(output$metadata[[1]])
    validate_nmr_dataset(output)
    return(output)
}


#' Print for nmr_dataset
#' @param x an [nmr_dataset] object
#' @param ... for future use
#' @family class helper functions
#' @return Print for nmr_dataset
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' print(dataset)
#'
print.nmr_dataset <- function(x, ...) {
    cat(format(x, ...), "\n")
    invisible(x)
}

#' Format for nmr_dataset
#' @param x an [nmr_dataset] object
#' @param ... for future use
#' @family class helper functions
#' @return Format for nmr_dataset
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' format(dataset)
#'
format.nmr_dataset <- function(x, ...) {
    paste0("An nmr_dataset (", x$num_samples, " samples)")
}

#' Validate nmr_dataset objects
#'
#' @param samples An nmr_dataset object
#' @family class helper functions
#' @export
#' @return Validate nmr_dataset objects
#' @name validate_nmr_dataset
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' validate_nmr_dataset(dataset)
#'
validate_nmr_dataset <- function(samples) {
    validate_nmr_dataset_family(samples)
    abort_if_not(
        inherits(samples, "nmr_dataset"),
        message = "Not an nmr_dataset object"
    )
    samples
}

#' Create an nmr_dataset object
#'
#' @param metadata A named list of data frames
#' @param data_fields A named list. Check the examples
#' @param axis A list. Check the examples
#' @family class helper functions
#' @name new_nmr_dataset
#' @return Create an nmr_dataset object
#' @export
#' @return Create an nmr_dataset object
#' @examples
#' #
#' metadata_1D <- list(external = data.frame(NMRExperiment = c("10", "20")))
#' # Sample 10 and Sample 20 can have different lengths (due to different setups)
#' data_fields_1D <- list(data_1r = list(runif(16), runif(32)))
#' # Each sample has its own axis list, with one element (because this example is 1D)
#' axis_1D <- list(list(1:16), list(1:32))
#' my_1D_data <- new_nmr_dataset(metadata_1D, data_fields_1D, axis_1D)
#'
#' # Example for 2D samples
#' metadata_2D <- list(external = data.frame(NMRExperiment = c("11", "21")))
#' data_fields_2D <- list(data_2rr = list(matrix(runif(16 * 3), nrow = 16, ncol = 3),
#'     runif(32 * 3),
#'     nrow = 32, ncol = 3
#' ))
#' # Each sample has its own axis list, with one element (because this example is 1D)
#' axis_2D <- list(list(1:16, 1:3), list(1:32, 1:3))
#' my_2D_data <- new_nmr_dataset(metadata_2D, data_fields_2D, axis_2D)
#'
new_nmr_dataset <- function(metadata, data_fields, axis) {
    samples <- list()
    samples[["metadata"]] <- metadata
    samples <- append(x = samples, values = data_fields)
    samples[["axis"]] <- axis
    samples[["num_samples"]] <- nrow(metadata[[1]])
    class(samples) <- c("nmr_dataset", "nmr_dataset_family")
    validate_nmr_dataset(samples)
    samples
}
