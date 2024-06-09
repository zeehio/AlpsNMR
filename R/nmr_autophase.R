#' @title Rephase 1D NMR data
#' @description
#' Use phasing algorithms to rephase data in the spectral domain.
#' 
#' This function may improve autophasing processing from instrument vendors. It
#' wraps the [NMRphasing::NMRphasing()] function, to automatically rephase spectra,
#' allowing you to choose from a number of algorithms
#' of which `NLS`, `MPC_DANM` and `SPC_DANM` are the most recent.
#' 
#' Rephasing should happen before any spectra interpolation.
#' 
#' Please use the `all_components = TRUE` when calling [nmr_read_samples()] in order
#' to load the complex spectra and fix NMR phasing correctly.
#' 
#' @param dataset An [nmr_dataset] object 
#' @param method The autophasing method. See [NMRphasing::NMRphasing()] for details. 
#' @param withBC `NMRphasing::NMRphasing` may perform a baseline correction using modified polynomial fitting. By default
#' AlpsNMR offers other baseline estimation methods and better visualization of its effect, so AlpsNMR by default
#' disables the baseline correction offered by NMRphasing.
#' @param ... Other parameters passed on to [NMRphasing::NMRphasing()].
#' @return A (hopefully better phased) [nmr_dataset] object, with updated real and imaginary parts. 
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset, all_components=TRUE)
#' dataset <- nmr_autophase(dataset, method = "NLS")
#' dataset <- nmr_interpolate_1D(dataset, axis = c(min = 1, max = 2, by = 0.002))
#' plot(dataset)
#' 
#' @export 
nmr_autophase <- function(dataset,
                          method = c("NLS", "MPC_DANM", "MPC_EMP", "SPC_DANM", "SPC_EMP", "SPC_AAM", "SPC_DSM"),
                          withBC = FALSE, ...) {
    require_pkgs("NMRphasing")
    if (!"data_1i" %in% names(c(dataset))) {
        cli::cli_warn(c(
            "!" = "nmr_autophase() performs better with access to the whole complex NMR spectra",
            "i" = "Please read the dataset using {.code all_components=TRUE}.",
            "i" = "See {.fun AlpsNMR::nmr_autophase} for a full example"
        ))
    }
    
    if (inherits(dataset, "nmr_dataset_1D")) {
        cli::cli_abort(c(
            "x" = "nmr_autophase() expects non-interpolated spectra",
            "i" = "Please use nmr_autophase() before calling nmr_interpolate_1D()",
            "i" = "See {.fun AlpsNMR::nmr_autophase} for a full example"
        ))
    }
    
    method <- match.arg(method)
    
    real_list_of_spectra <- dataset$data_1r
    imag_list_of_spectra <- dataset$data_1i
    absorptionOnly <- FALSE
    if (is.null(imag_list_of_spectra)) {
        imag_list_of_spectra <- vector(mode="list", length=dataset$num_samples)
    } else {
        any_imag_missing <- purrr::map_lgl(imag_list_of_spectra, is.null)
        any_imag_missing <- any_imag_missing[any_imag_missing]
        if (length(any_imag_missing) > 0) {
            if (length(any_imag_missing < 7)) {
                miss_sample_names <- paste0(names(any_imag_missing), collapse = ", ")
                msg <- "Samples without imaginary component: {miss_sample_names}"
            } else {
                miss_sample_names <- paste0(names(head(any_imag_missing, n=5)), collapse = ", ")
                msg <- "Samples without imaginary component: {miss_sample_names} and {length(any_imag_missing)-5} more"
            }
            cli::cli_warn(c(
                "!" = "{length(any_imag_missing)}/{dataset$num_samples} samples have a missing imaginary spectrum",
                "i" = msg,
                "i" = "Estimating autophase using only the absorption"
            ))
        }
    }
    
    real_imag_lists <- BiocParallel::bpmapply(
        FUN = function(real, imag, ...) {
            if (!is.null(imag)) {
                absorptionOnly <- FALSE
                to_phase <- complex(
                    real = real, 
                    imaginary = imag
                )
            } else {
                absorptionOnly <- TRUE
                to_phase <- real
            }
            phased <- NMRphasing::NMRphasing(to_phase, absorptionOnly = TRUE, ...)
            list(real = Re(phased), imag = Im(phased))
        },
        real_list_of_spectra, imag_list_of_spectra,
        MoreArgs = list(method = method, withBC = withBC, ...),
        SIMPLIFY = FALSE
    )
    dataset[["data_1r"]] <- purrr::map(real_imag_lists, "real")
    if ("data_1i" %in% c(dataset)) {
        dataset[["data_1i"]] <- purrr::map(real_imag_lists, "imag")
    }
    dataset
}
