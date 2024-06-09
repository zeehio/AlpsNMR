#' @title Export data for the spectral quantification library ASICS  
#' @description
#' A simple helper function for mangling data in the right format for the 
#' spectral quantification library ASICS. 
#' 
#' @param samples An `nmr_dataset_family` 1D object 
#' @param ... Additional arguments are passed 
#'            directly to `ASICS::createSpectra`, which (in theory) provide an 
#'            opportunity to use distinct normalisation methods. 
#' 
#' @return An `ASICS::Spectra` object 
#' @examples
#' # forAsics <- alps_asics(dataset)
#' # ASICS(forAsics)
#' 
#' @export 

alps_asics <- function(samples, ...) {
  forAsics <- t(nmr_data(samples))
  forAsics <- ASICS::createSpectra(forAsics, ...)
  forAsics@sample.name <-samples$metadata$external$NMRExperiment
  return(forAsics)
}


