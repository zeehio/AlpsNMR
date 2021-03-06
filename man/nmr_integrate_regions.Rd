% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nmr_integrate_regions.R
\name{nmr_integrate_regions}
\alias{nmr_integrate_regions}
\alias{nmr_integrate_regions.nmr_dataset_1D}
\title{Integrate regions}
\usage{
nmr_integrate_regions(samples, regions, ...)

\method{nmr_integrate_regions}{nmr_dataset_1D}(
  samples,
  regions,
  fix_baseline = TRUE,
  excluded_regions_as_zero = FALSE,
  set_negative_areas_to_zero = FALSE,
  ...
)
}
\arguments{
\item{samples}{A \link{nmr_dataset} object}

\item{regions}{A named list. Each element of the list is a region,
given as a named numeric vector of length two with the range
to integrate. The name of the region will be the name of the
column}

\item{...}{Keep for compatibility}

\item{fix_baseline}{A logical. If \code{TRUE} it removes the baseline. See details
below}

\item{excluded_regions_as_zero}{A logical. It determines the behaviour of the
integration when integrating regions that have been excluded. If \code{TRUE},
it will treat those regions as zero. If \code{FALSE} (the default) it will return
NA values.

If \code{fix_baseline} is \code{TRUE}, then the region boundaries are used to estimate
a baseline. The baseline is estimated "connecting the boundaries with a straight
line". Only when the spectrum is above the baseline the area is integrated
(negative contributions due to the baseline estimation are ignored).}

\item{set_negative_areas_to_zero}{A logical. Ignored if \code{fix_baseline} is \code{FALSE}.
When set to \code{TRUE} negative areas are set to zero.}
}
\value{
An \link{nmr_dataset_peak_table} object
}
\description{
Integrate given regions and return a data frame with them
}
\examples{
#Creating a dataset
dataset <- new_nmr_dataset_1D(ppm_axis = 1:10,
                              data_1r = matrix(sample(0:99,replace = TRUE), nrow = 10),
                              metadata = list(external = data.frame(NMRExperiment = c("10", 
                              "20", "30", "40", "50", "60", "70", "80", "90", "100"))))

# Integrating selected regions
peak_table_integration = nmr_integrate_regions(
                                   samples = dataset,
                                   regions = list(ppm = c(2,5)),
                                   fix_baseline = TRUE)

#Creating a dataset
dataset <- new_nmr_dataset_1D(ppm_axis = 1:10,
                              data_1r = matrix(sample(0:99,replace = TRUE), nrow = 10),
                              metadata = list(external = data.frame(NMRExperiment = c("10",
                               "20", "30", "40", "50", "60", "70", "80", "90", "100"))))

# Integrating selected regions
peak_table_integration = nmr_integrate_regions(
                                   samples = dataset,
                                   regions = list(ppm = c(2,5)),
                                   fix_baseline = TRUE)
        
}
\seealso{
Other peak detection functions: 
\code{\link{Pipelines}},
\code{\link{nmr_baseline_threshold}()},
\code{\link{nmr_identify_regions_blood}()},
\code{\link{nmr_identify_regions_cell}()},
\code{\link{nmr_identify_regions_urine}()},
\code{\link{regions_from_peak_table}()},
\code{\link{validate_nmr_dataset_peak_table}()}

Other peak integration functions: 
\code{\link{Pipelines}},
\code{\link{computes_peak_width_ppm}()},
\code{\link{nmr_identify_regions_blood}()},
\code{\link{nmr_identify_regions_cell}()},
\code{\link{nmr_identify_regions_urine}()},
\code{\link{validate_nmr_dataset_peak_table}()}

Other nmr_dataset_1D functions: 
\code{\link{[.nmr_dataset_1D}()},
\code{\link{computes_peak_width_ppm}()},
\code{\link{file_lister}()},
\code{\link{files_to_rDolphin}()},
\code{\link{format.nmr_dataset_1D}()},
\code{\link{is.nmr_dataset_1D}()},
\code{\link{load_and_save_functions}},
\code{\link{new_nmr_dataset_1D}()},
\code{\link{nmr_align_find_ref}()},
\code{\link{nmr_baseline_removal}()},
\code{\link{nmr_baseline_threshold}()},
\code{\link{nmr_exclude_region}()},
\code{\link{nmr_interpolate_1D}()},
\code{\link{nmr_meta_add}()},
\code{\link{nmr_meta_export}()},
\code{\link{nmr_meta_get_column}()},
\code{\link{nmr_meta_get}()},
\code{\link{nmr_normalize}()},
\code{\link{nmr_pca_build_model}()},
\code{\link{nmr_pca_outliers_filter}()},
\code{\link{nmr_pca_outliers_plot}()},
\code{\link{nmr_pca_outliers_robust}()},
\code{\link{nmr_pca_outliers}()},
\code{\link{nmr_ppm_resolution}()},
\code{\link{plot.nmr_dataset_1D}()},
\code{\link{plot_webgl}()},
\code{\link{print.nmr_dataset_1D}()},
\code{\link{rdCV_PLS_RF_ML}()},
\code{\link{rdCV_PLS_RF}()},
\code{\link{save_files_to_rDolphin}()},
\code{\link{to_ChemoSpec}()},
\code{\link{validate_nmr_dataset_peak_table}()},
\code{\link{validate_nmr_dataset}()}
}
\concept{nmr_dataset_1D functions}
\concept{peak detection functions}
\concept{peak integration functions}
