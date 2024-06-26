% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nmr_detect_peaks_align.R
\name{nmr_detect_peaks_tune_snr}
\alias{nmr_detect_peaks_tune_snr}
\title{Diagnose SNR threshold in peak detection}
\usage{
nmr_detect_peaks_tune_snr(
  ds,
  NMRExperiment = NULL,
  SNR_thresholds = seq(from = 2, to = 6, by = 0.1),
  ...
)
}
\arguments{
\item{ds}{An \link{nmr_dataset_1D} dataset}

\item{NMRExperiment}{A string with the single NMRExperiment used explore the SNR thresholds. If not given, use the first one.}

\item{SNR_thresholds}{A numeric vector with the SNR thresholds to explore}

\item{...}{
  Arguments passed on to \code{\link[=nmr_detect_peaks]{nmr_detect_peaks}}
  \describe{
    \item{\code{nmr_dataset}}{An \link{nmr_dataset_1D}.}
    \item{\code{nDivRange_ppm}}{Segment size, in ppms, to divide the spectra and search
for peaks.}
    \item{\code{baselineThresh}}{All peaks with intensities below the thresholds are excluded. Either:
\itemize{
\item A numeric vector of length the number of samples. Each number is a threshold for that sample
\item A single number. All samples use this number as baseline threshold.
\item \code{NULL}. If that's the case, a default function is used (\code{\link[=nmr_baseline_threshold]{nmr_baseline_threshold()}}), which assumes
that there is no signal in the region 9.5-10 ppm.
}}
    \item{\code{range_without_peaks}}{A numeric vector of length two with a region without peaks, only used when \code{baselineThresh = NULL}}
    \item{\code{fit_lorentzians}}{If \code{TRUE}, fit a lorentzian to each detected peak, to infer its inflection points. For now disabled for backwards compatibility.}
    \item{\code{verbose}}{Logical (\code{TRUE} or \code{FALSE}). Show informational messages, such as the estimated baseline}
    \item{\code{scales}}{The parameter of peakDetectionCWT function of MassSpecWavelet package, look it up in the original function.}
    \item{\code{SNR.Th}}{The parameter of peakDetectionCWT function of MassSpecWavelet package, look it up in the original function. If you set -1, the function will itself re-compute this value.}
  }}
}
\value{
A list with the following elements:
\itemize{
\item \code{peaks_detected}: A data frame with the columns from the \link{nmr_detect_peaks} output and an additional column
\code{SNR_threshold} with the threshold used on each row.
\item \code{num_peaks_per_region}: A summary of the \code{peaks_detected} table, with the number of peaks detected on
each chemical shift region
\item \code{plot_num_peaks_per_region}: A visual representation of \code{num_peaks_per_region}
\item \code{plot_spectrum_and_detections}: A visual representation of the spectrum and the peaks detected with each
SNR threshold. Use \link[plotly:ggplotly]{plotly::ggplotly} or \link{plot_interactive} on this to zoom and explore the results.
}
}
\description{
Diagnose SNR threshold in peak detection
}
\seealso{
nmr_detect_peaks

Other peak detection functions: 
\code{\link{Pipelines}},
\code{\link{nmr_baseline_threshold}()},
\code{\link{nmr_detect_peaks}()},
\code{\link{nmr_detect_peaks_plot}()},
\code{\link{nmr_detect_peaks_plot_overview}()},
\code{\link{nmr_identify_regions_blood}()},
\code{\link{nmr_identify_regions_cell}()},
\code{\link{nmr_identify_regions_urine}()},
\code{\link{nmr_integrate_regions}()}
}
\concept{peak detection functions}
