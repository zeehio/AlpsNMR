% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nmr_detect_peaks_align.R
\name{nmr_detect_peaks}
\alias{nmr_detect_peaks}
\title{Peak detection for NMR}
\usage{
nmr_detect_peaks(
  nmr_dataset,
  nDivRange_ppm = 0.1,
  scales = seq(1, 16, 2),
  baselineThresh = NULL,
  SNR.Th = 3,
  range_without_peaks = c(9.5, 10),
  fit_lorentzians = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{nmr_dataset}{An \link{nmr_dataset_1D}.}

\item{nDivRange_ppm}{Segment size, in ppms, to divide the spectra and search
for peaks.}

\item{scales}{The parameter of peakDetectionCWT function of MassSpecWavelet package, look it up in the original function.}

\item{baselineThresh}{All peaks with intensities below the thresholds are excluded. Either:
\itemize{
\item A numeric vector of length the number of samples. Each number is a threshold for that sample
\item A single number. All samples use this number as baseline threshold.
\item \code{NULL}. If that's the case, a default function is used (\code{\link[=nmr_baseline_threshold]{nmr_baseline_threshold()}}), which assumes
that there is no signal in the region 9.5-10 ppm.
}}

\item{SNR.Th}{The parameter of peakDetectionCWT function of MassSpecWavelet package, look it up in the original function. If you set -1, the function will itself re-compute this value.}

\item{range_without_peaks}{A numeric vector of length two with a region without peaks, only used when \code{baselineThresh = NULL}}

\item{fit_lorentzians}{If \code{TRUE}, fit a lorentzian to each detected peak, to infer its inflection points. For now disabled for backwards compatibility.}

\item{verbose}{Logical (\code{TRUE} or \code{FALSE}). Show informational messages, such as the estimated baseline}
}
\value{
A data frame with the NMRExperiment, the sample index, the position
in ppm and index and the peak intensity
}
\description{
The function detects peaks on an \link{nmr_dataset_1D} object, using
\link[speaq:detectSpecPeaks]{speaq::detectSpecPeaks}. \code{detectSpecPeaks} divides the whole spectra into
smaller segments and uses \link[MassSpecWavelet:peakDetectionCWT]{MassSpecWavelet::peakDetectionCWT} for peak
detection.
}
\details{
Optionally afterwards, the peak apex and the peak inflection points are used to efficiently
adjust a lorentzian to each peak, and compute the peak area and width, as well as
the error of the fit. These peak features can be used afterwards to reject false
detections.
}
\seealso{
\link{nmr_align} for peak alignment with the detected peak table

Peak_detection

Other peak detection functions: 
\code{\link{Pipelines}},
\code{\link{nmr_baseline_threshold}()},
\code{\link{nmr_detect_peaks_plot}()},
\code{\link{nmr_detect_peaks_plot_overview}()},
\code{\link{nmr_detect_peaks_tune_snr}()},
\code{\link{nmr_identify_regions_blood}()},
\code{\link{nmr_identify_regions_cell}()},
\code{\link{nmr_identify_regions_urine}()},
\code{\link{nmr_integrate_regions}()}
}
\concept{peak detection functions}
