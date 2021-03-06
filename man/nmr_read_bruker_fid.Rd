% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bruker.R
\name{nmr_read_bruker_fid}
\alias{nmr_read_bruker_fid}
\title{Read Free Induction Decay file}
\usage{
nmr_read_bruker_fid(sample_name, endian = "little")
}
\arguments{
\item{sample_name}{A single sample name}

\item{endian}{Endianness of the fid file ("little" by default, use "big" if acqus$BYTORDA == 1)}
}
\value{
A numeric vector with the free induction decay values
}
\description{
Reads an FID file. This is a very simple function.
}
\examples{
fid <- nmr_read_bruker_fid("sample.fid")
}
\seealso{
Other import/export functions: 
\code{\link{Pipelines}},
\code{\link{files_to_rDolphin}()},
\code{\link{load_and_save_functions}},
\code{\link{nmr_data}()},
\code{\link{nmr_meta_export}()},
\code{\link{nmr_read_samples}()},
\code{\link{nmr_zip_bruker_samples}()},
\code{\link{save_files_to_rDolphin}()},
\code{\link{save_profiling_output}()},
\code{\link{to_ChemoSpec}()}
}
\concept{import/export functions}
