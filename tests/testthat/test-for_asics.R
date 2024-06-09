test_that("to_ASICS works", {
    skip_if_not_installed("ASICS")
    # Create a random spectra matrix
    nsamp <- 3
    npoints <- 300
    dummy_ppm_axis <- seq(from = 0.2, to = 10, length.out = npoints)
    dummy_spectra_matrix <- matrix(runif(nsamp * npoints), nrow = nsamp, ncol = npoints)
    metadata <- list(external = data.frame(
        NMRExperiment = paste0("Sample", seq_len(nsamp)),
        DummyClass = c("a", "b")[sample.int(2, nsamp, replace=TRUE)]
    ))
    dummy_nmr_dataset_1D <- new_nmr_dataset_1D(
        ppm_axis = dummy_ppm_axis,
        data_1r = dummy_spectra_matrix,
        metadata = metadata
    )
    spec_obj <- to_ASICS(dummy_nmr_dataset_1D)
    expect_true(class(spec_obj)[1] == "Spectra")
})