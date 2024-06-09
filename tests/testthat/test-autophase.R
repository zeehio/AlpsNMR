test_that("nmr_dataset_autophase works", {
    skip_if_not_installed("NMRphasing")
    lorentzian <- function(x, x0, gamma, A) {
        A * (1 / (pi * gamma)) * ((gamma^2) / ((x - x0)^2 + gamma^2))
    }
    
    x <- seq(from=1, to=2, length.out = 300)
    y <- lorentzian(x, 1.3, 0.01, 1) + lorentzian(x, 1.6, 0.01, 1)
    dataset <- new_nmr_dataset(
        metadata = list(external = data.frame(NMRExperiment = "10")),
        data_fields = list(
            data_1r = list(y)
        ),
        axis = list(list(x))
    )
    expect_warning(
        nmr_autophase(dataset, method="NLS"),
        "all_components=TRUE"
    )
})
