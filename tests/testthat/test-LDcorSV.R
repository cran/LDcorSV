test_that("Biloci measures", {
    data(data.test)
    Geno <- data.test[[1]]
    V.WAIS <- data.test[[2]]
    S.2POP <- data.test[[3]]
    expect_equal(Measure.R2(Geno), 0.0954443)
    expect_equal(Measure.R2S(Geno, S.2POP), 0.1081003)
    expect_equal(round(Measure.R2V(Geno, V.WAIS), 8), 0.1056732, tolerance = 1e-6)
    expect_equal(as.numeric(round(Measure.R2VS(Geno, V.WAIS, S.2POP), 8)), 0.01578605, tolerance = 1e-6)
})

test_that("LD measures", {
    data(data.test)
    Geno <- data.test[[1]]
    V.WAIS <- data.test[[2]]
    S.2POP <- data.test[[3]]
    LD <- LD.Measures(Geno, V = V.WAIS, S = S.2POP, data = "G", supinfo = TRUE, na.presence = TRUE)
    expect_equal(LD[173, 3], 0.00736549685246069, tolerance = 1e-6)
    expect_equal(LD[173, 4], 0.0096994662565902, tolerance = 1e-6)
    expect_equal(LD[173, 5], 0.0115996465046284, tolerance = 1e-6)
    expect_equal(LD[173, 6], 0.0103672763202665, tolerance = 1e-6)
    expect_equal(LD[173, 7], 0.021978021978022, tolerance = 1e-6)
    expect_equal(LD[173, 8], 0.043956043956044, tolerance = 1e-6)
    expect_equal(LD[173, 9], 0, tolerance = 1e-6)
    expect_equal(LD[173, 10], 0.0934065934065934, tolerance = 1e-6)
    expect_equal(LD[173, 11], 0.120879120879121, tolerance = 1e-6)
    expect_equal(LD[173, 12], 0, tolerance = 1e-6)
})
