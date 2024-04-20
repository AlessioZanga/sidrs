test_that("Test extendr wrappers", {
    # Test from the paper example 1.
    g <- matrix(c(
        # X1, X2, Y1, Y2, Y3
        FALSE, TRUE, TRUE, TRUE, TRUE, # X1
        FALSE, FALSE, TRUE, TRUE, TRUE, # X2
        FALSE, FALSE, FALSE, FALSE, FALSE, # Y1
        FALSE, FALSE, FALSE, FALSE, FALSE, # Y2
        FALSE, FALSE, FALSE, FALSE, FALSE # Y3
    ), nrow = 5, ncol = 5)

    h_1 <- matrix(c(
        # X1, X2, Y1, Y2, Y3
        FALSE, TRUE, TRUE, TRUE, TRUE, # X1
        FALSE, FALSE, TRUE, TRUE, TRUE, # X2
        FALSE, FALSE, FALSE, TRUE, FALSE, # Y1
        FALSE, FALSE, FALSE, FALSE, FALSE, # Y2
        FALSE, FALSE, FALSE, FALSE, FALSE # Y3
    ), nrow = 5, ncol = 5)

    h_2 <- matrix(c(
        # X1, X2, Y1, Y2, Y3
        FALSE, FALSE, TRUE, TRUE, TRUE, # X1
        TRUE, FALSE, TRUE, TRUE, TRUE, # X2
        FALSE, FALSE, FALSE, FALSE, FALSE, # Y1
        FALSE, FALSE, FALSE, FALSE, FALSE, # Y2
        FALSE, FALSE, FALSE, FALSE, FALSE # Y3
    ), nrow = 5, ncol = 5)

    expect_equal(sid(g, h_1), 0)
    expect_equal(sid(g, h_2), 8)
})
