test_that("Test extendr wrappers using `matrix`", {
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

test_that("Test extendr wrappers using `igraph` with `sparse = FALSE`", {
    # Test from the paper example 1.
    g <- igraph::graph_from_adjacency_matrix(
        matrix(c(
            # X1, X2, Y1, Y2, Y3
            FALSE, TRUE, TRUE, TRUE, TRUE, # X1
            FALSE, FALSE, TRUE, TRUE, TRUE, # X2
            FALSE, FALSE, FALSE, FALSE, FALSE, # Y1
            FALSE, FALSE, FALSE, FALSE, FALSE, # Y2
            FALSE, FALSE, FALSE, FALSE, FALSE # Y3
        ), nrow = 5, ncol = 5),
        mode = "directed"
    )

    print(g)

    h_1 <- igraph::graph_from_adjacency_matrix(
        matrix(c(
            # X1, X2, Y1, Y2, Y3
            FALSE, TRUE, TRUE, TRUE, TRUE, # X1
            FALSE, FALSE, TRUE, TRUE, TRUE, # X2
            FALSE, FALSE, FALSE, TRUE, FALSE, # Y1
            FALSE, FALSE, FALSE, FALSE, FALSE, # Y2
            FALSE, FALSE, FALSE, FALSE, FALSE # Y3
        ), nrow = 5, ncol = 5),
        mode = "directed"
    )

    h_2 <- igraph::graph_from_adjacency_matrix(
        matrix(c(
            # X1, X2, Y1, Y2, Y3
            FALSE, FALSE, TRUE, TRUE, TRUE, # X1
            TRUE, FALSE, TRUE, TRUE, TRUE, # X2
            FALSE, FALSE, FALSE, FALSE, FALSE, # Y1
            FALSE, FALSE, FALSE, FALSE, FALSE, # Y2
            FALSE, FALSE, FALSE, FALSE, FALSE # Y3
        ), nrow = 5, ncol = 5),
        mode = "directed"
    )

    expect_equal(
        sid(
            igraph::as_adjacency_matrix(g, sparse = FALSE),
            igraph::as_adjacency_matrix(h_1, sparse = FALSE)
        ),
        0
    )
    expect_equal(
        sid(
            igraph::as_adjacency_matrix(g, sparse = FALSE),
            igraph::as_adjacency_matrix(h_2, sparse = FALSE)
        ),
        8
    )
})

test_that("Test extendr wrappers using `igraph` with `sparse = TRUE`", {
    # Test from the paper example 1.
    g <- igraph::graph_from_adjacency_matrix(
        matrix(c(
            # X1, X2, Y1, Y2, Y3
            FALSE, TRUE, TRUE, TRUE, TRUE, # X1
            FALSE, FALSE, TRUE, TRUE, TRUE, # X2
            FALSE, FALSE, FALSE, FALSE, FALSE, # Y1
            FALSE, FALSE, FALSE, FALSE, FALSE, # Y2
            FALSE, FALSE, FALSE, FALSE, FALSE # Y3
        ), nrow = 5, ncol = 5),
        mode = "directed"
    )

    print(g)

    h_1 <- igraph::graph_from_adjacency_matrix(
        matrix(c(
            # X1, X2, Y1, Y2, Y3
            FALSE, TRUE, TRUE, TRUE, TRUE, # X1
            FALSE, FALSE, TRUE, TRUE, TRUE, # X2
            FALSE, FALSE, FALSE, TRUE, FALSE, # Y1
            FALSE, FALSE, FALSE, FALSE, FALSE, # Y2
            FALSE, FALSE, FALSE, FALSE, FALSE # Y3
        ), nrow = 5, ncol = 5),
        mode = "directed"
    )

    h_2 <- igraph::graph_from_adjacency_matrix(
        matrix(c(
            # X1, X2, Y1, Y2, Y3
            FALSE, FALSE, TRUE, TRUE, TRUE, # X1
            TRUE, FALSE, TRUE, TRUE, TRUE, # X2
            FALSE, FALSE, FALSE, FALSE, FALSE, # Y1
            FALSE, FALSE, FALSE, FALSE, FALSE, # Y2
            FALSE, FALSE, FALSE, FALSE, FALSE # Y3
        ), nrow = 5, ncol = 5),
        mode = "directed"
    )

    expect_equal(
        sid(
            igraph::as_adjacency_matrix(g),
            igraph::as_adjacency_matrix(h_1)
        ),
        0
    )
    expect_equal(
        sid(
            igraph::as_adjacency_matrix(g),
            igraph::as_adjacency_matrix(h_2)
        ),
        8
    )
})
