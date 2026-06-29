test_that("non-string query throws error", {
  expect_error(fetch_fasta_query(123, "out.fasta"), "single non-empty string")
})

test_that("bad filename throws error", {
  expect_error(fetch_fasta_query("WNV[Organism]", ""), "non-empty string")
})

test_that("unnamed filters throw error", {
  expect_error(
    fetch_fasta_query("WNV[Organism]", "out.fasta", filters = c("complete genome")),
    "named character vector"
  )
})

test_that("invalid retmax throws error", {
  expect_error(
    fetch_fasta_query("WNV[Organism]", "out.fasta", retmax = -1),
    "positive number"
  )
})

test_that("filters build correct AND query", {
  expect_equal(
    tiphydr:::build_entrez_query(
      "West Nile virus[Organism]",
      c(Title = "complete genome", PDAT = "2010:2024")
    ),
    "West Nile virus[Organism] AND complete genome[Title] AND 2010:2024[PDAT]"
  )
})

test_that("NULL filters returns query unchanged", {
  expect_equal(
    tiphydr:::build_entrez_query("West Nile virus[Organism]", NULL),
    "West Nile virus[Organism]"
  )
})

test_that("fetches sequences and writes file", {
  skip_if_offline()
  tmp <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp))

  result <- suppressMessages(
    fetch_fasta_query(
      query   = "West Nile virus[Organism]",
      filename = tmp,
      filters  = c(Title = "complete genome"),
      retmax   = 5
    )
  )

  expect_s4_class(result, "DNAStringSet")
  expect_gte(length(result), 1L)
  expect_true(file.exists(tmp))
})
