#' Extract metadata from FASTA sequence headers into a tibble
#'
#' Splits delimited metadata embedded in sequence header names (e.g. GISAID-style
#' `hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05`) into a tidy tibble with
#' one row per sequence and one column per field. Returns metadata only — pair
#' with \code{fasta_nest()} if you also need the sequences.
#'
#' @param biostring An \code{XStringSet} (e.g. \code{DNAStringSet}) whose
#'   \code{names()} contain delimited metadata fields.
#' @param names Character vector of column names for the output tibble. Length
#'   controls how many fields are extracted; extra header fields are silently
#'   dropped. The first element is conventionally the sequence ID or strain name.
#'   Default \code{NULL} extracts all fields from the widest header and names
#'   them \code{V1}, \code{V2}, \ldots
#' @param delim Single string used as the field delimiter. Default \code{"|"} for
#'   GISAID-style headers. Treated as a fixed string, not a regular expression.
#'
#' @return A \code{tibble} with a leading \code{taxa} column (the full original
#'   header string) followed by one column per extracted field (all fields when
#'   \code{names = NULL}, otherwise \code{length(names)} columns), and one row
#'   per sequence. All columns are \code{character}; missing fields are
#'   \code{NA_character_}. Note: empty fields produced by consecutive delimiters
#'   (e.g. \code{"WNV||CO"}) are also converted to \code{NA} without a warning,
#'   as they are indistinguishable from padding.
#'
#' @export
#'
#' @examples
#' library(Biostrings)
#' seqs <- DNAStringSet(c(
#'   "hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05" = "ATCG",
#'   "hCoV-19/USA/CA-CDPH/2021|EPI_ISL_999999|2021-06-15" = "GCTA"
#' ))
#' extract_metadata(seqs, names = c("strain", "accession", "date"), delim = "|")
extract_metadata <- function(biostring, names = NULL, delim = "|") {

  # ── Input validation ─────────────────────────────────────────────────────────

  # biostring must be a Biostrings XStringSet (DNAStringSet, AAStringSet, etc.)
  if (!methods::is(biostring, "XStringSet")) {
    stop("biostring must be an XStringSet (e.g. DNAStringSet)")
  }

  # names, when provided, must be a non-empty character vector without NAs or dupes
  if (!is.null(names)) {
    if (!is.character(names) || length(names) == 0L) {
      stop("names must be a non-empty character vector")
    }
    if (anyNA(names)) {
      stop("names must not contain NA values")
    }
    dupe_names <- names[duplicated(names)]
    if (length(dupe_names) > 0L) {
      stop("names contains duplicate column names: ", paste(dupe_names, collapse = ", "))
    }
  }

  # delim must be a single non-empty string — passed to stringr::fixed()
  if (!is.character(delim) || length(delim) != 1L || is.na(delim) || !nzchar(delim)) {
    stop("delim must be a single non-empty string")
  }

  # ── Core logic ───────────────────────────────────────────────────────────────

  # Guard against empty input — max() on integer(0) returns -Inf with a warning
  if (length(biostring) == 0L) {
    stop("biostring must contain at least one sequence")
  }

  # Pull raw header strings from the XStringSet names slot.
  # NOTE: called as base::names() because the `names` argument shadows base::names()
  headers  <- base::names(biostring)
  n_delims <- stringr::str_count(headers, stringr::fixed(delim))

  # Determine split width and effective output column names.
  # When names = NULL: extract all fields from the widest header, auto-name V1, V2, ...
  # When names is supplied: split wide enough to cover all headers OR all requested
  # fields, then subset — str_split_fixed puts remainder in the last column so we
  # must split at least as wide as the widest header before subsetting.
  if (is.null(names)) {
    split_n   <- max(n_delims) + 1L
    names_out <- paste0("V", seq_len(split_n))
  } else {
    split_n   <- max(max(n_delims) + 1L, length(names))
    names_out <- names
  }

  mat <- stringr::str_split_fixed(
    string  = headers,
    pattern = stringr::fixed(delim),
    n       = split_n
  )

  # Subset to the requested column count; irrelevant when names = NULL (keep all).
  mat <- mat[, seq_along(names_out), drop = FALSE]

  # Warn when named extraction requests more fields than a header contains.
  # Skipped for NULL names since we always extract exactly what is present.
  if (!is.null(names)) {
    short_idx <- which(n_delims < (length(names_out) - 1L))
    if (length(short_idx) > 0L) {
      warning(
        "The following sequence headers have fewer fields than `names` ",
        "— missing fields set to NA:\n",
        paste(headers[short_idx], collapse = "\n"),
        call. = FALSE
      )
    }
  }

  # Convert "" padding produced by str_split_fixed for short headers to NA
  mat[mat == ""] <- NA_character_

  result           <- tibble::as_tibble(mat, .name_repair = "minimal")
  colnames(result) <- names_out
  result           <- tibble::add_column(result, taxa = headers, .before = 1L)

  result
}
