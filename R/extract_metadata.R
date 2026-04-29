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
#' @param delim Single string used as the field delimiter. Default \code{"|"} for
#'   GISAID-style headers. Treated as a fixed string, not a regular expression.
#'
#' @return A \code{tibble} with \code{length(names)} columns and one row per
#'   sequence. Columns are typed as \code{character}; missing fields are
#'   \code{NA_character_}.
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
extract_metadata <- function(biostring, names, delim = "|") {

  # ── Input validation ─────────────────────────────────────────────────────────

  # biostring must be a Biostrings XStringSet (DNAStringSet, AAStringSet, etc.)
  if (!methods::is(biostring, "XStringSet")) {
    stop("biostring must be an XStringSet (e.g. DNAStringSet)")
  }

  # names must be a non-empty character vector — controls output column count
  if (!is.character(names) || length(names) == 0L) {
    stop("names must be a non-empty character vector")
  }

  # NA values in names would produce silently-unnamed columns in the output tibble
  if (anyNA(names)) {
    stop("names must not contain NA values")
  }

  # Duplicate column names would produce an ambiguous tibble
  dupe_names <- names[duplicated(names)]
  if (length(dupe_names) > 0L) {
    stop("names contains duplicate column names: ", paste(dupe_names, collapse = ", "))
  }

  # delim must be a single non-empty string — passed to stringr::fixed()
  if (!is.character(delim) || length(delim) != 1L || is.na(delim) || !nzchar(delim)) {
    stop("delim must be a single non-empty string")
  }

  # ── Core logic ───────────────────────────────────────────────────────────────

  # Placeholder — implemented in Task 2
  NULL
}
