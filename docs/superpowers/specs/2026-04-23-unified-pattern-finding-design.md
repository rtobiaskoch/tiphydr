# Design: Unified Pattern Finding with Smart Translation

**Date:** 2026-04-23  
**Status:** Approved  
**Scope:** Extract generalizable functions from examples/, create public API for pattern finding in DNA/protein sequences

---

## Problem Statement

Examples contain three overlapping patterns:
- `find_aa_loc()`, `get_start_met()`, `find_nucleotide_patterns()` all loop through XStringSet, find patterns, return dataframes
- Logic is duplicated across functions
- Translation start position detection is hardcoded or manual
- No unified way to search for residues (DNA/AA) across a sequence set

Goal: Create one composable `find_residue()` function that handles both DNA and protein sequences, with smart translation for DNA-to-AA searches.

---

## Design

### Public API

#### **`find_residue(xstring_set, pattern, start_position = NULL, verbose = TRUE)`**

**Purpose:** Search for one or more DNA or amino acid patterns in an `XStringSet` (DNA or protein).

**Input:**
- `xstring_set`: `DNAStringSet` or `AAStringSet` (or `RNAStringSet`)
- `pattern`: character vector of patterns to search (e.g., `"ATG"` or `c("K", "M")`)
- `start_position`: integer; optional override for translation start (default: auto-detect)
- `verbose`: logical; if TRUE, messages the detected start position for DNA→AA translation

**Behavior:**
1. Detect sequence type (DNA vs protein) and pattern type
2. If DNA + AA pattern:
   - Call `detect_sequence_start()` to find modal ATG position
   - Call `translate_from_position()` to translate
   - Search in translated sequence
   - Message user of detected position
3. If DNA + DNA pattern or protein + AA pattern:
   - Search directly
4. Return long-format tibble: one row per match

**Output:** tibble with columns:
- `sequence_name` (character): from `names(xstring_set)`
- `location` (integer): position of match in the searched sequence (DNA positions for DNA patterns, AA positions for AA patterns)
- `pattern` (character): the pattern searched for

**Rationale:**
- Long format makes it easy to filter, group, and compose with tidyverse
- One row per match (not nested) simplifies downstream analysis
- Auto-translation with modal detection handles real-world data: missing start codons, sequencing errors, incomplete sequences

---

### Internal Helpers

#### **`detect_sequence_start(xstring_set, pattern = "ATG")`**

**Purpose:** Find the modal (most common) start position for a pattern across all sequences in a DNA StringSet.

**Input:**
- `xstring_set`: `DNAStringSet`
- `pattern`: character; pattern to search (default: "ATG")

**Output:** list with elements:
- `position` (integer): the modal start position
- `table` (table): count of each position (for transparency)

**Behavior:**
- Use `Biostrings::matchPattern()` to find all matches of pattern in each sequence
- For each sequence, record the start position of the first match
- Ignore sequences with no match (NAs excluded from mode calculation)
- Return mode and message showing: detected position, count of sequences with that position, count of sequences missing the pattern

**Rationale:**
- Robust to real-world data: incomplete sequences, sequencing errors that truncate start codons
- Transparent: messages show why a position was chosen
- Reusable: can be called directly by users or by `define_lineage()`

---

#### **`translate_from_position(xstring_set, start_position)`**

**Purpose:** Translate a DNA StringSet starting from a specific position.

**Input:**
- `xstring_set`: `DNAStringSet`
- `start_position`: integer; 1-indexed position to start translation

**Output:** `AAStringSet` with the same names as input

**Behavior:**
- Extract subsequence from `start_position` to end
- Use `Biostrings::translate()` with `if.fuzzy.codon = "X"`
- Return translated sequences

**Rationale:**
- Pure function: no side effects, testable in isolation
- Reusable: can be called standalone for any translation task
- Consistent with RSE principle: small, composable units

---

### Internal Helpers (Private)

#### **`is_aa_pattern(pattern)`**

**Purpose:** Detect whether a pattern is amino acids or nucleotides.

**Input:**
- `pattern`: character string

**Output:** logical; TRUE if pattern contains only standard amino acid codes (ACDEFGHIKLMNPQRSTVWY), FALSE if nucleotides (ACGTN) or ambiguous

**Rationale:**
- Enables `find_residue()` to auto-detect pattern type
- Edge case: some patterns are ambiguous (e.g., "A" is both AA and nucleotide) — prioritize AA

---

#### **`search_pattern_in_xstring(xstring, pattern)`**

**Purpose:** Generic search for a single pattern in a single XString (DNA or protein).

**Input:**
- `xstring`: `XString` (DNAString, AAString, or similar)
- `pattern`: character; one pattern to search

**Output:** integer vector of start positions of matches (empty if no matches)

**Behavior:**
- Wrapper around `Biostrings::matchPattern()`
- Use `start()` to extract positions

**Rationale:**
- Isolates Biostrings API so `find_residue()` is clean and testable
- Can be vectorized with `purrr::map()` to search multiple patterns

---

## Updates to Existing Functions

### **`define_lineage()`**

Current implementation should use `detect_sequence_start()` instead of duplicating start-position-detection logic.

**Before:** Hardcoded or manual start position specification  
**After:** Call `detect_sequence_start(xstring_set, pattern = "ATG")` to auto-detect modal position, message user, allow override via parameter

**Rationale:**
- Single source of truth for start position detection
- Consistent behavior with `find_residue()`
- Easier to test and maintain

---

## Testing Strategy

### Unit Tests (each function in isolation)

**`detect_sequence_start()`:**
- Sequences with single ATG at same position → returns that position
- Sequences with ATG at different positions → returns mode
- Sequences with no ATG → returns mode of remaining sequences, messages count of missing
- Empty DNAStringSet → error handling
- Frame calculation edge cases

**`translate_from_position()`:**
- Translate from position 1 (full sequence)
- Translate from position 2, 3 (offset starts)
- Fuzzy codons handled correctly (X returned)
- Different DNAStringSet sizes

**`is_aa_pattern()`, `search_pattern_in_xstring()`:**
- Standard test cases for pattern detection and search

### Integration Tests

**`find_residue()`:**
- DNA + DNA pattern (e.g., "ATG", "TAA")
- Protein + AA pattern (e.g., "K", "M")
- DNA + AA pattern (should auto-translate):
  - All sequences have ATG at same position
  - Sequences have ATG at different positions (modal detection)
  - Some sequences missing ATG (robust handling)
- Multiple patterns in one call
- Verbose message output captures start position
- User override of start_position parameter
- Output format (tibble, long format, one row per match)

### Edge Cases

- Empty patterns → error or empty result?
- Pattern longer than sequence → should return nothing
- Case sensitivity (should be case-insensitive for AA codes)
- Ambiguous patterns (e.g., "A" that could be AA or nucleotide) — document priority

---

## Removed/Deprecated Examples

These examples in `examples/` will be replaced by composable calls to `find_residue()`:

- `find_aa_loc()` → `find_residue(xstring_set, pattern = "K")` after translation
- `get_start_met()` → `detect_sequence_start(xstring_set)`  
- `find_nucleotide_patterns()` → `find_residue(xstring_set, pattern = "ATG")`
- `get_stop_codons.R` → `find_residue(xstring_set, pattern = c("TAA", "TAG", "TGA"))`

Cleaner example scripts will demonstrate composition.

---

## File Structure

New files in `R/`:
- `find_residue.R` — main public function
- `detect_sequence_start.R` — internal helper
- `translate_from_position.R` — internal helper (could be public if useful)
- `is_aa_pattern.R`, `search_pattern_in_xstring.R` — private helpers (same file as main or separate)

Updated file:
- `define_lineage.R` — refactor to use `detect_sequence_start()`

New tests in `tests/testthat/`:
- `test-find_residue.R`
- `test-detect_sequence_start.R`
- `test-translate_from_position.R`

Example scripts in `examples/`:
- Rewrite to use `find_residue()` and new helpers
- Show common workflows (find stop codons, find start codon, find AA in translated sequence)

---

## Implementation Order

1. **Phase 1: Core helpers**
   - `translate_from_position()` (pure, testable, no dependencies on others)
   - `is_aa_pattern()` (pure, standalone)
   - `search_pattern_in_xstring()` (wrapper, testable)

2. **Phase 2: Start detection**
   - `detect_sequence_start()` (depends on Phase 1)
   - Unit tests

3. **Phase 3: Main function**
   - `find_residue()` (composes all helpers)
   - Integration tests

4. **Phase 4: Update existing**
   - Refactor `define_lineage()` to use `detect_sequence_start()`
   - Rewrite example scripts

---

## Notes on Design Decisions

**Why long format (one row per match)?**
- Easier to filter, group, summarize with tidyverse
- Composable with downstream analysis
- Avoids nested list-columns (simpler API)

**Why auto-detect start position?**
- Real-world sequences often have missing/shifted start codons (sequencing errors, incomplete reads)
- Modal approach is robust: tolerates some failures
- User can override when needed

**Why detect sequence type and pattern type?**
- Reduces user cognitive load
- Enables composition: `find_residue(fasta_read(...), "K")` just works
- Error on ambiguity (e.g., "A" defaults to AA, documented)

**Why return modal position instead of per-sequence position?**
- Simpler API: one start position per run
- Appropriate for typical use (homogeneous sequence set)
- Still respects variation (modal captures most common behavior)

---

## Open Questions / Future Extensions

1. **Frames:** Should `find_residue()` support frame-aware searching (e.g., find pattern only in frame 1)? Deferred — not needed by examples, can add later.
2. **Region restriction:** Should users be able to search only a subsequence (e.g., "codons 1-100")? Deferred.
3. **Make `translate_from_position()` public?** Yes — it's useful on its own. Include in exports.
4. **Performance:** For very large StringSets, consider vectorization with `purrr::map_df()` instead of loops. Deferred until profiling shows need.
