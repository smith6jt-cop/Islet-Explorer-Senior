# Three defects an adversarial contract review found in the Excel path, pinned so they stay fixed.
#
# They share a shape: each is silent. The app does not fail to start, it produces a wrong number or
# crashes a tab nobody has open, which is why none had been noticed.

.workbook_domains <- function() {
  out <- list()
  for (cfg in c("apps/islet/config/islet_default.yml", "apps/islet/config/run20_default.yml")) {
    path <- tryCatch(.find_repo_file(cfg), error = function(e) NULL)
    if (!is.null(path)) out[[cfg]] <- path
  }
  out
}

# domain_data_file() returns a path relative to the APP's working directory (apps/islet), while
# test_dir() runs from tests/testthat. Re-root it, or every assertion below silently skips --
# which is how a test that guards nothing still shows green.
.domain_workbook <- function() {
  rel <- domain_data_file("master", "master_results.xlsx")
  if (is.null(rel)) return(NULL)
  if (file.exists(rel)) return(normalizePath(rel))
  # .find_repo_file already walks the ../ prefixes, so hand it the repo-relative tail.
  tryCatch(.find_repo_file(sub("^(\\.\\./)+", "", rel)), error = function(e) NULL)
}

test_that("every workbook carries Age and the AAb flags on Islet_Composition", {
  # get_donor_metadata() (data_loading.R:24-43) returns on the FIRST sheet carrying
  # `Case ID` + `Donor Status`, which is always Islet_Composition when it exists. AAb columns
  # placed only on Islet_Targets are then deleted by select(-any_of(aab_cols_all)) at :118 and
  # never re-attached, so the Plot AAb filter and the Statistics AAb cards go dark silently.
  for (cfg in names(.workbook_domains())) {
    ie_set_domain(ie_load_domain(.workbook_domains()[[cfg]]))
    path <- .domain_workbook()
    skip_if(is.null(path), paste("no workbook for", cfg))

    sheets <- readxl::excel_sheets(path)
    expect_true("Islet_Composition" %in% sheets, info = cfg)
    comp <- readxl::read_excel(path, sheet = "Islet_Composition", n_max = 50)
    expect_true(all(c("Case ID", "Donor Status") %in% names(comp)), info = cfg)
    expect_true(any(grepl("^AAb_", names(comp))),
                info = paste(cfg, "- AAb flags must be on Islet_Composition, not only Targets"))
  }
})

test_that("Islet_Targets cannot fan out the islet table", {
  # region_um2 must be constant across the structure classes of one (islet, type). When it is
  # not, the distinct() at data_loading.R:58-63 leaves one key per class and the join at :256
  # multiplies every banded islet -- measured at 11,587 -> 33,783 rows on a perturbed workbook,
  # silently triplicating every plot and every statistical test.
  for (cfg in names(.workbook_domains())) {
    ie_set_domain(ie_load_domain(.workbook_domains()[[cfg]]))
    path <- .domain_workbook()
    skip_if(is.null(path), paste("no workbook for", cfg))

    # guess_max: readxl guesses a column's type from the first 1,000 rows. The CODEX
    # Islet_Targets sheet has a column that is empty until row ~48,000 and then holds AAb
    # names, so the default guess (logical) emits one warning per row -- 27,600 of them,
    # which buries every real result in the summary reporter.
    targets <- readxl::read_excel(path, sheet = "Islet_Targets", guess_max = 1e6)
    # base R on purpose: the dplyr path here emitted tens of thousands of warnings under
    # test_dir() but none standalone, i.e. it depended on package load order. A check that is
    # noisy in one harness and quiet in another is a check nobody reads.
    key <- paste(targets[["Case ID"]], targets[["region"]], targets[["type"]], sep = "\u001f")
    spread <- tapply(targets[["region_um2"]], key, function(v) length(unique(v)))
    varying <- sum(spread > 1, na.rm = TRUE)
    expect_equal(varying, 0L,
                 info = paste(cfg, "-", varying, "(islet, type) groups disagree on region_um2"))
  }
})

test_that("an age column reaches the donor cohort table on every data path", {
  # mod_inflammation_server.R:73 sorts on d$age. Absent, as.numeric(NULL) is numeric(0) and base
  # order() throws "argument lengths differ" -- and the Excel path supplied no age at all, so this
  # crashed for every Excel-backed dataset. The tab is hidden by default, which is why it went
  # unseen rather than unnoticed.
  for (cfg in names(.workbook_domains())) {
    ie_set_domain(ie_load_domain(.workbook_domains()[[cfg]]))
    path <- .domain_workbook()
    skip_if(is.null(path), paste("no workbook for", cfg))

    comp <- readxl::read_excel(path, sheet = "Islet_Composition", n_max = 100)
    donors <- comp[!duplicated(comp[["Case ID"]]), ]
    expect_true("Age" %in% names(donors) || "age" %in% names(donors),
                info = paste(cfg, "- no age column reaches the donor cohort table"))
    # and the sort itself must not throw
    expect_error(donors[order(donors[["Donor Status"]],
                              suppressWarnings(as.numeric(donors[["Age"]]))), ], NA)
  }
})
