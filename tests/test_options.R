library(testthat)
library(stringr)

source("lib/options.R")

# grab binding to this function so it can be overwritten in the tests
commandArgs <- NULL

test_that("initialize_parser works", {
  expect_true(TRUE)

  version <- "999"
  program <- "TEST123"
  description <- "testing parser"

  parser <- initialize_parser(version, program, description)

  test_that("parser prints version number when passed -v", {
    flags <- c("-v")

    expect_error(parser$parse_args(flags), version)
  })

  test_that("parser prints help message when passed -h", {
    flags <- c("-h")

    expect_error(parser$parse_args(flags), "METADATA DATA OUTPUT")
    expect_error(parser$parse_args(flags), program)
    expect_error(parser$parse_args(flags), description)
  })

  test_that("parser requires 3 positional arguments for the files", {
    not_enough_positional_args <- list(c(), c("m"), c("m", "d"))
    for (bad_args in not_enough_positional_args) {
      expect_error(parser$parse_args(bad_args), "error: the following arguments are required:")
    }
  })

  test_that("positional arguments for files are parsed correctly", {
    flags <- c("m", "d", "o")
    opts <- parser$parse_args(flags)
    expect_equal(opts$METADATA, flags[1])
    expect_equal(opts$DATA, flags[2])
    expect_equal(opts$OUTPUT, flags[3])
  })

  test_that("data_dir flag behaves correctly", {
    # default
    flags <- str_split_1("m d o", " ")
    opts <- parser$parse_args(flags)
    expect_equal(opts$data_dir, ".")

    # set another value
    flags <- str_split_1("--data_dir /data m d o", " ")
    opts <- parser$parse_args(flags)
    expect_equal(opts$data_dir, "/data")
  })
})

# for each flag, provide the optiions the flag can be as well as a value to set the flag to that isn't the default
# this is kept separate from the list in the options.R file so flag changes will be detected in tests
test_flag_values <- list(
  taxa_hfe_base_flags=list(
    subject_identifier=list(flags=list("-s", "--subject_identifier"), value="sid"),
    label=list(flags=list("-l", "--label"), value="label_factor"),
    feature_type=list(flags=list("-t", "--feature_type"), value="ftype"),
    abundance=list(flags=list("-a", "--abundance"), value=0.9),
    prevalence=list(flags=list("-p", "--prevalence"), value=0.02),
    lowest_level=list(flags=list("-L", "--lowest_level"), value=0),
    max_level=list(flags=list("-m", "--max_level"), value=20),
    cor_level=list(flags=list("-c", "--cor_level"), value=.99),
    disable_super_filter=list(flags=list("-d", "--disable_super_filter"), value=TRUE),
    write_old_files=list(flags=list("-w", "--write_old_files"), value=TRUE),
    write_flattened_tree=list(flags=list("-W", "--write_flattened_tree"), value=TRUE),
    write_both_outputs=list(flags=list("-D", "--write_both_outputs"), value=TRUE),
    nperm=list(flags=list("--nperm"), value=100),
    ncores=list(flags=list("-n", "--ncores"), value=1),
    seed=list(flags=list("--seed"), value=314159)
  ),
  taxa_hfe_ml_flags=list(
    train_split=list(flags=list("--train_split"), value=0.7),
    model=list(flags=list("--model"), value="enet"),
    folds=list(flags=list("--folds"), value=11),
    metric=list(flags=list("--metric"), value="accuracy"),
    tune_length=list(flags=list("--tune_length"), value=70),
    tune_time=list(flags=list("--tune_time"), value=1),
    tune_stop=list(flags=list("--tune_stop"), value=9),
    permute=list(flags=list("--permute"), value=2),
    shap=list(flags=list("--shap"), value=TRUE)
  )
)

test_that("loading parsers work", {
  expect_true(TRUE)

  version <- "999"
  program <- "TEST123"
  description <- "testing parser"
  base_flags <- c("m", "d", "o")

  test_that("parsers load with default values", {
    parsers <- list(load_taxa_hfe_parser(version), load_taxa_hfe_ml_parser(version))

    for (parser in parsers) {
      expect_no_error(parser$parse_args(base_flags))
    }
  })

  test_that("parsers set flags as expected", {
    # mapping of parser to the flags it handles, including a value to test setting the flag with
    parser_flag_values_map <- list(
      # base taca hfe flags
      list(
        parser=load_taxa_hfe_parser(version),
        flag_values=test_flag_values$taxa_hfe_base_flags
      ),
      # taxa hfe ml flags, includes base taxa hfe flags
      list(
        parser=load_taxa_hfe_ml_parser(version),
        flag_values=c(test_flag_values$taxa_hfe_base_flags, test_flag_values$taxa_hfe_ml_flags)
      )
    )

    for (parser_flag_values in parser_flag_values_map) {
      parser <- parser_flag_values$parser
      flags_to_test <- parser_flag_values$flag_values

      # get default opts for all flags to compare against set values
      default_opts <- parser$parse_args(base_flags)

      for (flag_name in names(flags_to_test)) {
        flag_test_obj <- flags_to_test[[flag_name]]

        # for each actual flag form (ex. -s, --seed), add to the base flags
        # also add in a as.character(value) for the flag if it accepts one
        # ensure the value is set correctly in the parsed options
        # and also does not match the default
        for (flag in flag_test_obj$flags) {
          flags <- c(base_flags, flag)
          if (!is.logical(flag_test_obj$value)) {
            flags <- c(flags, as.character(flag_test_obj$value))
          }
          opts <- parser$parse_args(flags)

          # ensure the value is set correctly in the parsed options
          # and also does not match the default
          # wrap the values in identical here because a flag without a default is NULL and can't be compared directly to its set value
          expect_equal(opts[[flag_name]], flag_test_obj$value)
          expect_false(identical(opts[[flag_name]], default_opts[[flag_name]]))
        }
      }
    }
  })
})

test_that("load_args works correctly", {
  expect_true(TRUE)

  version <- "999"
  program <- "TEST123"
  description <- "testing parser"

  test_that("defaults to version 0 when no version is set in the env", {
    commandArgs <<- function(x) {
      c("-v")
    }
    Sys.unsetenv("TAXA_HFE_VERSION")

    expect_error(load_args(load_taxa_hfe_parser), "version requested:\n0")
  })

  test_that("use env var version when set", {
    commandArgs <<- function(x) {
      c("-v")
    }
    Sys.setenv(TAXA_HFE_VERSION=version)

    expect_error(load_args(load_taxa_hfe_parser), paste("version requested:\n", version, sep=""))
  })

  test_that("it sets the paths using --data_dir", {
    commandArgs <<- function(x) {
      c("--data_dir", "/path", "m.txt", "d.txt", "o.txt")
    }

    expect_no_error(opts <- load_args(load_taxa_hfe_parser))
    expect_equal(opts$METADATA, "/path/m.txt")
    expect_equal(opts$DATA, "/path/d.txt")
    expect_equal(opts$OUTPUT, "/path/o.txt")
  })

  test_that("ignores --data_dir when paths are already linked to real files", {
    # need to create tmp files so they will be detected and no resolved against
    tmp_files <- list("/tmp/m.txt", "/tmp/d.txt", "/tmp/o.txt")
    for (f in tmp_files) {
      file.create(f)
    }

    commandArgs <<- function(x) {
      c("--data_dir", "/path", tmp_files)
    }

    expect_no_error(opts <- load_args(load_taxa_hfe_parser))
    expect_equal(opts$METADATA, "/tmp/m.txt")
    expect_equal(opts$DATA, "/tmp/d.txt")
    expect_equal(opts$OUTPUT, "/tmp/o.txt")
  })
})

test_that("load fuction don't error", {
  commandArgs <<- function(x) {
    c("m", "d", "o")
  }

  expect_no_error(load_taxa_hfe_args())
  expect_no_error(load_taxa_hfe_ml_args())
})
