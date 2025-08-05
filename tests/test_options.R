library(testthat)
library(stringr)

source("lib/options.R")

# grab binding to these functions so they can be overwritten in the tests
commandArgs <- NULL
quit <- NULL

test_that("initialize_parser works as expected", {
  expect_true(TRUE)

  version <- "998"
  program <- "TEST123"
  description <- "testing parser"

  parser <- initialize_parser(version, program, description, list())

  test_that("parser prints version number when passed -v", {
    flags <- c("-v")
    # mock this method from the argparse package so it doesn't call quit but still raises the error it would have with stop()
    local_mocked_bindings(
      print_message_and_exit = function(message, x) stop(message),
      .package = 'argparse'
    )

    expect_error(parser$parse_args(flags), version)
  })

  test_that("parser prints help message when passed -h", {
    flags <- c("-h")
    # mock this method from the argparse package so it doesn't call quit but still raises the error it would have with stop()
    local_mocked_bindings(
      print_message_and_exit = function(message, x) stop(message),
      .package = 'argparse'
    )

    expect_error(parser$parse_args(flags), "METADATA DATA")
    expect_error(parser$parse_args(flags), program)
    expect_error(parser$parse_args(flags), description)
  })

  test_that("parser requires 2 positional arguments for the files", {
    not_enough_positional_args <- list(c(), c("m"), c("d"))
    # mock this method from the argparse package so it doesn't call quit but still raises the error it would have with stop()
    local_mocked_bindings(
      print_message_and_exit = function(message, x) stop(message),
      .package = 'argparse'
    )

    for (bad_args in not_enough_positional_args) {
      expect_error(parser$parse_args(bad_args), "error: the following arguments are required:")
    }
  })

  test_that("positional arguments for files are parsed correctly", {
    flags <- c("m", "d")
    opts <- parser$parse_args(flags)
    expect_equal(opts$METADATA, flags[1])
    expect_equal(opts$DATA, flags[2])
  })

  test_that("data_dir flag behaves correctly", {
    # default
    flags <- str_split_1("m d", " ")
    opts <- parser$parse_args(flags)
    expect_equal(opts$data_dir, ".")

    # set another value
    flags <- str_split_1("--data_dir /data m d", " ")
    opts <- parser$parse_args(flags)
    expect_equal(opts$data_dir, "/data")
  })

  test_that("argument group loading works", {
    arg_group <- list(
      name="arg group name",
      desc="arg group description",
      args=list(
        foo=list("-f", "--foo", type="character", metavar="<string>", default="bar", help="test arg in arg group")
      )
    )

    parser <- initialize_parser(version, program, description, list(arg_group))

    test_that("help text includes arg group info", {
      # mock this method from the argparse package so it doesn't call quit but still raises the error it would have with stop()
      local_mocked_bindings(
        print_message_and_exit = function(message, x) stop(message),
        .package = 'argparse'
      )
      flags <- c("-h")

      expect_error(parser$parse_args(flags), arg_group$name)
      expect_error(parser$parse_args(flags), arg_group$desc)
      expect_error(parser$parse_args(flags), arg_group$args$foo$help)
    })

    flags <- c("m", "d", "-f", "baz")
    opts <- parser$parse_args(flags)
    expect_equal(opts$foo, "baz")
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
    # mock this method from the argparse package so it doesn't call quit but still raises the error it would have with stop()
    local_mocked_bindings(
      print_message_and_exit = function(message, x) stop(message),
      .package = 'argparse'
    )

    Sys.unsetenv("TAXA_HFE_VERSION")

    expect_error(load_args(program, description, argument_groups = list()), "0")
  })

  test_that("use env var version when set", {
    commandArgs <<- function(x) {
      c("-v")
    }
    # mock this method from the argparse package so it doesn't call quit but still raises the error it would have with stop()
    local_mocked_bindings(
      print_message_and_exit = function(message, x) stop(message),
      .package = 'argparse'
    )

    Sys.setenv(TAXA_HFE_VERSION=version)

    expect_error(load_args(program, description, argument_groups = list()), version)
  })

  test_that("it sets the paths using --data_dir, trims slash and creates output dir", {
    commandArgs <<- function(x) {
      c("--data_dir", "/path", "m.txt", "d.txt", "-o", "custom_outputs/")
    }

    expect_no_error(opts <- load_args(program, description, argument_groups = list()))
    expect_equal(opts$METADATA, "/path/m.txt")
    expect_equal(opts$DATA, "/path/d.txt")
    # also tests that the output directory gets the trailing slash trimmed
    expect_equal(opts$output_dir, "/path/custom_outputs")
    dir.exists("/path/custom_outputs")
  })

  test_that("ignores --data_dir when paths are already linked to real files", {
    # need to create tmp files so they will be detected and no resolved against
    tmp_files <- list("/tmp/m.txt", "/tmp/d.txt", "-o", "/tmp/custom_outputs")
    file.create("/tmp/m.txt")
    file.create("/tmp/d.txt")

    commandArgs <<- function(x) {
      c("--data_dir", "/path", tmp_files)
    }

    expect_no_error(opts <- load_args(program, description, argument_groups = list()))
    expect_equal(opts$METADATA, "/tmp/m.txt")
    expect_equal(opts$DATA, "/tmp/d.txt")
    expect_equal(opts$output_dir, "/tmp/custom_outputs")
  })

  test_that("load_args() always sets the global seed", {
    expect_true(TRUE)

    test_that("sets the seed when no seed is provided", {
      # no seed set
      # will use default_seed()
      commandArgs <<- function(x) {
        c("m", "d")
      }
      
      # cache the random seed vector and ensure that it is different after running load args
      rng_state <- .Random.seed
      expect_no_error(opts <- load_args(program, description, argument_groups = list()))
      expect_equal(all(rng_state == .Random.seed), FALSE)
    })

    test_that("sets the seed from the flag", {
      # seed set
      seed <- 42
      # will be used directly in set.seed()
      commandArgs <<- function(x) {
        c("--seed", as.character(seed), "m", "d")
      }

      # cache the random seed vector and ensure that it is different after running load args
      rng_state <- .Random.seed
      expect_no_error(opts <- load_args(program, description, argument_groups = list()))
      expect_equal(opts$seed, seed)
      expect_equal(all(rng_state == .Random.seed), FALSE)
    })
  })
})

# for each flag, provide the options the flag can be as well as a value to set the flag to that isn't the default
# this is kept separate from the list in the options.R file so flag changes will be detected in tests
# additionally, any errors=list(), warnings=list() values will be checked for that output
test_flag_values <- list(
  taxa_hfe_base_args=list(
    subject_identifier=list(flags=list("-s", "--subject_identifier"), value="sid"),
    label=list(flags=list("-l", "--label"), value="label_factor"),
    feature_type=list(flags=list("-t", "--feature_type"), value="ftype"),
    random_effects=list(flags=list("-R", "--random_effects"), value=TRUE),
    k_splits=list(flags=list("-k", "--k_splits"), value=4, errors=list(-1), warnings=list(10)),
    abundance=list(flags=list("-a", "--abundance"), value=0.9, errors=list(-1)),
    prevalence=list(flags=list("-p", "--prevalence"), value=0.02, errors=list(-1,2)),
    lowest_level=list(flags=list("-L", "--lowest_level"), value=4, warnings=list(1), errors=list(-1)),
    max_level=list(flags=list("-m", "--max_level"), value=10, warnings=list(30, 500)),
    cor_level=list(flags=list("-c", "--cor_level"), value=.99, errors=list(-1, 2), warnings=list(.4)),
    disable_super_filter=list(flags=list("-d", "--disable_super_filter"), value=TRUE),
    write_old_files=list(flags=list("-w", "--write_old_files"), value=TRUE),
    write_flattened_tree=list(flags=list("-W", "--write_flattened_tree"), value=TRUE),
    write_both_outputs=list(flags=list("-D", "--write_both_outputs"), value=TRUE),
    nperm=list(flags=list("--nperm"), value=100),
    ncores=list(flags=list("-n", "--ncores"), value=1, errors=list(-4,0)),
    # this is technically part of the initialized parser args but is tested as a part of this group
    seed=list(flags=list("--seed"), value=314159, errors=list("-1000000000000000", "1000000000000000"))
  ),
  taxa_hfe_ml_args=list(
    train_split=list(flags=list("--train_split"), value=0.7, errors=list(-1,2), warnings=list(0.4)),
    model=list(flags=list("--model"), value="enet"),
    folds=list(flags=list("--folds"), value=5, errors=list(-1,0), warnings=list(12)),
    cv_repeats=list(flags=list("--cv_repeats"), value=2, errors=list(-1,0), warnings=list(6)),
    metric=list(flags=list("--metric"), value="accuracy"),
    tune_length=list(flags=list("--tune_length"), value=70),
    tune_time=list(flags=list("--tune_time"), value=1, errors=list(-1,0), warnings=list(30)),
    tune_stop=list(flags=list("--tune_stop"), value=9),
    permute=list(flags=list("--permute"), value=2, errors=list(-1,0), warnings=list(50)),
    shap=list(flags=list("--shap"), value=TRUE),
    summarized_levels=list(flags=list("--summarized_levels"), value=TRUE)
  )
)

test_that("program arg loaders work", {
  expect_true(TRUE)

  version <- "999"
  program <- "TEST123"
  description <- "testing parser"
  commandArgs <<- function(x) {
    c("m", "d")
  }

  # mapping of parser to the flags it handles, including a value to test setting the flag with
  # when adding new argument groups, add a new item to this list
  parser_flag_values_map <- list(
    # base taxa hfe flags
    taxa_hfe_base_args=list(
      load_arg_function=load_taxa_hfe_args,
      flag_values=test_flag_values$taxa_hfe_base_args
    ),
    # taxa hfe ml flags, includes base taxa hfe flags
    taxa_hfe_ml_args=list(
      load_arg_function=load_taxa_hfe_ml_args,
      flag_values=c(test_flag_values$taxa_hfe_base_args, test_flag_values$taxa_hfe_ml_args)
    )
  )

  test_that("every flag is tested", {
    for (group_name in names(test_flag_values)) {
      expect_true(group_name %in% names(argument_groups), info = sprintf("%s flag group not tested in test_flag_values", group_name))

      test_flags <- test_flag_values[[group_name]]
      actual_flags <- argument_groups[[group_name]]
      for (flag in names(actual_flags$args)) {
        expect_true(flag %in% names(test_flags), info = sprintf("%s is not tested in test_flag_values, add flag to correct test argument group", flag))
      }
    }
  })

  test_that("parsers load with default values", {
    for (parser in parser_flag_values_map) {
      expect_no_error(parser$load_arg_function())
    }
  })

  test_that("parsers set flags as expected", {
    expect_true(TRUE)

    base_flags <- c("m.txt", "d.txt")

    for (parser_name in names(parser_flag_values_map)) {
      parser_flag_values <- parser_flag_values_map[[parser_name]]
      load_arg_function <- parser_flag_values$load_arg_function
      flags_to_test <- parser_flag_values$flag_values

      # get default opts for all flags to compare against set values
      commandArgs <<- function(x) {
        base_flags
      }
      default_opts <- load_arg_function()

      for (flag_name in names(flags_to_test)) {
        flag_test_obj <- flags_to_test[[flag_name]]

        # for each actual flag form (ex. -s, --seed), add to the base flags
        # also add in a as.character(value) for the flag if it accepts one
        # ensure the value is set correctly in the parsed options
        # and also does not match the default
        for (flag in flag_test_obj$flags) {
          test_that("valid non-default value is parsed for flag", {
            # test a good non-default value works
            flags <- c(base_flags, flag)
            if (!is.logical(flag_test_obj$value)) {
              flags <- c(flags, as.character(flag_test_obj$value))
            }
            commandArgs <<- function(x) {
              flags
            }
            quit_called <<- FALSE
            quit <<- function(...) {
              quit_called <<- TRUE
            }

            expect_no_warning(opts <- load_arg_function(), message = sprintf("parser: %s, flag_name: %s, flag: %s, value: %s - parsing valid value produces warning", parser_name, flag_name, flag, flag_test_obj$value))

            # ensure the value is set correctly in the parsed options
            # and also does not match the default
            # wrap the values in identical here because a flag without a default is NULL and can't be compared directly to its set value
            expect_equal(opts[[flag_name]], flag_test_obj$value, info = sprintf("parser: %s, flag_name: %s, flag: %s, value: %s - flag not set to expected value", parser_name, flag_name, flag, flag_test_obj$value))
            expect_false(identical(opts[[flag_name]], default_opts[[flag_name]]), info = sprintf("parser: %s, flag_name: %s, flag: %s, value: %s - flag set to default value", parser_name, flag_name, flag, flag_test_obj$value))
            expect_false(quit_called, info = sprintf("parser: %s, flag_name: %s, flag: %s, value: %s - quit called during valid value test", parser_name, flag_name, flag, flag_test_obj$value))
          })

          test_that("flag parse generates error/quit when expected", {
            expect_true(TRUE)

            # test that errors produce errors
            for (error_value in flag_test_obj$errors) {
              # add flag to base values
              if (is.logical(error_value)) {
                if (error_value) {
                  flags <- c(base_flags, flag)
                }
              } else {
                flags <- c(base_flags, flag)
                flags <- c(flags, as.character(error_value))
              }

              # overload the quit() function
              # instead of exiting, it will call stop() so the overall load_args() function still doesn't proceed but we can also catch the error
              # this is why we expect_error around expect_message
              quit <<- function(...) {
                stop("quit")
              }
              commandArgs <<- function(x) {
                flags
              }

              # expect message for the actual flag error, but also expect an error from the new stop() call in the mocked quit function above
              expect_error(
                expect_message(load_arg_function(), regexp = "Flag .+? must be", info = sprintf("parser: %s, flag_name: %s, flag: %s, value: %s - expected flag error message", parser_name, flag_name, flag, error_value)),
                regexp = "quit",
                info = sprintf("parser: %s, flag_name: %s, flag: %s, value: %s - expected flag error did not call quit() on error", parser_name, flag_name, flag, error_value)
              )
            }
          })

          test_that("flag parse generates warnings when expected", {
            expect_true(TRUE)

            # test that warning values are produced
            for (warning_value in flag_test_obj$warnings) {
              # add flag to base values
              if (is.logical(warning_value)) {
                if (warning_value) {
                  flags <- c(base_flags, flag)
                }
              } else {
                flags <- c(base_flags, flag)
                flags <- c(flags, as.character(warning_value))
              }

              quit_called <<- FALSE
              quit <<- function(...) {
                quit_called <<- TRUE
              }
              commandArgs <<- function(x) {
                flags
              }

              expect_warning(load_arg_function(), info = sprintf("parser: %s, flag_name: %s, flag: %s, value: %s - expected flag warning did not produce warning message", parser_name, flag_name, flag, warning_value))
              expect_equal(quit_called, FALSE, info = sprintf("parser: %s, flag_name: %s, flag: %s, value: %s - expected warning produced a quit instead", parser_name, flag_name, flag, warning_value))
            }
          })
        }
      }
    }
  })
})
