library(testthat)
library(stringr)

source("lib/options.R")

# grab binding to these functions so they can be overwritten in the tests
commandArgs <- NULL

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

    expect_error(parser$parse_args(flags), "METADATA DATA OUTPUT")
    expect_error(parser$parse_args(flags), program)
    expect_error(parser$parse_args(flags), description)
  })

  test_that("parser requires 3 positional arguments for the files", {
    not_enough_positional_args <- list(c(), c("m"), c("m", "d"))
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

    flags <- c("m", "d", "o", "-f", "baz")
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

  test_that("it sets the paths using --data_dir", {
    commandArgs <<- function(x) {
      c("--data_dir", "/path", "m.txt", "d.txt", "o.txt")
    }

    expect_no_error(opts <- load_args(program, description, argument_groups = list()))
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

    expect_no_error(opts <- load_args(program, description, argument_groups = list()))
    expect_equal(opts$METADATA, "/tmp/m.txt")
    expect_equal(opts$DATA, "/tmp/d.txt")
    expect_equal(opts$OUTPUT, "/tmp/o.txt")
  })
})

# for each flag, provide the optiions the flag can be as well as a value to set the flag to that isn't the default
# this is kept separate from the list in the options.R file so flag changes will be detected in tests
test_flag_values <- list(
  taxa_hfe_base_args=list(
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
  taxa_hfe_ml_args=list(
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

test_that("program arg loaders work", {
  expect_true(TRUE)

  version <- "999"
  program <- "TEST123"
  description <- "testing parser"
  commandArgs <<- function(x) {
    c("m", "d", "o")
  }

  # mapping of parser to the flags it handles, including a value to test setting the flag with
  # when adding new argument groups, add a new item to this list
  parser_flag_values_map <- list(
    # base taxa hfe flags
    list(
      load_arg_function=load_taxa_hfe_args,
      flag_values=test_flag_values$taxa_hfe_base_args
    ),
    # taxa hfe ml flags, includes base taxa hfe flags
    list(
      load_arg_function=load_taxa_hfe_ml_args,
      flag_values=c(test_flag_values$taxa_hfe_base_args, test_flag_values$taxa_hfe_ml_args)
    )
  )

  test_that("every flag is tested", {
    for (group_name in names(test_flag_values)) {
      test_flags <- test_flag_values[[group_name]]
      actual_flags <- argument_groups[[group_name]]
      for (flag in names(actual_flags$args)) {
        expect_true(flag %in% names(test_flags), info = flag)
      }
    }
  })

  test_that("parsers load with default values", {
    for (parser in parser_flag_values_map) {
      expect_no_error(parser$load_arg_function())
    }
  })

  test_that("parsers set flags as expected", {
    base_flags <- c("m.txt", "d.txt", "o.txt")

    for (parser_flag_values in parser_flag_values_map) {
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
          flags <- c(base_flags, flag)
          if (!is.logical(flag_test_obj$value)) {
            flags <- c(flags, as.character(flag_test_obj$value))
          }
          commandArgs <<- function(x) {
            flags
          }
          opts <- load_arg_function()

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