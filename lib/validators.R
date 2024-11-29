# validation functions against flags
# values for min and max are numeric, and values for min_warning and max_warning are lists with a numeric and warning message
# ex. list(1, "are you crazy???")
validate_numeric <- function (min=NULL, max=NULL, min_warning=NULL, max_warning=NULL) {
  return(function(flag_name, flag_value, all_flags) {
    if (!is.null(min)) {
      if (!is.numeric(min)) {
        stop(sprintf("'min' validate_numeric value for %s must be numeric", flag_name))
      }
      if (flag_value < min) {
        stop(sprintf("Flag %s must be greater than or equal to %s (current value: %s)", flag_name, min, flag_value)) 
      }
    }
    
    if (!is.null(max)) {
      if (!is.numeric(max)) {
        stop(sprintf("'max' validate_numeric value for %s must be numeric", flag_name))
      }
      if (flag_value > max) {
        stop(sprintf("Flag %s must be less than or equal to %s (current value: %s)", flag_name, max, flag_value)) 
      }
    }

    if (!is.null(min_warning)) {
      if (!(is.list(min_warning) && length(min_warning) == 2 && is.numeric(min_warning[[1]]))) {
        stop(sprintf("'min_warning' validate_numeric value for %s must be a length 2 list of (numeric, character)", flag_name))
      }
      if (flag_value < min_warning[[1]]) {
        warning(sprintf("Warning: for best results, flag %s should be greater than or equal to %s (current value: %s); %s", flag_name, min_warning[[1]], flag_value, min_warning[[2]])) 
      }
    }

    if (!is.null(max_warning)) {
      if (!(is.list(max_warning) && length(max_warning) == 2 && is.numeric(max_warning[[1]]))) {
        stop(sprintf("'max_warning' validate_numeric value for %s must be a length 2 list of (numeric, character)", flag_name))
      }
      if (flag_value > max_warning[[1]]) {
        warning(sprintf("Warning: for best results, flag %s should be less than or equal to %s (current value: %s); %s", flag_name, max_warning[[1]], flag_value, max_warning[[2]])) 
      }
    }
  })
}