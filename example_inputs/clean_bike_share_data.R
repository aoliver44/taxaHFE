## massage bikeshare data

## download data from https://doi.org/10.24432/C5W894

## read in hours data
hours <- readr::read_csv(file = "/home/rstudio/example_inputs/bike_share_hour.csv")
days <- readr::read_csv(file = "/home/rstudio/example_inputs/bike_share_day.csv")

## modify columns
hours_modified <- hours %>%
  ## create subject_identifier name, convert existing instant column
  dplyr::mutate(., instant = paste0("record_", instant)) %>%
  ## remove date
  dplyr::select(., -dteday) %>%
  ## convert season to prettier factor
  dplyr::mutate(., season = dplyr::case_when(season == 1 ~ "spring",
                                             season == 2 ~ "summer",
                                             season == 3 ~ "fall", 
                                             season == 4 ~ "winter", 
                                             .default = "no_season_info")) %>%
  ## make year prettier
  dplyr::mutate(., yr = dplyr::case_when(yr == 0 ~ "2011yr",
                                             yr == 1 ~ "2012yr",
                                             .default = "no_year_info")) %>%
  ## make month prettier
  dplyr::mutate(., mnth = dplyr::case_when(mnth == 1 ~ "jan",
                                           mnth == 2 ~ "feb",
                                           mnth == 3 ~ "mar",
                                           mnth == 4 ~ "apr",
                                           mnth == 5 ~ "may",
                                           mnth == 6 ~ "jun",
                                           mnth == 7 ~ "jul",
                                           mnth == 8 ~ "aug",
                                           mnth == 9 ~ "sep",
                                           mnth == 10 ~ "oct",
                                           mnth == 11 ~ "nov",
                                           mnth == 12 ~ "dec", 
                                           .default = "no_month_info")) %>%
  ## make holiday prettier
  dplyr::mutate(., holiday = dplyr::case_when(holiday == 0 ~ "no_holiday",
                                              holiday == 1 ~ "holiday",
                                         .default = "no_holiday_info")) %>%
  ## make hours prettier
  dplyr::mutate(., hr = paste0("hour_", hr)) %>%
  ## make weekday prettier
  dplyr::mutate(., weekday = paste0("weekday_", weekday)) %>%
  ## make workday prettier
  dplyr::mutate(., workingday = dplyr::case_when(workingday == 0 ~ "not_workday",
                                                 workingday == 1 ~ "workday",
                                              .default = "no_workday_info")) %>%
  ## make weathersit prettier
  dplyr::mutate(., weathersit = dplyr::case_when(weathersit == 1 ~ "Clear",
                                                 weathersit == 2 ~ "Mist",
                                                 weathersit == 3 ~ "Light Snow-Light Rain",
                                                 weathersit == 4 ~ "Heavy Rain",
                                                 .default = "no_weather_info"))
  
  
## modify columns
days_modified <- days %>%
  ## create subject_identifier name, convert existing instant column
  dplyr::mutate(., instant = paste0("record_", instant)) %>%
  ## remove date
  dplyr::select(., -dteday) %>%
  ## convert season to prettier factor
  dplyr::mutate(., season = dplyr::case_when(season == 1 ~ "spring",
                                             season == 2 ~ "summer",
                                             season == 3 ~ "fall", 
                                             season == 4 ~ "winter", 
                                             .default = "no_season_info")) %>%
  ## make year prettier
  dplyr::mutate(., yr = dplyr::case_when(yr == 0 ~ "2011yr",
                                         yr == 1 ~ "2012yr",
                                         .default = "no_year_info")) %>%
  ## make month prettier
  dplyr::mutate(., mnth = dplyr::case_when(mnth == 1 ~ "jan",
                                           mnth == 2 ~ "feb",
                                           mnth == 3 ~ "mar",
                                           mnth == 4 ~ "apr",
                                           mnth == 5 ~ "may",
                                           mnth == 6 ~ "jun",
                                           mnth == 7 ~ "jul",
                                           mnth == 8 ~ "aug",
                                           mnth == 9 ~ "sep",
                                           mnth == 10 ~ "oct",
                                           mnth == 11 ~ "nov",
                                           mnth == 12 ~ "dec", 
                                           .default = "no_month_info")) %>%
  ## make holiday prettier
  dplyr::mutate(., holiday = dplyr::case_when(holiday == 0 ~ "no_holiday",
                                              holiday == 1 ~ "holiday",
                                              .default = "no_holiday_info")) %>%
  ## make weekday prettier
  dplyr::mutate(., weekday = paste0("weekday_", weekday)) %>%
  ## make workday prettier
  dplyr::mutate(., workingday = dplyr::case_when(workingday == 0 ~ "not_workday",
                                                 workingday == 1 ~ "workday",
                                                 .default = "no_workday_info")) %>%
  ## make weathersit prettier
  dplyr::mutate(., weathersit = dplyr::case_when(weathersit == 1 ~ "Clear",
                                                 weathersit == 2 ~ "Mist",
                                                 weathersit == 3 ~ "Light Snow-Light Rain",
                                                 weathersit == 4 ~ "Heavy Rain",
                                                 .default = "no_weather_info"))

write.csv(x = hours_modified, file = "/home/rstudio/example_inputs/bike_share_hour.csv", quote = F, row.names = F)
write.csv(x = days_modified, file = "/home/rstudio/example_inputs/bike_share_day.csv", quote = F, row.names = F)

