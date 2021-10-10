capture_error <- function(code, otherwise = NULL, quiet = TRUE) {
  tryCatch(list(result = code, error = NULL), error = function(e) {
    if (!quiet) {
      message("Error: ", e$message)
    }
    list(result = otherwise, error = e)
  }, interrupt = function(e) {
    stop("Terminated by user", call. = FALSE)
  })
}

safely <- function(.f, otherwise = NULL, quiet = TRUE) {
  function(...) capture_error(.f(...), otherwise, quiet)
}

possibly <- function(.f, otherwise, quiet = TRUE) {
  force(otherwise)
  function(...) capture_error(.f(...), otherwise, quiet)$result
}

quietly <- function(.f) {
  function(...) suppressMessages(suppressWarnings(.f(...)))
}
