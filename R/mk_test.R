#' test for directional gradient
#'
#' @param df data frame representing least cost path
#' @param yvar column name for y variable
#' @param mk_alpha significance threshold for MK test
#' @param mk_tau_threshold threshold for tau
#'
#' @returns results of MK test as list
#' @export
#'
#' @examples
#' # test_mannkendall(lc_path, yvar = 'TSM')
test_mannkendall <- function(df, yvar, mk_alpha = 0.05, mk_tau_threshold = 0.5){
  mk_test <- Kendall::MannKendall(df[[yvar]])
  mk_true <- mk_test$sl < mk_alpha
  mk_pos <- mk_test$tau > 0
  mk_gt50 <- mk_test$tau >= mk_tau_threshold
  mk_truepos <- mk_true & mk_pos & mk_gt50
  return(list(mk_truepos = mk_truepos,
              mk_tau = mk_test$tau))
}
