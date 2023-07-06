#' @useDynLib fastperm C_wrapper_log_sum_exp


#' @export
TEST_logsumexp= function(x){
  res = .Call(C_wrapper_log_sum_exp, x, as.integer(length(x)))
  return(res)
}
