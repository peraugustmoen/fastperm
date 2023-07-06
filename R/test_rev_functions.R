#' @useDynLib fastperm C_test_functions_with_example


#' @export
test_rev_functions= function(){
  
  res = .Call(C_test_functions_with_example)
  
  return
}
