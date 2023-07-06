#' @useDynLib fastperm C_test_sparse_functions_with_example


#' @export
test_rev_sparse_functions= function(){
  
  res = .Call(C_test_sparse_functions_with_example)
  
  return
}
