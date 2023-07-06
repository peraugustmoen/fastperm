#' @useDynLib fastperm C_test_get_alphabetagamma


#' @export
test_get_alphabetagamma= function(){
  
  res = .Call(C_test_get_alphabetagamma)
  
  return
}
