#' @useDynLib fastperm C_test_get_alphabetagamma
#' @useDynLib fastperm C_test_get_alphabetagamma_input


#' @export
test_get_alphabetagamma= function(){
  
  res = .Call(C_test_get_alphabetagamma)
  
  return
}


#' @export
test_get_alphabetagamma_input= function(x,a,b,debug){
  n = length(a)
  res = .Call(C_test_get_alphabetagamma_input, as.numeric(x), as.numeric(a), as.numeric(b), 
              as.integer(n),as.integer(debug))
  
  return(res)
}