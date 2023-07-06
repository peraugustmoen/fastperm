#' @useDynLib fastperm C_get_log_permanents

#' @export
get_log_permanents= function(X, a, b,debug){
  
  # X is n times T

  
  if(!is.matrix(X)){
    if(debug){
      print("input is vector:")
      print(X)
    }
    n = length(X)
    Tt = 1
    if(debug){
      print(n)
      print(Tt)
    }
    
    res = .Call(C_get_log_permanents, (X[]), as.numeric(a[]), as.numeric(b[]), as.integer(n), as.integer(Tt), as.integer(debug))
    #res = 1
    
  }else{
    print("input is matrix:")
    print(X)
    n = dim(X)[2]
    Tt = dim(X)[1]
    if(debug){
      print(n)
      print(Tt)
    }
    #res = 1
    res = .Call(C_get_log_permanents, t(X[,]), as.numeric(a[]), as.numeric(b[]), as.integer(n), as.integer(Tt), as.integer(debug))
  }
  
 
  
  #return(0)
  return(res)
}