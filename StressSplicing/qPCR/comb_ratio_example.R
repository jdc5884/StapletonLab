### developing code for combination ratios for qPCR ###
name = c('x','x','y','y','z','z')
a = c(1,2,3,4,5,6)
b = c(7,8,9,10,11,12)
m = data.frame(cbind(name,a,b), stringsAsFactors = FALSE)

divide <- function(col1, col2){
  r = c()
  for (i in col1){
    r = cbind(r,col1[i]/col2)
  }
  return(r) 
}

m$a = as.numeric(m$a)
m$b = as.numeric(m$b)


group = split(m, name)
for (k in group){
  #print(divide(k$a, k$b))
  print(k)
}
  
  
  
  
  
  
  
  
  
            