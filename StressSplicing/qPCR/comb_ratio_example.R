### developing code for combination ratios for qPCR ###
name = c('x','x','y','y','z','z')
a = c(1,2,3,4,5,6)
b = c(7,8,9,10,11,12)
m = data.frame(cbind(name,a,b), stringsAsFactors = FALSE)

# this function divides every item in column 1 by every item in column 2
divide <- function(col1, col2){
  r = c()
  for (i in col1){
    r = cbind(r,col1[i]/col2)
  }
  return(r) 
}

# change the components of the matrix to numeric
m$a = as.numeric(m$a)
m$b = as.numeric(m$b)

# split the matrix by each group
group = split(m, name)
print(group)

#for each group in k, divide column a by column b

for (k in group){
  print(divide(k$a, k$b))
  #print(k)
}
 

  
  
  
            