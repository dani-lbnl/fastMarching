#R has no builtin meshgrid() function but you can write one:
#Ex: call the function using yourvar=meshgrid(1:5,10:12);

meshgrid <- function(a,b) {
  list(
       x=outer(b*0,a,FUN="+"),
       y=outer(b,a*0,FUN="+")
       )
} 
