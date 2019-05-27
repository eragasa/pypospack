
J=function(h){
+ fhat=Vectorize(function(x) density(X,from=x,to=x,n=1,bw=h)$y)
+ fhati=Vectorize(function(i) density(X[-i],from=X[i],to=X[i],n=1,bw=h)$y)
+ F=fhati(1:length(X))
+ return(integrate(function(x) fhat(x)^2,-Inf,Inf)$value-2*mean(F))
+ }
vx=seq(.1,1,by=.01)
vy=Vectorize(J)(vx)
df=data.frame(vx,vy)
library(ggplot2)
qplot(vx,vy,geom="line",data=df)
