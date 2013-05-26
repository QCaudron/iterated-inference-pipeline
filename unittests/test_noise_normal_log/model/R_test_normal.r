
require(coda)

d<-read.csv("best_0.csv",header=TRUE)

ess = round(effectiveSize(d$r0.city1__all))

step = floor(length(d$r0.city1__all)/ess)
x = d$r0.city1__all[seq(1,length(d$r0.city1__all),step)]
x2 = hist(x,12+qnorm(seq(0,1,0.1)),plot=FALSE)$counts
t= chisq.test(x2)
print(t)

if(t$p.value <0.001){
	cat("0",file="outfile.txt")
} else {
	cat("1",file="outfile.txt")
}
