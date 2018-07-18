

modfit<-function(times,gexp,plot=T){
	#estimate linear regression model - y=a+bx 
	lin<-lm(gexp~times)
	#estimate non-linear regression model y = a+b^x
	nonlin<-nls(gexp~a*exp(b*times),start=list(a=0.1,b=0.01),
		control=nls.control(maxiter=1e3,minFactor=1e-20))
	
	#compute rho for each model fit
	lincor<-cor(gexp,predict(lin))
	nonlincor<-cor(gexp,predict(nonlin))
		
	#return slope, rho, p-value
	col<-c('slope','rho','p')
	lintab<-c(summary(lin)$coefficients[2,1],lincor,summary(lin)$coefficients[2,4])
	names(lintab)<-col
	nonlintab<-c(summary(nonlin)$coefficients[2,1],nonlincor,summary(nonlin)$coefficients[2,4])
	names(nonlintab)<-col
	
		res<-list(linear_fit=lintab,nonlinear_fit=nonlintab,lin,nonlin)
	#plot best-fit linear and non-linear slopes to data
	if(plot==T){
			df<-data.frame(times=times,gexp=gexp)
			xci<-seq(min(df$times),max(df$times),length.out=50)
			par(mfrow=c(1,2),mar=c(4,4,1,1))
			plot(times,gexp,xlab="time",ylab="gexp")
				lines(predict(lin),col=colors(1)[runif(1,10,502)],lwd=2)
			plot(times,gexp,xlab="time",ylab="gexp")
				lines(predict(nonlin),col=colors(1)[runif(1,10,502)],lwd=2)
		}
		return(res)
}	
