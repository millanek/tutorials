# ----------- Probability activity - code for section 4, "estimating Fst" ----------------
# Author: Milan Malinsky
# Last edit: 5th June 2022

# ----------- FUNCTIONS ----------------------

# These function implement the Hudson estimator as defined in equation 10 of Bhatia et al. (2013) Genome Research
calculateFstNumerator <- function(p_1, p_2, n1, n2) {
	power <- (p_1-p_2)^2
	fraction1 <- (p_1*(1-p_1))/(n1-1)
	fraction2 <- (p_2*(1-p_2))/(n2-1);
	numerator <- power - fraction1 - fraction2;
	return(numerator)
}
calculateFstDenominator <- function(p_1, p_2) {
	denominator <- (p_1*(1-p_2))+(p_2*(1-p_1));
	return(denominator)
}

# Calculating F_ST with true population allele frequency
calculateFstNumeratorNoNs <- function(p_1, p_2) {
	power <- (p_1-p_2)^2
	numerator <- power;
	return(numerator)
}


# ----------- ANALYSES AND PLOTS ----------------------

# 1) Sample allele frequencies (AFs) when population minor AF is 0.45
#pdf("~/WPSGsampleAF_0_45.pdf",width=9,height=5.5)
plot((0:8)/8,dbinom(0:8, size = 8, prob=0.45),pch=16,col="red",xlab="Sample allele frequency",ylab="Probability",cex=1.2,type='b',main="Population allele frequency = 0.45")
points((0:20)/20,dbinom(0:20, size = 20, prob=0.45),pch=16,col="blue",cex=1.2,type='b')
points((0:40)/40,dbinom(0:40, size = 40, prob=0.45),pch=16,col="orange",type='b',cex=1.2)
legend("topright",legend=c("4 diploid individuals","10 diploid individuals", "20 diploid individuals"),pch=16,col=c("red","blue","orange"))
#dev.off()

sum(dbinom(c(0,1,2,6,7,8), size = 8, prob=0.45))
sum(dbinom(c(1:5,15:20), size = 20, prob=0.45))
sum(dbinom(c(1:10,30:40), size = 40, prob=0.45))

# 2) Sample allele frequencies (AFs) when population minor AF is 0.15
#pdf("~/WPSGsampleAF_0_15.pdf",width=9,height=5.5)
plot((0:8)/8,dbinom(0:8, size = 8, prob=0.15),pch=16,col="red",xlab="Sample allele frequency",ylab="Probability",cex=1.2,type='b',main="Population allele frequency = 0.15")
points((0:20)/20,dbinom(0:20, size = 20, prob=0.15),pch=16,col="blue",cex=1.2,type='b')
points((0:40)/40,dbinom(0:40, size = 40, prob=0.15),pch=16,col="orange",type='b',cex=1.2)
legend("topright",legend=c("4 diploid individuals","10 diploid individuals", "20 diploid individuals"),pch=16,col=c("red","blue","orange"))
#dev.off()



# 3) Get pairwise Fst estimates for a single SNP
p1 <- 0.45; p2 <- 0.15        # Population allele frequencies in populations 1 and 2

# To get Fst based on population allele frequencies,
# we can either assume a huge sample size (e.g. 1 million) and feed it into the Hudson estimator:
calculateFstNumerator(p1,p2,1000000,1000000)/calculateFstDenominator(p1,p2)
# Alternatively, even better, we can ignore the terms adjusting for sample sizes (n1,n2) altogether:
populationFst <- calculateFstNumeratorNoNs(p1,p2)/calculateFstDenominator(p1,p2)

# Simulate the single-SNP Fst distrubutions based on sampling 4, 10, and 20 individuals from each population 
n <- c(8,20,40)
#pdf("~/WPSGsampleFST_oneSNP.pdf")
par(mfrow=c(3,1),mar=c(4.1,4.1,3.1,2.1))
Fst_samples <- list(3)
for (i in 1:3) {
	p_1_samples <- replicate(100000, sum(rbinom(n[i], size = 1, prob=p1))/n[i])
	p_2_samples <- replicate(100000, sum(rbinom(n[i], size = 1, prob=p2))/n[i])
	Fst_samples[[i]] <- calculateFstNumerator(p_1_samples, p_2_samples,n[i],n[i])/calculateFstDenominator(p_1_samples, p_2_samples)
	Fst_samples[[i]][which(Fst_samples[[i]] < 0)] <- 0
	h <- hist(Fst_samples[[i]],breaks=seq(0,1,length=40),xlab="Fst estimates",freq=T,main=paste0("Sample Fst distribution; ",n[i]/2," diploid individuals"))
	abline(v=populationFst,col="red",lwd=2)
	text(0.8,max(h$counts)*0.9,paste0("population Fst: ", round(populationFst,3)))
	text(0.8,max(h$counts)*0.8,paste0("mean sample Fst: ", round(mean(Fst_samples[[i]],na.rm=T),3)))
	text(0.8,max(h$counts)*0.7,paste0("sd of sample Fst: ", round(sd(Fst_samples[[i]],na.rm=T),3)))
	#which(is.na(calculateFstNumerator(p_1_samples, p_2_samples,n[i],n[i])/calculateFstDenominator(p_1_samples, p_2_samples)))
}
#dev.off()


# 4) Simulate Fst distributions when averaging over 5 SNPs: 
# Population allele frequencies:
p1s_5 <- c(0.375,0.018,0.159,0.273,0.064); p2s_5 <- c(0.393,0.308,0.188,0.278,0.49)
populationFstNumerator_5 <- sum(calculateFstNumeratorNoNs(p1s_5,p2s_5))
populationFstDenominator_5 <- sum(calculateFstDenominator(p1s_5,p2s_5))
populationFst_5 <- populationFstNumerator_5/populationFstDenominator_5

#pdf("~/WPSGsampleFST_fiveSNPs.pdf")
par(mfrow=c(3,1),mar=c(4.1,4.1,3.1,2.1))
Fst_samples_5 <- list(3)
for (i in 1:3) {
	FstSampleNums_5 <- numeric(0); FstSampleDenoms_5 <- numeric(0); 
	for (SNPnum in 1:5) {
		p_1_samples <- replicate(100000, sum(rbinom(n[i], size = 1, prob=p1s_5[SNPnum]))/n[i])
		p_2_samples <- replicate(100000, sum(rbinom(n[i], size = 1, prob=p2s_5[SNPnum]))/n[i])
		FstSampleNums_5 <- cbind(FstSampleNums_5,calculateFstNumerator(p_1_samples, p_2_samples,n[i],n[i]))
		FstSampleDenoms_5 <- cbind(FstSampleDenoms_5,calculateFstDenominator(p_1_samples, p_2_samples))
	}
	Fst_samples_5[[i]] <-apply(FstSampleNums_5,1,sum)/apply(FstSampleDenoms_5, 1,sum)
	Fst_samples_5[[i]][which(Fst_samples_5[[i]] < 0)] <- 0
	h <- hist(Fst_samples_5[[i]],breaks=seq(0,1,length=40),xlab="Fst estimates from 5 independent SNPs",freq=T,main=paste0("Sample Fst distribution; 5 SNPs; ",n[i]/2," diploid individuals"))
	abline(v= populationFst_5,col="red",lwd=2)
	text(0.8,max(h$counts)*0.9,paste0("population Fst: ", round(populationFst_5,3)))
	text(0.8,max(h$counts)*0.8,paste0("mean sample Fst: ", round(mean(Fst_samples_5[[i]],na.rm=T),3)))
	text(0.8,max(h$counts)*0.7,paste0("sd of sample Fst: ", round(sd(Fst_samples_5[[i]],na.rm=T),3)))
}
#dev.off()


# 5) Simulate Fst distributions when averaging over 10 SNPs: 
# Population allele frequencies:
p1s_10 <- c(0.421, 0.154, 0.055, 0.197, 0.062, 0.043, 0.094, 0.343, 0.314, 0.024); p2s_10 <- c(0.03, 0.423, 0.052, 0.407, 0.157, 0.038, 0.285, 0.4, 0.172, 0.352)
populationFstNumerator_10 <- sum(calculateFstNumeratorNoNs(p1s_10,p2s_10))
populationFstDenominator_10 <- sum(calculateFstDenominator(p1s_10,p2s_10))
populationFst_10 <- populationFstNumerator_10/populationFstDenominator_10

#pdf("~/WPSGsampleFST_tenSNPs.pdf")
par(mfrow=c(3,1),mar=c(4.1,4.1,3.1,2.1))
Fst_samples_10 <- list(3)
for (i in 1:3) {
	FstSampleNums_10 <- numeric(0); FstSampleDenoms_10 <- numeric(0); 
	for (SNPnum in 1:10) {
		p_1_samples <- replicate(100000, sum(rbinom(n[i], size = 1, prob=p1s_10[SNPnum]))/n[i])
		p_2_samples <- replicate(100000, sum(rbinom(n[i], size = 1, prob=p2s_10[SNPnum]))/n[i])
		FstSampleNums_10 <- cbind(FstSampleNums_10,calculateFstNumerator(p_1_samples, p_2_samples,n[i],n[i]))
		FstSampleDenoms_10 <- cbind(FstSampleDenoms_10,calculateFstDenominator(p_1_samples, p_2_samples))
	}
	Fst_samples_10[[i]] <- apply(FstSampleNums_10,1,sum)/apply(FstSampleDenoms_10, 1,sum)
	Fst_samples_10[[i]][which(Fst_samples_10[[i]] < 0)] <- 0
	h <- hist(Fst_samples_10[[i]],breaks=seq(0,1,length=40),xlab="Fst estimates from 10 independent SNPs",freq=T,main=paste0("Sample Fst distribution; 10 SNPs; ",n[i]/2," diploid individuals"))
	abline(v= populationFst_5,col="red",lwd=2)
	text(0.8,max(h$counts)*0.9,paste0("population Fst: ", round(populationFst_10,3)))
	text(0.8,max(h$counts)*0.8,paste0("mean sample Fst: ", round(mean(Fst_samples_10[[i]],na.rm=T),3)))
	text(0.8,max(h$counts)*0.7,paste0("sd of sample Fst: ", round(sd(Fst_samples_10[[i]],na.rm=T),3)))
}
#dev.off()


# 6) Simulate Fst distributions when averaging over 20 SNPs: 
# Population allele frequencies:
p1s_20 <- c(0.421, 0.154, 0.055, 0.197, 0.062, 0.043, 0.094, 0.343, 0.314, 0.024, 0.35, 0.344, 0.076, 0.406, 0.325, 0.369, 0.413, 0.033, 0.242, 0.027); p2s_20 <- c(0.03, 0.423, 0.052, 0.407, 0.157, 0.038, 0.285, 0.4, 0.172, 0.352,0.013, 0.300, 0.200, 0.095, 0.425, 0.104, 0.154, 0.365, 0.184, 0.060)
populationFstNumerator_20 <- sum(calculateFstNumeratorNoNs(p1s_20,p2s_20))
populationFstDenominator_20 <- sum(calculateFstDenominator(p1s_20,p2s_20))
populationFst_20 <- populationFstNumerator_20/populationFstDenominator_20

#pdf("~/WPSGsampleFST_twentySNPs.pdf")
par(mfrow=c(3,1),mar=c(4.1,4.1,3.1,2.1))
Fst_samples_20 <- list(3)
for (i in 1:3) {
	FstSampleNums_20 <- numeric(0); FstSampleDenoms_20 <- numeric(0); 
	for (SNPnum in 1:20) {
		p_1_samples <- replicate(100000, sum(rbinom(n[i], size = 1, prob=p1s_20[SNPnum]))/n[i])
		p_2_samples <- replicate(100000, sum(rbinom(n[i], size = 1, prob=p2s_20[SNPnum]))/n[i])
		FstSampleNums_20 <- cbind(FstSampleNums_20,calculateFstNumerator(p_1_samples, p_2_samples,n[i],n[i]))
		FstSampleDenoms_20 <- cbind(FstSampleDenoms_20,calculateFstDenominator(p_1_samples, p_2_samples))
	}
	Fst_samples_20[[i]] <-apply(FstSampleNums_20,1,sum)/apply(FstSampleDenoms_20, 1,sum)
	Fst_samples_20[[i]][which(Fst_samples_20[[i]] < 0)] <- 0
	h <- hist(Fst_samples_20[[i]],breaks=seq(0,1,length=40),xlab="Fst estimates from 20 independent SNPs",freq=T,main=paste0("Sample Fst distribution; 20 SNPs; ",n[i]/2," diploid individuals"))
	abline(v= populationFst_5,col="red",lwd=2)
	text(0.8,max(h$counts)*0.9,paste0("population Fst: ", round(populationFst_20,3)))
	text(0.8,max(h$counts)*0.8,paste0("mean sample Fst: ", round(mean(Fst_samples_20[[i]],na.rm=T),3)))
	text(0.8,max(h$counts)*0.7,paste0("sd of sample Fst: ", round(sd(Fst_samples_20[[i]],na.rm=T),3)))
}
#dev.off()


#pdf("~/WPSGsampleFST_densities.pdf")
par(mfrow=c(3,1),mar=c(4.1,4.1,3.1,2.1))

for (i in 1:3) {
	if (length(which(is.na(Fst_samples[[i]]) > 0))) {
		Fst_samples[[i]] <- Fst_samples[[i]][-which(is.na(Fst_samples[[i]]))]
	}
	plot(density(Fst_samples[[i]]),xlim=c(0,1),ylim = c(0,max(c(density(Fst_samples_20[[i]])$y,density(Fst_samples[[i]])$y))),lwd=2, main=paste0("Sample Fst estimates; ",n[i]/2," diploid individuals"),xlab="Fst estimates")	
	lines(density(Fst_samples_5[[i]]),col="orange",lwd=2)
	lines(density(Fst_samples_10[[i]]),col="green",lwd=2)
	lines(density(Fst_samples_20[[i]]),col="blue",lwd=2)
	legend("topright",c("single SNP","5 SNPs","10 SNPs","20 SNPs"),lwd=2,col=c("black","orange","green","blue"))
	abline(v= populationFst_5,col="red",lwd=2)
}
#dev.off()



