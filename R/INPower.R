INPower<-function(MAFs, betas, pow, sample.size, signif.lvl, k, span=0.5,
 binary.outcome=TRUE, multi.stage.option=NULL, tgv=NULL){
# This function uses the effect sizes and allele frequencies for a set of 
# known susceptibility SNPs and the power of detection of these SNPs from the 
# original discovery samples to obtain an estimate of the total number of underlying 
# susceptibility SNPs for athat trait and the distribution of their effect sizes. 
# The function can further use the estimated number of loci and distribution of effect 
# sizes to evaluate the power for discovery of a future GWAS study (up to three-stage). 
# MAFs : Vector of minor allele frequencies associated with the set of known loci
# betas : Vector of regression effects for the set of known loci under an additive genetic model. 
# For a continuous phenotype analyzed with linear regression model, it is assumed that the 
# outcome has been standardized so that the coefficients correspond to mean change in outcome 
# per unit of s.d. for each copy of the given allele. 
# For a binary outcome analyzed with logistic regression, the regression coefficients 
# should correspond to change in log-odds-ratio per copy of the given allele.
# pow : A vector representing the powers for the known loci  in the original studies that 
#  led to their discoveries. Note these power calculations should be carefully
#  done to avoid winner's curse (it is best to obtain effect size estimates from an independent 
# replication study) and to take into consideration all complexities of the designs of the original study. 
# If the total SNP set is obtained from a group of studies for a given trait, then the power 
# for an individual marker should reflect the probability of its detection in at least one of the studies.
# span : The parameter which controls the degree of smoothing in Loess. 
#  It specifies the fraction of SNPs that are used in local linear regression to obtain the 
#        estimated number of loci at each effect size. The default is set at 0.5, but we 
#    recommend the user to set it at a value depending on the total size of the 
#        SNP set so that about 10-20 SNPs are used for local smoothing at each effect size. 
#  The total size of the SNP set should be reasonably large (e.g. at least 20 and preferably more) 
#  for application of Loess.
# binary.outcome : logical. Is the outcome binary or continuous?
# sample.size : Sample size for a future study for which integrated power calculation is desired. 
#  For case-control studies, half of the subjects are assumed to be cases and half to be controls. 
#  It can take a vector of several sample sizes for the same study as shown in the example below. 
# signif.lvl : The required genome-wide significance level for a future study 
# multi.stage.option: This option allows to set-up design parameters for the future study if it 
#  would be done in multiple stages (up to three). The option has a list of two arguments: alpha and pi, 
#  where alpha specifies the significance level(s) used for each stage to select markers for the
# subsequent stage and pi specifies the fraction of subjects who are included in the corresponding stages. 
# The default for the option is NULL, that is, the study is assumed to be single-stage.
# tgv : An optional argument using which the user can input an estimate of the known total genetic 
# variance (TGV) of the trait that may be available from familial aggregation studies. 
# For a continuous outcome, this could be an estimate of the fraction of the total variance of the 
# trait attributed to heritability. For a binary outcome, this could logarithm of squared 
#  sibling-relative-risk that is known to approximate total genetic variance under log-normal model for risk.
# k : A vector of integer values for which the user would like to calculate probabilities of the 
#  type Pr(X>= k) to evaluate the probability of detection of at least a specified number of loci in 
# future studies. In addition, the function automatically finds nine values for "k", for which the 
# probabilities are close to 0.1 to 0.9 with an increment of 0.1.
# 
### Output
# A list of two sublists with names esdist.summary and future.study.summary is returned as results.
# esdist.summary: 1) t.n.loci: the total of estimated numbers of loci, 
#                 2) gve: the genetic variance explained (GVE) by the estimated number of loci. 
#                     Note for linear regression, the GVE is expressed as a percentage of the 
#                    total variance of the outcome, since it assumed that outcome has been standerdized.
#                      Further, if an estimate of total genetic variance (TGV) 
#                    is provided by the user, then the estimate for GVE will be automatically expressed 
#                    as a percentage of TGV.
#                 3) es.dist: estimated number of loci at each different effect size
# future.study.summary: 1) e.discov: Expected number of loci to be discovered in the future study,
#                       2) e.gve: Expected genetic variance explained and its expression is similar to GVE above.
#                       3) prob.k: A table for probability of discovering at least k loci for different values of k.
#
### Note
# The projections are only shown in the range of effect size for which the original studies had at least 1% power. 
# The loess fitting procedure, however, may include additional SNPs with smaller effect sizes for local linear smoothing.  
# The user is recommended to remove SNPs that may seem clearly outliers compared to the rest 
# in terms of their effect sizes. By default the program currently remove all SNPs with power less than 0.1% from 
# the analysis to avoid undue influence of potentially outlying observations. 

#library(mvtnorm)

# check whether lengths of MAFs, betas, and pow are the same
if ((length(MAFs)!=length(betas))||(length(MAFs)!=length(pow))) {
stop("Lengths of MAFs, betas, and pow are not consistent")}
n.x<-length(MAFs)

# lower.bound for projection set at 1%
lower.bd<-0.01

# Define effect sizes as a function of MAFs and betas
es<-betas^2*2*MAFs*(1-MAFs)
# sort the input data by es
ord.es<-order(es)
MAFs<-MAFs[ord.es]
betas<-betas[ord.es]
pow<-pow[ord.es]
es<-es[ord.es]


data.tmp<-data.frame(es, 1/pow)
attr(data.tmp,'names')<-c("es","Mk.raw")
loess.fit<-loess(Mk.raw~es,data=data.tmp,span=span, degree=1)
Mk.est<-round(loess.fit$fitted,1)


# The results are shown for SNPs with power of at least lower.bd.
# Total of estimated numbers of loci by the non-parametric method (Loess). 
total.loci.nonparm<-sum(Mk.est[pow>=lower.bd])

# GV explained by the estimated number of loci by the non-parametric method (Loess)
if (is.null(tgv)&(binary.outcome)){
    total.GVE.nonparm<-sum(Mk.est[pow>=lower.bd]*es[pow>=lower.bd])
} else if (is.null(tgv)&(binary.outcome!=1)){
    total.GVE.nonparm<-paste(round(sum(Mk.est[pow>=lower.bd]*es[pow>=lower.bd])*100,2),
            "% of the total variance of the outcome",sep="")
} else if (is.null(tgv)!=1){
    total.GVE.nonparm<-paste(round(sum(Mk.est[pow>=lower.bd]*es[pow>=lower.bd])/tgv*100,2),
            "% of the total genetic variance",sep="")
}


datas<-data.frame(MAFs,betas,pow,es, Mk.est)
#show the data with estimated number of loci only for SNPs with power bigger than lower.bd
datas<-datas[pow>=lower.bd,]         
attr(datas,'names')[c(4,5)]<-c("effect.size", "est.num.loci")
nonparm.summary<-list(total.loci.nonparm, total.GVE.nonparm, datas)
attr(nonparm.summary,'names')<-c("t.n.loci","gve","es.dist")

# A function used to calculaate power for given samples sizes and signif.lvlicance levels 
#  for multi-stage design
pow.stages<-function(sig.lvls, Ns, nci, binary){
if (length(Ns)==2){
    clev1 <- qnorm(1-sig.lvls[1]/2,0,1)
    clev2 <- qnorm(1-sig.lvls[2]/2,0,1)
    pows <- rep(0,length(nci))
    for (i in (1:length(nci))){
        if (binary){ 
            mu_vec <- c(sqrt(Ns[1]), sqrt(sum(Ns)))*sqrt(nci[i]/4)
        }else{
            mu_vec <- c(sqrt(Ns[1]), sqrt(sum(Ns)))*sqrt(nci[i])
        }
        var_mat <- matrix(c(1, sqrt(Ns[1])/sqrt(sum(Ns)),sqrt(Ns[1])/sqrt(sum(Ns)), 1),2,2)
        pows[i] <- (pmvnorm(c(clev1,-Inf),c(Inf,-clev2),mu_vec,sigma=var_mat)[[1]] +
                    pmvnorm(c(clev1,clev2),c(Inf,Inf),mu_vec,sigma=var_mat)[[1]] +
                    pmvnorm(c(-Inf,-Inf),c(-clev1,-clev2),mu_vec,sigma=var_mat)[[1]] +
                    pmvnorm(c(-Inf,clev2),c(-clev1,Inf),mu_vec,sigma=var_mat)[[1]])
    }
}else{
    clev1 <- qnorm(1-sig.lvls[1]/2,0,1)
    clev2 <- qnorm(1-sig.lvls[2]/2,0,1)
    clev3 <- qnorm(1-sig.lvls[3]/2,0,1)
    pows <- rep(0,length(nci))
    for (i in (1:length(nci))){
        if (binary){ 
            mu_vec <- c(sqrt(Ns[1]), sqrt(sum(Ns[1:2])), sqrt(sum(Ns)))*sqrt(nci[i]/4)
        }else{
            mu_vec <- c(sqrt(Ns[1]), sqrt(sum(Ns[1:2])), sqrt(sum(Ns)))*sqrt(nci[i])
        }

        var_mat <- matrix(c(1, sqrt(Ns[1])/sqrt(sum(Ns[1:2])),sqrt(Ns[1])/sqrt(sum(Ns)), 
        sqrt(Ns[1])/sqrt(sum(Ns[1:2])),1, sqrt(sum(Ns[1:2]))/sqrt(sum(Ns)),
        sqrt(Ns[1])/sqrt(sum(Ns)),sqrt(sum(Ns[1:2]))/sqrt(sum(Ns)),1),3,3)
        pows[i]<-pmvnorm(c(clev1,clev2,clev3),c(Inf, Inf, Inf),mu_vec,sigma=var_mat)[[1]]+
                  pmvnorm(c(-Inf,clev2,clev3),c(-clev1, Inf, Inf),mu_vec,sigma=var_mat)[[1]]+
                  pmvnorm(c(clev1,-Inf,clev3),c(Inf, -clev2, Inf),mu_vec,sigma=var_mat)[[1]]+
                  pmvnorm(c(-Inf, -Inf,clev3),c(-clev1, -clev2, Inf),mu_vec,sigma=var_mat)[[1]]+
                  pmvnorm(c(clev1,clev2,-Inf),c(Inf, Inf, -clev3),mu_vec,sigma=var_mat)[[1]]+
                  pmvnorm(c(-Inf,clev2,-Inf),c(-clev1, Inf, -clev3),mu_vec,sigma=var_mat)[[1]]+
                  pmvnorm(c(clev1,-Inf,-Inf),c(Inf, -clev2, -clev3),mu_vec,sigma=var_mat)[[1]]+
                  pmvnorm(c(-Inf,-Inf,-Inf),c(-clev1, -clev2, -clev3),mu_vec,sigma=var_mat)[[1]]
    }
}
return(pows)
}


# multi.stage.option<-list(al=c(10^(-5),10^(-7)),pi=c(0.2, 0.3))
# multi.stage.option<-NULL

# check whether lengths of al and pi in multi.stage.option are not the same
if (!is.null(multi.stage.option)){
    if (length(multi.stage.option[[1]])!=length(multi.stage.option[[2]])) {
      stop("Lengths of al and pi in multi.stage.option are not consistent")
    }
}


pow.pr<-matrix(0,length(sample.size),sum(pow>=lower.bd))

if (is.null(multi.stage.option)){
print(binary.outcome)
    if (binary.outcome){
        for (j in (1:length(sample.size))){
            pow.pr[j,] <- pchisq(qchisq(signif.lvl,1,lower.tail=FALSE),1,sample.size[j]/4*es[pow>lower.bd],lower.tail=FALSE)
        }
    }else{
        for (j in (1:length(sample.size))){
            pow.pr[j,] <- pchisq(qchisq(signif.lvl,1,lower.tail=FALSE),1,sample.size[j]*es[pow>=lower.bd],lower.tail=FALSE)
        }
    }
}else{
    for (j in (1:length(sample.size))){
        signif.lvls<-c(multi.stage.option[[1]],signif.lvl)
        sample.sizes<-c(multi.stage.option[[2]]*sample.size[j], (1-sum(multi.stage.option[[2]]))*sample.size[j])
        pow.pr[j,]<-pow.stages(signif.lvls,sample.sizes,es[pow>=lower.bd],binary=binary.outcome)
    }
}


# Find the probability distribution of discovering k loci to get its cumulative distribution
x.k <-max(ceiling(sum(Mk.est[pow>=lower.bd])),max(k))
Mk.est.up<-ceiling(Mk.est[pow>=lower.bd])
Mk.est.dn<-floor(Mk.est[pow>=lower.bd])
n.x<-sum(pow>=lower.bd)
Pr.Z<-matrix(0,j,x.k+1)
for (j in (1:length(sample.size))){
    bino.prs.up <- dbinom(matrix(0:x.k, x.k+1, n.x), matrix(Mk.est.up, x.k+1, n.x, byrow=TRUE), matrix(pow.pr[j,],x.k+1,n.x, byrow=TRUE))
    bino.prs.dn <- dbinom(matrix(0:x.k, x.k+1, n.x), matrix(Mk.est.dn, x.k+1, n.x, byrow=TRUE), matrix(pow.pr[j,],x.k+1,n.x, byrow=TRUE))
# interpolate two binomial distributions based on the distances of Mk.est to neiboring integers.
    bino.prs <-matrix(Mk.est[pow>=lower.bd]-Mk.est.dn, x.k+1, n.x, byrow=TRUE)*bino.prs.up + 
            (1-matrix(Mk.est[pow>=lower.bd]-Mk.est.dn, x.k+1, n.x, byrow=TRUE))*bino.prs.dn
    idx <- matrix(1:(x.k+1),x.k+1,x.k+1,byrow=TRUE)-matrix(0:x.k, x.k+1, x.k+1)+lower.tri(matrix(1, x.k+1,x.k+1))*(x.k+1)
    Pr.Zn <- bino.prs[,1]
    for (j2 in 2:n.x){
        Pr.Zn1 <- Pr.Zn
        Pr.Zn <- apply(matrix(Pr.Zn1[idx],x.k+1,x.k+1)*matrix(bino.prs[,j2],x.k+1,x.k+1)*upper.tri(matrix(1,x.k+1,x.k+1),diag=TRUE),2,sum)
    }
    Pr.Z[j,] <-  Pr.Zn
}
# cumulative probability distribution of observing at least k loci
cum.Pr.Z <- 1-apply(t(Pr.Z),2,cumsum)
cum.Pr.Z[cum.Pr.Z<0]<-0     # set minus values due to compuatational error to zero.


pld<-matrix(0,length(k)+9,0)
e.discov<-rep(0,length(sample.size))
e.gve<-rep(0,length(sample.size))
for (j in (1:length(sample.size))){
    # find k's such that the probabilities of at least k discoveries are 0.1 to 0.9 with an increment of 0.1
    tmp.mat<-matrix(cum.Pr.Z[,j],dim(Pr.Z)[2],9)-matrix(seq(0.1,0.9,by=0.1),dim(Pr.Z)[2],9,byrow=TRUE)
    k.list<-sort(apply(abs(tmp.mat),2,which.min))
    ks<-sort(c(k,k.list))
    # probalities of at least k discoveries with different sample sizes
    pld<-data.frame(pld,ks, round(cum.Pr.Z[ks,j],2))
    # expected genetic variance explained
    if (is.null(tgv)&(binary.outcome)){
        e.gve[j]<-sum(Mk.est[pow>=lower.bd]*pow.pr[j,]*es[pow>=lower.bd])
    } else if (is.null(tgv)&(binary.outcome!=1)){
        e.gve[j]<-paste(round(sum(Mk.est[pow>=lower.bd]*pow.pr[j,]*es[pow>=lower.bd])*100,2),
                  "% of the total variance of the outcome",sep="")
    } else if (is.null(tgv)!=1){
        e.gve[j]<-paste(round(sum(Mk.est[pow>=lower.bd]*pow.pr[j,]*es[pow>=lower.bd])/tgv*100,2),
                  "% of the total genetic variance",sep="")
    }
    # expected number of discoveries
    e.discov[j]<-sum(Mk.est[pow>=lower.bd]*pow.pr[j,])
}
e.discov<-data.frame(sample.size, round(e.discov,1))
attr(e.discov,"names")[2]<-"e.discov"
e.gve<-data.frame(sample.size,e.gve)
attr(pld,"names")[seq(1,length(sample.size)*2-1,by=2)]<-rep("k",length(sample.size))
attr(pld,"names")[seq(2,length(sample.size)*2,by=2)]<-paste("Pr(X>=k).n=",round(sample.size),sep="")

future.study.summary<-list(e.discov, e.gve, pld)
attr(future.study.summary,'names')<-c("e.discov", "e.gve", "prob.k")
outcome<-list(nonparm.summary,future.study.summary)
attr(outcome,'names')<-c("esdist.summary","future.study.summary")
return(outcome)
}
