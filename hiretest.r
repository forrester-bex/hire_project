# load HIRE
library("HIREewas")

################################################################################################
#Generate the EWAS data
################################################################################################
set.seed(05222018)
#define a function to draw samples from a Dirichlet distribution
rDirichlet <- function(alpha_vec){
num <- length(alpha_vec)

temp <- rgamma(num, shape = alpha_vec, rate = 1)
return(temp / sum(temp))
}
#number of samples
n <- 180

#assignment of discrete phenotype 

#number of controls (1/3 of total)
n1 <- (n/3)
#number of cases (2/3 of total)
n2 <- (2*(n/3))

#################################################################################################
# Generating cell lineages
#################################################################################################
m <- 10000   #number of CpG sites
K <- 3       #underlying cell type number

#methylation profiles
#assume cell type 1 and cell type 2 are from the same lineage
#cell type 1
methy1 <- rbeta(m,3,6)
#cell type 2 (from the same lineage as cell type 1)
methy2 <- methy1 + rnorm(m, sd=0.01)
ind <- sample(seq_len(m), m/5)
methy2[ind] <- rbeta(length(ind),3,6)
#cell type 3
methy3 <- rbeta(m,3,6)
mu <- cbind(methy1, methy2, methy3)

#################################################################################################
#Discrete and continuous phenotypes
#################################################################################################

#number of covariates
p <- 2
#covariates / phenotype (continuous phenotype age randomly drawn btw 20 and 50)
X <- rbind(c(rep(0, n1),rep(1, n2)), runif(n, min=20, max=50))


#set risk-CpG sites under each cell type for each phenotype
beta <- array(0, dim=c(m,K,p))

#adding phenotype effects for case/control status (phenotype 1)
# effect will be assigned as pos or neg effect with equal chance
m_common <- 10
max_signal <- 0.15
min_signal <- 0.07
signs <- sample(c(-1,1), m_common*K, replace=TRUE)
beta[seq_len(m_common),seq_len(K),1] <- signs * runif(m_common*K, min=min_signal, max=max_signal)
m_seperate <- 10
signs <- sample(c(-1,1), m_seperate*2, replace=TRUE)
beta[m_common+(seq_len(m_seperate)),seq_len(2),1] <- signs *
runif(m_seperate*2, min=min_signal, max=max_signal)
signs <- sample(c(-1,1), m_seperate, replace=TRUE)
beta[m_common+m_seperate+(seq_len(m_seperate)),K,1] <- signs *
runif(m_seperate, min=min_signal, max=max_signal)


#age (bases 21 to 40) <<<<<<<<<<< check this part
base <- 20 # should be 21?
m_common <- 10

max_signal <- 0.015
min_signal <- 0.007
signs <- sample(c(-1,1), m_common*K, replace=TRUE)
beta[base+seq_len(m_common),seq_len(K),2] <- signs *
runif(m_common*K, min=min_signal, max=max_signal)
m_seperate <- 10
signs <- sample(c(-1,1), m_seperate*2, replace=TRUE)
beta[base+m_common+seq_len(m_seperate),seq_len(2),2] <- signs *
runif(m_seperate*2, min=min_signal, max=max_signal)
signs <- sample(c(-1,1), m_seperate, replace=TRUE)
beta[base+m_common+m_seperate+seq_len(m_seperate),seq_len(K),2] <- signs *
runif(m_seperate, min=min_signal, max=max_signal)


#generate the cellular compositions

# assign dirichlet distribution based on size of K 

#ddist_3cont <- (4, 4, 2+X[2,i]/10)
#ddist_3case <- (4, 4, 5+X[2,i]/10)
#ddist_5cont <- (3, 3, 3, 3, 2+X[2,i]/10)
#ddist_5case <- (3, 3, 3, 3, 5+X[2,i]/10)
#ddist_7cont <- (1, 3, 3, 3, 2, 2, 2+X[2,i]/10)
#ddist_7case <- (1, 3, 3, 3, 2, 2, 5+X[2,i]/10)


P <- vapply(seq_len(n), function(i){
if(X[1,i]==0){ #if control
rDirichlet(c(4,4, 2+X[2,i]/10))
}else{
rDirichlet(c(4,4, 5+X[2,i]/10))
}
}, FUN.VALUE = rep(-1, 3))

#generate the observed methylation profiles
Ometh <- NULL
for(i in seq_len(n)){
utmp <- t(vapply(seq_len(m), function(j){
tmp1 <- colSums(X[ ,i] * t(beta[j, , ]))
rnorm(K,mean=mu[j, ]+tmp1,sd=0.01)
}, FUN.VALUE = rep(-1, K)))
tmp2 <- colSums(P[ ,i] * t(utmp))
Ometh <- cbind(Ometh, tmp2 + rnorm(m, sd = 0.01))
}
sum(Ometh > 1)
Ometh[Ometh > 1] <- 1
sum(Ometh < 0)
Ometh[Ometh < 0] <- 0

################################################################################################
#Apply HIRE to the simulated EWAS data
################################################################################################
#return list by HIRE
ret_list <- HIRE(Ometh, X, num_celltype=K)

#case vs control
#Visualize the association pattern with the case/control status in the first 100 CpG sites 
#write output to pdf
pdf(K+"K"+n+"n_casecontrolstatus.pdf")
riskCpGpattern(ret_list$pvalues[seq_len(100), c(2,1,3)],
main_title="Detected association pattern\n with disease status", hc_row_ind = FALSE) #c(2,1,3) was used because of the label switching
dev.off()
#age
#Visualize the association pattern with the age in the first 100 CpG sites
#write output to pdf
pdf(K+"K"+n+"n_age.pdf")
riskCpGpattern(ret_list$pvalues[seq_len(100), K+c(2,1,3)],
    main_title="Detected association pattern\n with age", hc_row_ind = FALSE)
dev.off()

#a p-value matrix from the uniform distribution
#write output to pdf
    pvalues <- matrix(runif(600), 100, 6)
    #Visualize this p-value matrix
    pdf(K+"K"+n+"n_pvalues.pdf")
    riskCpGpattern(pvalues,
    main_title="An example", hc_row_ind = FALSE)
    dev.off()