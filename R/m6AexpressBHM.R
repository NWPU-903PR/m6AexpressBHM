Joint_MCMC_estimate <- function(initial_parameters_infor,gene_expre_methy,
                                size_factor,it_num,num_group,ar_lower,ar_up,
                                geweke_pvalue,prop_burn,MCMC_output_path,MCMC_output_name){
  
  #load(paste0(initial_data_path,initial_data_name))
  beta_infor <- initial_parameters_infor$coeff_infor 
  b0_infor <- initial_parameters_infor$b_infor
  initial_sigma2 <- as.numeric(as.character(b0_infor$variance_b0))
  initial_beta <- beta_infor[,1:2]
  initial_alpha <- as.numeric(as.character(beta_infor$alpha))
  initial_b0s <- b0_infor[,-ncol(b0_infor)]
  initial_b0 <- sapply(as.matrix(initial_b0s), as.numeric)
  initial_b0 <- matrix(initial_b0, nrow=nrow(initial_b0s), ncol = ncol(initial_b0s))
  colnames(initial_b0) <- colnames(initial_b0s)
  initial_b0 <- as.data.frame(initial_b0)
  #fa <-  paste0(counts_methy_path,counts_methy_name)
  gene_count_methy <-gene_expre_methy
  gene_name <- (as.character(output_data[[3]]))
  counts_methy <- data.frame()
  for (k in 1:length(gene_name)) {
    one_select <- gene_count_methy[which(!is.na(match(rownames(gene_count_methy),as.character(gene_name[k])))),]
    counts_methy <- rbind(counts_methy,one_select)
  }
  gene_reads <- counts_methy[,1:(ncol(counts_methy)/2)]
  methy_level <- counts_methy[,((ncol(counts_methy)/2)+1):(ncol(counts_methy))]
  groups <- rep((1:num_group),rep((ncol(gene_reads)/num_group),num_group))
  
  ##estimate parmater prior 
  initial_beta$beta0 <- as.numeric(as.character(initial_beta$beta0))
  initial_beta$beta1 <- as.numeric(as.character(initial_beta$beta1))
  # beta_mean <- as.numeric(as.character((round(apply(initial_beta,2,mean),2))))
  # beta_mean <- c(beta_mean[1],0)
  # beta_sd <- as.numeric(as.character(round(apply(initial_beta,2,sd),2)))
  intersecpt <- as.numeric(as.character(initial_beta$beta0))
  beta <- as.numeric(as.character(initial_beta$beta1))
  # beta0_lower <- quantile(intersecpt,0.05)
  # beta0_up <- quantile(intersecpt,0.95)
  # beta1_lower <- quantile(beta,0.05)
  # beta1_up <- quantile(beta,0.95)
  # select_beta0 <- intersecpt[(intersecpt>beta0_lower)&(intersecpt<beta0_up)]
  # select_beta1 <- beta[(beta>beta1_lower)&(beta<beta1_up)]
  # beta_mean <- c(round(mean(select_beta0),3),0)
  # beta_sd <- c(round(sd(select_beta0),3),round(sd(select_beta1),3))
  beta_mean <- c(round(mean(intersecpt),3),0)
  beta_sd <- c(round(sd(intersecpt),3),round(sd(beta)*2,3))
  ##alpha prior
  disp_fitted_value <- .disp_fitted(raw_alpha=initial_alpha,gene_reads=gene_reads,numCoef=ncol(initial_beta),numSample=ncol(gene_reads))
  disp_var_prior <- .calculate_varPrior(raw_alpha=initial_alpha,dispFitted=disp_fitted_value,numCoef=ncol(initial_beta),numSample=ncol(gene_reads))
  disp_var_prior <- (disp_var_prior)
  log_disp_fit <- log(disp_fitted_value)
  ##sigma2 prior parameter
  sigma2_prior <- .sigma2_parior_param(sigma2=initial_sigma2)
  est_alpha <- sigma2_prior[1]
  est_beta <- sigma2_prior[2]
  converage_just <- matrix(data = NA,nrow = nrow(initial_beta),ncol = (ncol(initial_beta)+4))
  ##ouput param
  betas_MCMC <- matrix(data = NA,nrow = nrow(initial_beta),ncol = ncol(initial_beta))
  alpha_MCMC <- vector()
  sigma2_MCMC <- vector()
  beta1_pvalue <- vector()
  b0s_MCMC <- matrix(data = NA,nrow = nrow(initial_b0),ncol = ncol(initial_b0))
  param_pvalue <- matrix(data = NA, nrow = nrow(initial_beta),ncol = ncol(initial_beta))
  beta1_BF <- vector()
  no_conver_gene <- vector()
  Rcpp::sourceCpp(system.file("extdata", "MCMC_process.cpp", package="m6AexpressBHM"))
  ##parameter estimation for each gene
  for (i in 1:nrow(counts_methy)) {
    y <- as.numeric(as.character(gene_reads[i,]))
    gene_methy <- as.numeric(as.character(methy_level[i,]))
    datas <- data.frame(y, gene_methy,groups)
    datas$groups <- as.factor(datas$groups)
    fixed <- y~gene_methy
    random <- ~0+groups
    kX <- model.matrix(fixed, data = datas)
    kZ <- model.matrix(random, data = datas)
    ##MC input parm
    beta_last <- as.numeric(as.character(initial_beta[i,]))
    b0_last <- as.numeric(as.character(initial_b0[i,]))
    alpha_last <- initial_alpha[i]
    sigma2_last <- initial_sigma2[i]
    log_dispersion_fit <- as.numeric(log_disp_fit[i])
    sdtune0 <- 1.0
    sdtune1 <- 1.0
    log_offset <- size_factor
    ##M-H MCMC to estimated parameter
    allparm_MCMC <- MH_MC_iter(counts=y,
                               kX,
                               kZ,
                               betas=beta_last,
                               b0s=b0_last,
                               beta_mean=beta_mean,
                               beta_sd=beta_sd,
                               log_offset = log_offset,
                               alpha=alpha_last,
                               sigma2=sigma2_last,
                               disp_variance=disp_var_prior,
                               log_dispersion_fit=log_dispersion_fit,
                               est_alpha=est_alpha,
                               est_beta=est_beta,
                               sd0=sdtune0,
                               sd1=sdtune1,
                               n_b0=num_group,
                               n_beta=2,
                               n_sample=ncol(gene_reads),
                               it_num=it_num)
    converage_justment <- .converge_justment(parm_MHMC=allparm_MCMC,it_num=it_num,ar_lower=ar_lower,
                                            ar_up=ar_up,geweke_p=geweke_pvalue,sdtune0 = sdtune0,sdtune1=sdtune1,n_beta = ncol(initial_beta),prop_burn=prop_burn)
    
    
    filed_converge <- converage_justment$failed_converge
    resample_num <- 1
    while ((filed_converge>0)&(resample_num<50)) {
      sdtune0 <- converage_justment$sdtune0
      sdtune1 <- converage_justment$sdtune1
      allparm_MCMC <- MH_MC_iter(counts=y,
                                 kX,
                                 kZ,
                                 betas=beta_last,
                                 b0s=b0_last,
                                 beta_mean=beta_mean,
                                 beta_sd=beta_sd,
                                 log_offset = log_offset,
                                 alpha=alpha_last,
                                 sigma2=sigma2_last,
                                 disp_variance=disp_var_prior,
                                 log_dispersion_fit=log_dispersion_fit,
                                 est_alpha=est_alpha,
                                 est_beta=est_beta,
                                 sd0=sdtune0,
                                 sd1=sdtune1,
                                 n_b0=num_group,
                                 n_beta=2,
                                 n_sample=ncol(gene_reads),
                                 it_num=it_num)
      converage_justment <- .converge_justment(parm_MHMC=allparm_MCMC,it_num=it_num,ar_lower=ar_lower,
                                              ar_up=ar_up,geweke_p=geweke_pvalue,sdtune0 = sdtune0,sdtune1=sdtune1,n_beta = ncol(initial_beta),prop_burn=prop_burn)
      
      filed_converge <- converage_justment$failed_converge
      resample_num <- resample_num+1
    }
    if(resample_num==50){
      no_coverage_gene <- paste0("The ",i," gene is not converage")
      print(no_coverage_gene)
      no_conver_gene[i] <- 1
    }else{
      no_conver_gene[i] <- 0
    }
    
    betas_b0s_sample <- allparm_MCMC$betas_b0s_est
    betas_sample <- (betas_b0s_sample[,1:ncol(initial_beta)])
    alpha_sample <- allparm_MCMC$alphas
    sigma2_sample <- allparm_MCMC$sigma2
    b0s_sample <- (betas_b0s_sample[,(ncol(betas_sample)+1):ncol(betas_b0s_sample)])
    ##MCMC estimated parameter value
    burn_bound <- it_num*prop_burn
    betas_MCMC[i,] <- colMedians(betas_sample[(burn_bound+1):nrow(betas_sample),],na.rm =T,hasNA = F)
    alpha_MCMC[i] <- median(alpha_sample[(burn_bound+1):length(alpha_sample)])
    sigma2_MCMC[i] <- median(sigma2_sample[(burn_bound+1):length(sigma2_sample)])
    b0s_MCMC[i,] <- colMedians(b0s_sample[(burn_bound+1):nrow(b0s_sample),],na.rm =T,hasNA = F)
    ##MCMC pvalue 
    select_betas <- betas_sample[(burn_bound+1):nrow(betas_sample),]
    select_parm <- (select_betas)
    select_alpha <- alpha_MCMC[i]
    select_b0s <- b0s_MCMC[i,] 
    param_pvalue[i,] <- .mcmc_pval(dat=select_parm,testlim=0,sided=1,ptype="z")
    # beta1_BF[i] <- BF(dat=select_betas,betas_sd = beta_sd)
    beta1_BF[i] <- .BF(readscount=y,betas=betas_MCMC[i,],alpha = select_alpha,b0s=select_b0s,kX=kX,kZ=kZ,log_offset = log_offset)
    # beta1_pvalue[i] <- mcmc_pval(dat=select_parm)
    ##converage justment result
    converage_just[i,] <- converage_justment[[1]]
    colnames(converage_just) <- c("beta0_geweke_pvalue","beta1_geweke_pvalue",
                                  "dispersion_geweke_pvalue","sigma2_geweke_pvalue",
                                  "beta_accept_rate","dispersion_accept_rate")
  }
  
  select_label <- which(no_conver_gene==0)
  colnames(betas_MCMC) <- c("beta0","beta1")
  betas_MCMC <- data.frame(gene_name=gene_name,betas_MCMC)
  betas_MCMC <- betas_MCMC[select_label,]
  colnames(b0s_MCMC) <- colnames(initial_b0)
  b0s_MCMC <- b0s_MCMC[select_label,]
  colnames(param_pvalue) <- c("beta0_pvalue","beta1_pvalue")
  param_pvalue<- data.frame(gene_name=gene_name,param_pvalue)
  param_pvalue <- param_pvalue[select_label,]
  alpha_MCMC <- alpha_MCMC[select_label]
  sigma2_MCMC <- sigma2_MCMC[select_label]
  converage_just <- converage_just[select_label,]
  beta1_BF <- beta1_BF[select_label]
  results <- list(betas_MCMC,alpha_MCMC,sigma2_MCMC,b0s_MCMC,param_pvalue,converage_just,beta1_BF)
  # results <- list(betas_MCMC,alpha_MCMC,sigma2_MCMC,b0s_MCMC,beta1_pvalue,converage_just)
  names(results) <- c("betas","dispersion","sigma2","b0s_MCMC","parameter_pvalue", "converage_result","Beta1_BF")
  save(results,file = paste0(MCMC_output_path,MCMC_output_name))
}

##parameter prior
##estimated sigma2 prior parameter given initial sigma2 value
.sigma2_parior_param <- function(sigma2){
  sigma2_lower <- quantile(sigma2,0.05)
  sigma2_uper <- quantile(sigma2,0.95)
  sigma2s <- sigma2[which((sigma2>sigma2_lower)&(sigma2<sigma2_uper))]
  mu_sigma2 <- mean(sigma2s) 
  est_var <- sum((sigma2s-mu_sigma2)^2)/(length(sigma2s)-1)
  est_alpha <- ((mu_sigma2^2)/est_var)+2
  est_beta <- mu_sigma2*(((mu_sigma2^2)/est_var)+1)
  prior_parm <- c(est_alpha,est_beta)
  return(prior_parm)
}
####estimated alpha prior by DEseq2 method given initial alpha value
.disp_fitted <- function(raw_alpha,gene_reads,numCoef,numSample){
  mean_count <- apply(gene_reads,1,mean)
  inve_meancount <- 1/mean_count
  alpha_lower <- quantile(raw_alpha,0.1)
  alpha_up <- quantile(raw_alpha,0.9)
  select_label <- which((raw_alpha>alpha_lower)&(raw_alpha<alpha_up))
  reg_data <- cbind(raw_alpha,inve_meancount)
  reg_data <-as.data.frame(reg_data)[select_label,]
  gamma_reg <- glm(raw_alpha~inve_meancount ,data = reg_data,family =  Gamma(link = "identity"))
  lambda <- as.numeric(gamma_reg$coefficients)
  alpha_fit <- (lambda[2]/mean_count)+lambda[1]
  return(alpha_fit)
}

.calculate_varPrior <- function(raw_alpha, dispFitted, numSample,numCoef){
  varlogdisp <- trigamma((numSample-numCoef)/2)
  logResidule = log(raw_alpha) - log(dispFitted)
  stdLogResidule = median(abs(logResidule - median(logResidule))) * 1.4826
  varLogResidule = stdLogResidule ** 2
  varPrior = varLogResidule - varlogdisp
  varPrior = max(varPrior, 0.25)
  return(varPrior)   
}
##converge justment
.converge_justment <- function(parm_MHMC,it_num,ar_lower,ar_up,geweke_p,sdtune0, sdtune1,n_beta,prop_burn){
  betas_b0s_sample <- parm_MHMC$betas_b0s_est
  betas_sample <- (betas_b0s_sample[(it_num*prop_burn+1):it_num,1:n_beta])
  alpha_sample <- parm_MHMC$alphas[(it_num*prop_burn+1):it_num]
  sigma2_sample <- parm_MHMC$sigma2[(it_num*prop_burn+1):it_num]
  b0s_sample <- (betas_b0s_sample[(it_num*prop_burn+1):it_num,(ncol(betas_sample)+1):ncol(betas_b0s_sample)])

  ####geweke pvalue
  beta_geweke_pvalue <- pnorm(abs(geweke.diag(mcmc(betas_sample),frac1 = 0.2)$z),lower.tail=FALSE)*2
  alpha_geweke_pvalue <- pnorm(abs(geweke.diag(mcmc(alpha_sample),frac1 = 0.2)$z),lower.tail=FALSE)*2
  sigma2_geweke_pvalue <- pnorm(abs(geweke.diag(mcmc(sigma2_sample),frac1 = 0.2)$z),lower.tail=FALSE)*2
  beta_accept_rate <- min(length(unique(betas_sample[,1]))/it_num,length(unique(betas_sample[,2]))/it_num)
  alpha_accept_rate <- length(unique(alpha_sample))/it_num
  ##geweke check 
  geweke_all <- c(beta_geweke_pvalue, alpha_geweke_pvalue,sigma2_geweke_pvalue)
  if(length(which(is.na(geweke_all)))>0){
    geweke_all[which(is.na(geweke_all))] <- rep(0,length(which(is.na(geweke_all))))
  }
  if(min(geweke_all)<geweke_p){
    failed_geweke <- 1
  }else{
    failed_geweke <- 0
  }
  failed_ar <- 0
  if (beta_accept_rate < ar_lower){
    failed_ar <- 1
    sdtune0 <- 0.8 * sdtune0
  }
  
  if (beta_accept_rate > ar_up){
    failed_ar <- 1
    sdtune0 <- 1.2 * sdtune0
  }
  
  if (alpha_accept_rate < ar_lower){
    sdtune1 <- 0.8 * sdtune1
    failed_ar <- 1
  }
  
  if (alpha_accept_rate > ar_up){
    failed_ar <- 1
    sdtune1 <- 1.2 * sdtune1
  }
  
  failed_converge <- sum(failed_geweke,failed_ar)
  justment_result <- c(beta_geweke_pvalue,alpha_geweke_pvalue,sigma2_geweke_pvalue,
                       beta_accept_rate,alpha_accept_rate)
  names(justment_result) <- c("beta0_geweke_pvalue","beta1_geweke_pvalue",
                              "dispersion_geweke_pvalue","sigma2_geweke_pvalue",
                              "beta_accept_rate", "dispersion_accept_rate")
  converage_result<- list(justment_result,failed_converge,sdtune0,sdtune1)
  names(converage_result) <- c("justment_result","failed_converge","sdtune0", "sdtune1")
  return(converage_result)
}
#####MCMC p-value
.mcmc_pval <-
  function(dat,testlim=0,sided=1,ptype="z"){
    dat=data.frame(dat)
    lim=1/length(dat[,1])
    bps=c()
    for (i in c(1:length(names(dat)))){
      ss=dat[,i]-testlim
      mm=mean(ss)
      if(ptype=="z"){
        zs=mean(ss)/sd(ss)
        tst=sided*(1-pnorm(abs(zs)))
      }
      else {
        if (sided==2) tst=sum(mm*ss<0)
        if (sided==1) tst=sum(ss<0)
        if (tst==0) tst=sided*lim else tst=sided*tst/length(ss)
        #		tst=as.numeric(tst)
      }
      bps=append(bps,tst)
    }
    return(bps)
  }
##Bayese factor BF
.BF<- function(readscount,betas,alpha,b0s,kX,kZ,log_offset){
  M1_fixed_terms <- as.matrix(kX)%*%as.numeric(t(betas))
  M0_fixed_terms <- as.matrix(kX)%*%c(betas[1],0)
  random_terms <- as.matrix(kZ)%*%as.numeric(t(b0s))
  M1_mu <- as.numeric(as.character(exp(log_offset)*round(as.numeric(exp(M1_fixed_terms+random_terms)),3)))
  M0_mu <- as.numeric(as.character(exp(log_offset)*round(as.numeric(exp(M0_fixed_terms+random_terms)),3)))
  M1_ll <- sum(dnbinom(readscount,mu=M1_mu,size=1/alpha,log=F))
  M0_ll <- sum(dnbinom(readscount,mu=M0_mu,size=1/alpha,log=F))
  BFs <- (M1_ll)/(M0_ll)
  return(BFs)
}
