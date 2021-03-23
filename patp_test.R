patp_test <- function(data, tmat, cid = "cid", id = "id", group = "group", h = 1,
                      j = 2, s = 0, weighted=FALSE, LMAJ=FALSE, B = 1000, ipw = 0,
                      method = "linear", tau = NULL){
  n <- length(unique(data[,cid]))
  groups <- unique(data[,group])
  groups <- sort(groups[!is.na(groups)])
  
  Pt <- list()
  tms <- list()
  S = seq(nrow(tmat))
  T_c = seq(nrow(tmat) - 1)
  
  sop_t <- function(data, tau=NULL, S, T_c, ipw=0, trans,
                    times= NULL){
    if(is.null(tau)){
      tau <- max(data$Tstop)
    }
    if(is.null(times)){
      
    }
    CTI <- list()
    counter <- 1
    for(h in T_c){
      for(j in S[trans[h,]]){
        data_h <- data[data$from==h,]
        data_h$delta <- 1*(data_h$to==j)
        if(ipw==0){
          fit <- coxph(Surv(Tstart,Tstop,delta,type="counting")~1, data=data_h,
                       control = coxph.control(timefix = FALSE))
        } else {
          fit <- coxph(Surv(Tstart,Tstop,delta,type="counting")~1,
                       weight=weightVL, data=data_h,
                       control = coxph.control(timefix = FALSE))
        }
        A <- basehaz(fit, centered=FALSE)
        A_t<-stepfun(A$time,c(0,A$hazard))
        CTI[[counter]] <- A_t
        if(counter==1){
          pointer <- c(h,j,counter)
        } else {
          pointer <- rbind(pointer,c(h,j,counter))
        }
        counter <- counter + 1
      }
    }
    
    tt<-sort(unique(data[data$to != data$from,"Tstop"]))
    tt <- tt[tt<=tau]
    
    dA <- sapply(seq_along(CTI), function(i) {
      diff(c(0, CTI[[i]](tt)), lag = 1)
    })
    
    ttrans <- t(trans)
    mat0 <- matrix(0, nrow = nrow(trans), ncol = ncol(trans))
    mat_list <- lapply(seq_len(nrow(dA)), function(i) {
      out <- mat0
      out[which(ttrans)] <- dA[i, ]
      out <- t(out)
      ## compute diagonal elements
      diag(out) <- - rowSums(out)
      out + diag(nrow(trans))
    })
    
    #Calculate the product integral estimator
    P_n<-Reduce("%*%",  mat_list, accumulate = TRUE)
    p_0 <- sapply(S, function(i){sum(data[data$Tstart==0,"from"] == i)/nrow(data[data$Tstart==0,])})
    p_n <- matrix(NA,nrow=length(tt),ncol= nrow(trans))
    for(i in 1:length(tt)){
      p_n[i,] <- p_0%*%P_n[[i]]
    }
    p_n <- rbind(c(0,p_0),
                 cbind(tt,p_n))
    p_n <- as.data.frame(p_n)
    colnames(p_n) <- c("t",paste("p",S,sep=""))
    rownames(p_n) <- 1:(length(tt)+1)
    if (!is.null(times)){
      p_nt <- sapply(times, function(t){tail(p_n[which(p_n$t < t), ], 1)})
      colnames(p_nt) <- times
      res <- t(p_nt)
      out <- do.call(c, res)
      attributes(out) <- attributes(res)
    }else{
      out <- p_n
    }
    return(out)
  }
  
  for(g in groups){
    if (!LMAJ) {
      group_data <- data[data[,group]==g,]
      P0 <- sop_t(data = group_data, tau=NULL, S = seq(nrow(tmat)), T_c = seq(nrow(tmat) - 1),
                  ipw=0, trans = tmat, times = NULL)
      
    }
    rep1 <- function(x){
      x[c(1, seq_along(x))]
    }
    if(length(Pt)==0){
      Pt[[1]] <- stepfun(P0$t, rep1(P0[, paste0("p", j)]))
      tms[[1]] <- P0$t[diff(rep1(P0[, paste0("p", j)]), lag=1)!=0]
      
    } else {
      Pt[[2]] <- stepfun(P0$t, rep1(P0[, paste0("p", j)]))
      tms[[2]] <- P0$t[diff(rep1(P0[, paste0("p", j)]), lag=1)!=0]
    }
  }
  
  tms <- sort(unique(c(tms[[1]], tms[[2]])))
  
  if(nrow(data[data$from==j,])>0){
    tS <- sort(unique(c(data[data$to==j,"from"],j)))
  } else {
    tS <- sort(unique(data[data$to==j,"from"]))
  }
  
  EY <- NULL
  EY.t <- function(t,dt){
    sum(dt$Tstart<t & dt$Tstop>=t)/n
  }
  for(i in tS){
    for(g in groups){
      dat <- unique(data[data$from==i & data[,group]==g,
                         c(id, "from", "Tstart", "Tstop")])
      if(is.null(EY)){
        EY <- sapply(tms, EY.t, dt=dat)
      } else {
        EY <- cbind(EY, sapply(tms, EY.t, dt=dat))
      }
    }
  }
  
  Wt <- rowProds(EY)/rowSums(EY)
  
  tms <- tms[!is.na(Wt)]
  if(length(tms)==0 | max(Wt, na.rm=TRUE)==0){
    stop("Weights NA or 0 for all timepoints")
  }
  Wt <- Wt[!is.na(Wt)]
  # D_hat=(Pt[[1]](tms) - Pt[[2]](tms))
  
  estimator <- function(data,cov,tau,S,T_c,ipw, trans, Wt, method){
    if (cov != 0 && cov != 1) {
      stop("The 'cov' has to be 0 or 1.")
    }
    
    P1 <- sop_t(data[cov==1,], tau=tau, S=S, T_c=T_c, ipw=ipw, trans = tmat)
    P0 <- sop_t(data[cov==0,], tau=tau, S=S, T_c=T_c, ipw=ipw, trans = tmat)
    
    # dAUC <- vector(length = length(S))
    res <- vector(length = length(S))
    for(j in S){
      p1 <- P1[,c("t",paste("p",j,sep=""))]
      p0 <- P0[,c("t",paste("p",j,sep=""))]
      
      tau0 <- max(data$Tstop)
      
      tmax <- min(max(p1[p1[,2]>0, "t"]),
                  max(p0[p0[,2] >0, "t"]))
      
      tau <- min(tau0, tmax)
      
      dm<-diff(c(tms,tau),lag=1)
      elemenTstart<-function(t){
        max(1:length(p1$t)*(p1$t<=t))
      }
      elementfrom<-sapply(tms,elemenTstart)
      elementfrom<-(elementfrom==0)+(elementfrom>0)*elementfrom
      
      element0<-function(t){
        max(1:length(p0$t)*(p0$t<=t))
      }
      elements0<-sapply(tms,element0)
      elements0<-(elements0==0)+(elements0>0)*elements0
      
      D_t <- p1[elementfrom,paste("p",j,sep="")]-p0[elements0,paste("p",j,sep="")]
      if (method == "linear"){
        res[j] <- sum(Wt * D_t * dm)
      }
    }
    return(res)
  }
  
  if (method == "linear") {
    
    Z0 <-  estimator( data = data, cov = data$group, tau = NULL, S = S,
                      T_c = T_c, ipw = 0, trans = tmat, Wt = Wt, method = "linear")
    dauc_boot <- function (data, cov, tau =NULL, S, T_c, ipw, trans,
                           Wt, B,  verbose = 0, id)
    {
      ids <- unique(data[[id]])
      n <- length(ids)
      res <- matrix(NA, length(Z0), B)
      for (b in 1:B) {
        if (verbose > 0) {
          cat("\nBootstrap replication", b, "\n")
          flush.console()
        }
        bootdata <- NULL
        bids <- sample(ids, replace = TRUE)
        bidxs <- unlist(sapply(bids, function(x) which(x == data[[id]])))
        bootdata <- data[bidxs, ]
        if (verbose > 0) {
          print(date())
          print(events(bootdata))
          cat("applying theta ...")
        }
        if (all(bootdata$from == bootdata$to)){
          next
        }
        thstar <-  try(estimator( data = bootdata, cov = bootdata$group, tau = tau, S = S,
                                  T_c = T_c, ipw = ipw, trans = tmat, Wt = Wt, method = "linear"))
        if (class(thstar)!="try-error"){
          res[, b] <- thstar
        }
        
      }
      if (verbose)
        cat("\n")
      return(res)
    }
    dauc_res <- dauc_boot( data = data, cov = data$group, S = S, T_c = T_c, 
                           ipw = ipw, trans = tmat, Wt = Wt, B = B, id = cid)
    dauc_res <- as.numeric(na.omit(dauc_res))
    dauc_sd <- sd(dauc_res)
    
    T_linear <- Z0/dauc_sd
    pval_lin <- 2*pnorm(abs(T_linear), lower.tail = FALSE)
    pval <- pval_lin[j]
    names(pval) <- paste0("p-value at State", j)
    return(pval)
  }
}

#######################################
#### Add group information#############
#######################################

add_group <- function(data){
  group <- by(data, data$cid, function(sub_dat){
    l <- length(unique(sub_dat$id))
    rbinom(l, 1, 0.5) 
  })
  g <- do.call(c, group)
  ## Add group information
  res <- do.call(rbind, by(data, data$id, function(sub_dat){
    sub_dat$group <- rep(g[sub_dat$id[1]], nrow(sub_dat))
    return(sub_dat)
  }))
  return(res)
}

trans <- function(nstate, state_names, from, to) {
  if (missing(nstate) && missing(state_names))
    stop("One of 'nstate' and 'state_names' has to be specified.")
  if (missing(state_names)) {
    state_names <- as.character(seq_len(nstate))
  } else {
    state_names <- unique(state_names)
    nstate <- length(state_names)
  }
  if (length(from) != length(to))
    stop("The length of 'from' and 'to' must be the same.")
  if (is.character(from)) {
    from <- match(from, state_names)
  } else {
    from <- as.integer(from)
  }
  if (is.character(to)) {
    to <- match(to, state_names)
  } else {
    to <- as.integer(to)
  }
  mat <- matrix(FALSE, ncol = nstate, nrow = nstate)
  dimnames(mat) <- list(state_names, state_names)
  mat[cbind(from, to)] <- TRUE
  mat
}


get_data <- function(n, M0 = 10, M1 = 20, tmat){
  dat <- simulate(n = n, M0 = M0, M1 = M1)
  dat <- reshape_long(dat, tmat = tmat)
  dat <- add_group(dat)
  return(dat)
}

