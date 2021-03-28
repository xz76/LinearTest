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