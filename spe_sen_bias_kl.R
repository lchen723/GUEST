spe_sen_bias_kl = function(Theta,Strue,p,singular=F){
  SPE = 0 #non-0
  SEN = 0 #0
  
  for(i in 1:p) {
    for(j in 1:p)   {
      if(Strue[i,j]!=0 && Theta[i,j] !=0) {
        SPE = SPE + 1  }
      if(Strue[i,j]==0 && Theta[i,j] ==0) {
        SEN = SEN + 1  }
      
    }
  }
  
  if (singular==F){
    result_list <- list(
      SPE = SPE / sum((Strue!=0)*1),
      SEN = SEN / (p*p - sum((Strue!=0)*1)),
      Bias = norm(Strue-Theta,"F"),
      kl_loss = abs(log(norm(Strue)/norm(Theta))+sum(diag(Theta%*%solve(Strue)))-p))
    
  }
  else
    result_list <- list(
      SPE = SPE / sum((Strue!=0)*1),
      SEN = SEN / (p*p - sum((Strue!=0)*1)),
      Bias = norm(Strue-Theta,"F"),
      kl_loss = "non")
  
  return(result_list)
  
}
