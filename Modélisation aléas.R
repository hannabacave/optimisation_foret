Modele <- function(p, nb_tranches, nb_iterations){
  Xn <- as.matrix(sample(nb_tranches, nb_tranches))
  A <- as.matrix(cbind(c(1, rep(0, nb_tranches-1)), rbind(diag(1, nrow=nb_tranches -1, ncol=nb_tranches - 1), rep(0, nb_tranches-1))))
  vecteur <- seq(1:nb_tranches)
  S <- sum(Xn)
  for( i in 1:nb_iterations){
    Modele <- Xn
    B <- as.matrix(c(rep(0, nb_tranches-1),1))
    proba <- Xn/S
    Bernouilli <- sample(c(1,0), 1, prob = c(p, 1-p))
    if(Bernouilli == 0){
      Xn <- A %*% Xn
    }else{
      lieu <- sample(vecteur, 1, prob = proba)
      B[lieu] <- B[lieu] - 1
      Xn <- A %*% Xn + B
    Modele <- cbind(Modele, Xn)
  Modele <- as.matrix(Modele)
  heatmap(Modele)
    }}
}  
  