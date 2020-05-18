#Elements utilises
Xn=c(3, 2, 2)
p = 0.1
nb_iterations = 50
rendunitaire = 1
delta = 0.87
k = 2

Modele <- function(p, Xn, nb_iterations){
  nb_tranches <- length(Xn)
  S <- sum(Xn)
  A <- as.matrix(cbind(c(1, rep(0, nb_tranches-1)), rbind(diag(1, nrow=nb_tranches -1, ncol=nb_tranches - 1), rep(0, nb_tranches-1))))
  Modelisation <- Xn
  for( i in 1:nb_iterations){
    B <- as.matrix(c(rep(0, nb_tranches-1),1))
    Xn <- A %*% as.matrix(Xn)
    if(S != 0){
      proba <- Xn/S
      Bernouilli <- sample(c(1,0), 1, prob = c(p, 1-p))
      if(Bernouilli == 1){
        lieu <- sample(1:nb_tranches, 1, prob = proba)
        B[lieu] <- B[lieu] - 1
        Xn <- Xn + B
      }
    }
    Modelisation <- cbind(Modelisation, Xn)
  }
  Modelisation <- as.matrix(Modelisation)
  return(Xn)
}

PolitiqueDurable <- function(p, Xn, nb_iterations, rendunitaire, k, delta){
  nb_tranches <- length(Xn)
  Xn <- as.matrix(Xn)
  rendementD <- rep(NA, nb_iterations)
  for( j in 1:nb_iterations){
    nb_parcelles_coupees <- min(Xn[1],k)
    Xn <- Xn - as.matrix(c(nb_parcelles_coupees, rep(0, nb_tranches - 1)))
    Xn <- Modele(p, Xn, 1)
    rendementD[j] <- (delta^j)*sqrt(rendunitaire*nb_parcelles_coupees)
    Xn <- Xn + as.matrix(c(rep(0, nb_tranches - 1), nb_parcelles_coupees))
  }
  plot(c(1:nb_iterations),rendementD,type="b",xlab="Temps (années)",ylab="Rendement")
  RendementD_final <- sum(rendementD)
  return(RendementD_final)
}


PolitiqueGloutonne <- function(p, Xn, nb_iterations, rendunitaire, delta){
  Xn <- as.matrix(Xn)
  nb_tranches <- length(Xn)
  rendementG <- rep(NA, nb_iterations)
  for( j in 1:nb_iterations){
    nb_parcelles_coupees <- Xn[1]
    Xn <- Xn - as.matrix(c(nb_parcelles_coupees, rep(0, nb_tranches - 1)))
    Xn <- Modele(p, Xn, 1)
    rendementG[j] <- (delta^j)*sqrt(rendunitaire*nb_parcelles_coupees) 
    Xn <- Xn + as.matrix(c(rep(0, nb_tranches - 1), nb_parcelles_coupees))
  }
  plot(c(1:nb_iterations),rendementG,type="b",xlab="Temps (années)",ylab="Rendement")
  RendementG_final <- sum(rendementG)
  
  return(RendementG_final)
}

PolitiqueMixte <- function(p, Xn, nb_iterations, rendunitaire, delta){
  Xn <- as.matrix(Xn)
  nb_tranches <- length(Xn)
  rendementM <- rep(NA, nb_iterations)
  for( j in 1:nb_iterations){
    nb_parcelles_coupees <- sample(0:Xn[1],1)
    Xn <- Xn - as.matrix(c(nb_parcelles_coupees, rep(0, nb_tranches - 1)))
    Xn <- Modele(p, Xn, 1)
    rendementM[j] <- (delta^j)*sqrt(rendunitaire*nb_parcelles_coupees)
    Xn <- Xn + as.matrix(c(rep(0, nb_tranches - 1), nb_parcelles_coupees))
  }
  plot(c(1:nb_iterations),rendementM,type="b",xlab="Temps (années)",ylab="Rendement")
  RendementM_final <- sum(rendementM)
  return(RendementM_final)
}

MonteCarlo <- function(p, Xn, nb_iterations, delta, k, rendunitaire, nb_simulations){
  MCG <- rep(NA, nb_simulations)
  MCD <- rep(NA, nb_simulations)
  MCM <- rep(NA, nb_simulations)
  
  for(i in 1:nb_simulations){
    MCG[i] <- PolitiqueGloutonne(p, Xn, nb_iterations, rendunitaire, delta)
    MCD[i] <- PolitiqueDurable(p, Xn, nb_iterations, rendunitaire, k, delta)
    MCM[i] <- PolitiqueMixte(p, Xn, nb_iterations, rendunitaire, delta)
  }
  
  EspG <- (1/nb_simulations)*sum(MCG)
  EspD <- (1/nb_simulations)*sum(MCD)
  EspM <- (1/nb_simulations)*sum(MCM)
  
  varG <- sum((MCG - EspG)^2)/nb_simulations
  varD <- sum((MCD - EspD)^2)/nb_simulations
  varM <- sum((MCM - EspM)^2)/nb_simulations
  
  b1ICG <- EspG - 1.96*sqrt(varG/nb_simulations)
  b1ICD <- EspD - 1.96*sqrt(varD/nb_simulations)
  b1ICM <- EspM - 1.96*sqrt(varM/nb_simulations)
  
  b2ICG <- EspG + 1.96*sqrt(varG/nb_simulations)
  b2ICD <- EspD + 1.96*sqrt(varD/nb_simulations)
  b2ICM <- EspM + 1.96*sqrt(varM/nb_simulations)

  Resultats <- data.frame(Gloutonne = c(EspG, b1ICG, b2ICG), Durable = c(EspD, b1ICD, b2ICD), Mixte =  c(EspM, b1ICM, b2ICM), row.names = c("Estimation de Monte-Carlo", "Borne infÃ©rieure de l'IC", "Borne supérieure de l'IC"))
  resultats <- data.frame(Estimation = c(EspG, EspD, EspM), Borne_Inf_IC = c(b1ICG, b1ICD, b1ICM), Borne_Sup_IC = c(b2ICG, b2ICD, b2ICM), row.names = c("Politique Gloutonne", "Politique Durable", "Politique Mixte"))
  View(resultats)
}
