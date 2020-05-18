########## Calculs utiles pour le code

#Calcul de la décompostion d'un nombre en base 10 en base quelconque
decomposition <- function(nombre, base){
  q <- nombre
  r <- NULL
  while(q > 0){
    r <- cbind(r, q %% base)
    q <- q %/% base
    }
  return(rev(r))
}

#Codage du vecteur (on rentre le vecteur et on obtient sa localisation dans l'ensemble)
code <- function(S, N, v){
  Card <- factorial(S+N)/(factorial(N)*factorial(S))
  localisation <- NA
  for(i in 1:Card){
    if(isTRUE(all(v == X(S,N)[,i]))){
      localisation <- i
    }
  }
  return(localisation)
}

#Decodage du vecteur (on rentre la localisation d'un vecteur et on obtient le vecteur)
decode <- function(S, N, localisation){
  return(as.matrix(X(S,N)[,localisation]))
}


########## Description de l'ensemble X pour S parcelles et N+1 tranches d'arbres
#Attention à bien rentre N-1

X <- function(S, N){
  N <- N+1
  nombre <- (S+1)^(N)
  Vect <- NULL
  for(i in 1:nombre){
    dec <- as.matrix(decomposition(i, S+1))
    if(length(dec) == N){
      if(sum(dec) == S){
        Vect <- cbind(Vect,dec)
      }
    }
    if(length(dec) < N){
      taille <- N - length(dec)
      complete <- as.matrix(sample(0, taille, replace=TRUE))
      dec <- rbind(complete, dec)
      if(sum(dec) == S){
        Vect <- cbind(Vect,dec)
      }
    }  
  }
  View(as.data.frame(Vect))
}


########## Calculs des évolutions possibles pour chaque vecteur

#Cacul des évolutions possibles à partir d'un vecteur initial
evol <- function(p, X0, e){
  X0 <- as.matrix(X0)
  N <- length(X0) - 1
  S <- sum(X0)
  A <- as.matrix(cbind(c(1, rep(0, N)), rbind(diag(1, nrow=N, ncol=N), rep(0, N))))
  Proba <- NULL
  Evolution <- NULL
  Val_de_e <- NULL
  if(X0[1] >= e){
    coupe <- as.matrix(c(e, rep(0, N)))
    replante <- as.matrix(c(rep(0, N), e))
    #On coupe les parcelles et on fait évoluer
    X0 <- X0 - coupe
    X0 <- A %*% X0
    #Si on n'a pas d'incendie
    X1 <- X0 + replante
    Evolution <- cbind(Evolution, code(S,N,X1))
    Proba <- cbind(Proba, 1 - p)
    Val_de_e <- cbind(Val_de_e, e)
    #Si on a un incendie
    if( S-e != 0){
      proba <- X0/(S-e)
      for(i in 1:length(proba)){
        B <- as.matrix(c(rep(0, N), 1))
        if(proba[i] != 0){
          B[i] <- B[i] - 1
          X1 <- X0 + B
          X1 <- X1 + replante
          Evolution <- cbind(Evolution, code(S,N,X1))
          Proba <- cbind(Proba, p*proba[i])
          Val_de_e <- cbind(Val_de_e, e)
        }
      }
    }
  }
  EvolProba <- as.data.frame(rbind(Evolution, Proba, Val_de_e), row.names = c("localisation de X1 dans l'ensemble", "valeur de la proba", "valeur de e"))
  return(EvolProba)
}

#Calculs de toutes les évolutions possibles
ensemble_des_evolutions <- function(S, N, p){
  Card <- factorial(S+N)/(factorial(N)*factorial(S))
  Ensemble_des_evols_possibles <- evol(p, X(S, N)[,1], 0)
  for(i in 1:Card){
    for(e in 0:as.matrix(X(S,N)[,i])[1]){
      Ensemble_des_evols_possibles <- cbind(Ensemble_des_evols_possibles, evol(p, X(S, N)[,i], e))
    }
  }
  Ensemble_des_evols_possibles <- as.data.frame(Ensemble_des_evols_possibles[,-(1:3)])
  return(Ensemble_des_evols_possibles)
}

########## Algorithme de programmation dynamique



#Programmation dynamique
prog_dyn <- function(S, N, p, delta, h){
  Card <- factorial(S+N)/(factorial(N)*factorial(S))
V <- matrix(data=0, nrow = Card, ncol = h+1)
  maximiseur <- matrix(data=0, nrow = Card, ncol = h)
  for(i in 1:h){
    print("un tour")
    for(x in 1:Card){
      calcul <- NULL
      for(e in 0:as.matrix(X(S,N)[,x])[1]){
        evolution <- evol(p, X(S,N)[,x], e)
        somme <- rep(NA, length(evolution))
        for(j in 1:length(evolution[1,])){
          somme[j] <- V[evolution[1,j], i]*evolution[2,j]
        }
        calcul[e+1] <- sqrt(e) + delta*sum(somme)
        }
      V[x, i+1] <- max(calcul)
      maximiseur[x,h-i+1] <- which.max(calcul) - 1
    }
  }
  
  nomlignes <- rep(NA, Card)
  nomcols <- rep(NA, h)
  
  for(i in 1:Card){
    nomlignes[i] <- paste0("X[,",i,"]")
  }
  for(i in 1:h){
    nomcols[i] <- paste0("h=", i)
  }
  
  
  V <- as.data.frame(V, row.names = as.vector(nomlignes), col.names= as.vector(nomcols))
  Maximiseur <- as.data.frame(maximiseur, row.names = as.vector(nomlignes), col.names=as.vector(nomcols[1:(h-1)]))
  
  View(V)
  View(Maximiseur)
  
  return(V)
}

#Fonction valeur du vecteur (3,2,2)'
V <- prog_dyn(7, 2, 0.5, 1, 50)
plot(seq(1,h),V[24,], ylab="Fonction valeur pour le vecteur (3,2,2)' ", xlab="Années", main = "Fonction valeur pour le vecteur (3,2,2)' en fonction des années", type = 'o')






