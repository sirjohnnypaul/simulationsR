#Symulacja przebiegu choroby o tendencji epidemicznej w schemacie:
#SIR - Susceptible, Infected, and Removed | Podatne, zainfekowane i usunięte (przechorowali, są odporni na ponowne zarażenie)

#funkcja symulująca zmiany poziomu osób podatnych na zakażenie, zainfekowanych i ozdrowiałych lub zmarłych
# T - liczba okresów
# a - wskaźnik infekcji
# b - wskaźnik ozdrowienia
# N - liczba początkowych potencjalnych osób do zakeżenia 
spreadSimulation <- function(a, b, N, T) { 
  S <- rep(0, T+1)
  I <- rep(0, T+1)
  R <- rep(0, T+1)
  S[1] <- N
  I[1] <- 1
  R[1] <- 0
  for (i in 1:T) {
    #szansa na to, że osoby zostaną niezainfekowane wynosi ((1-a)^I(t)
    #każda z zainfekowanych osób zaraża z prawdopodobieństwem a zatem szansa na niezainfekowanie 
    #w danym momencie wynosi 1-a do potęgi I(t) ilości osób zakażonych w danym momencie
    #stąd poziom niezainfekowanych jest wartością rozkładu dwumianowego 
    S[i+1] <- rbinom(1, S[i], (1 - a)^I[i])
    #wartość osób w grupie R (ozdrowiałe/zmarłe) jest zatem wartością ilości osób w grupie R w aktualnym
    #momencie + wartością rozkładu dwumianowego
    R[i+1] <- R[i] + rbinom(1, I[i], b)
    #wartość ilości zainfekowanych to natomiast wartość populacji + 1 pomniejszona o Ilość osób wyzdrowiałych
    #w czasie t+1 oraz pomniejszona o wartość osób niezainfekowanych w czasie t+1
    I[i+1] <- N + 1 - R[i+1] - S[i+1]
  }
  # returns a matrix
  return(matrix(c(S, I, R), ncol = 3))
}

#Wartości początkowe
#populacja całkowita średniego miasta - 100 000 osób
N <- 100000
#ilość okresów
T <- 100  
#prawdopodobieństwo infekcji
a <- 0.0005
#prawodpodobieństwo do ozdrowienia/śmierci
b <- 0.1

#wykonanie funkcji dla wybranych wartości
Z <- spreadSimulation(a,b,N,T)
colnames(Z) <- c("S", "I", "R")

#wizualizacja przebiegu epidemii 
par(mfrow=c(3,1)) 
plot(Z[,"S"], type = "n",
     xlab = "Generation", 
     ylab = "Population")
lines(Z[,"S"], col="blue")
plot(Z[,"I"], type = "n",
     xlab = "Generation", 
     ylab = "Infected")
lines(Z[,"I"], col="red")
plot(Z[,"R"], type = "n",
     xlab = "Generation", 
     ylab = "Removed/Recovered")
lines(Z[,"R"])


#wizualizacja wielokrotnej symulacji przebiegu epidemii (bardzo zbliżone symulacje)
par(mfrow=c(1,1))
plot(Z[,"S"], type = "n",
     xlab = "Generation", 
     ylab = "POPULATION | INFECTED | REMOVED")
for (i in 1:50){
  Z <- spreadSimulation(a,b,N,T)
  colnames(Z) <- c("S", "I", "R")
  lines(Z[,"S"], col="blue")
  lines(Z[,"I"], col="red")
  lines(Z[,"R"])
}


#Wizualizacja rozwoju epidemii w populacji np. miasta - forest fire epidemic model
rm(list = ls())
infectionPlate <- function(A, i, j) {
  contacted <- 0
  
  if (i > 1) {
    if (j > 1) contacted <- contacted + (A[i-1, j-1] == 1)
    contacted <- contacted + (A[i-1, j] == 1)
    if (j < ncol(A)) contacted <- contacted + (A[i-1, j+1] == 1)
  }
  if (j > 1) contacted <- contacted + (A[i, j-1] == 1)
  contacted <- contacted + (A[i, j] == 1)
  if (j < ncol(A)) contacted <- contacted + (A[i, j+1] == 1)

  if (i < nrow(A)) {
    if (j > 1) contacted <- contacted + (A[i+1, j-1] == 1)
    contacted <- contacted + (A[i+1, j] == 1)
    if (j < ncol(A)) contacted <- contacted + (A[i+1, j+1] == 1)
  }
  return(contacted)
}

forest.fire.plot <- function(X) {
  #wykres zainfekowanych i należących do grupy removed
  for (i in 1:nrow(X)) {
    for (j in 1:ncol(X)) {
      if (X[i,j] == 1) points(i, j, col = "red", pch = 19)
      else if (X[i,j] == 0) points(i, j, col = "gray48", pch = 19)
    }
  }
}

forest.fire <- function(X, a, b, 
                        pausing = FALSE) {
  plot(c(1,nrow(X)), c(1,ncol(X)), 
       type = "n", xlab = "", ylab = "")
  forest.fire.plot(X)
  burning <- TRUE
  while (burning) {
    burning <- FALSE
    if (pausing) {
      input <- readline("hit any key to see next infections")
    }
    B <- X
    for (i in 1:nrow(X)) {
      for (j in 1:ncol(X)) {
        if (X[i, j] == 2) {
          if (runif(1) > (1 - a)^infectionPlate(X, i, j)) {
            B[i, j] <- 1
          }
        } else if (X[i, j] == 1) {
          burning <- TRUE
          if (runif(1) < b) {
            B[i, j] <- 0
          }
        }
      }
    }
    X <- B
    forest.fire.plot(X)
  }
  return(X)
}

set.seed(3)
X <- matrix(2, 21, 21)
X[11, 11] <- 1
X <- forest.fire(X, .2, .4, TRUE)
X <- matrix(2, 21, 21)
X[21,] <- 1
