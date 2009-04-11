############################################################
## class methods
############################################################

CheckInput     <- function(x, ...) UseMethod("CheckInput")
QueueingModel  <- function(x, ...) UseMethod("QueueingModel")
Inputs         <- function(x, ...) UseMethod("Inputs")
RO             <- function(x, ...) UseMethod("RO")
Lq             <- function(x, ...) UseMethod("Lq")
Wq             <- function(x, ...) UseMethod("Wq")
L              <- function(x, ...) UseMethod("L")
W              <- function(x, ...) UseMethod("W")
WWq            <- function(x, ...) UseMethod("WWq")
Pn             <- function(x, ...) UseMethod("Pn")
LLq            <- function(x, ...) UseMethod("LLq")
Throughput     <- function(x, ...) UseMethod("Throughput")
Throughputi    <- function(x, ...) UseMethod("Throughputi")
WWs            <- function(x, ...) UseMethod("WWs")
SP             <- function(x, ...) UseMethod("SP")
Li             <- function(x, ...) UseMethod("Li")
Wi             <- function(x, ...) UseMethod("Wi")

############################################################
## Auxiliary functions
############################################################
is.anomalous <- function(x)
{
  is.null(x) || is.na(x) || is.nan(x)
}

C_erlang2 <- function(c, r)
{
  ro <- r / c

  totr <- 1
  totn <- 1
  total <- totr / totn

	i <- 1
	while (i <= c-1)
  {
		totr <- totr * r
		totn <- totn * i
		total <- total + (totr / totn)
		i <- i + 1
  }

  totr <- totr * r
  totn <- totn * c
  numerator <- totr / totn
  denominator <- (1 - ro) * (total + (numerator / (1 - ro)))
  numerator / denominator  
}

nodes <- function(...)
{
  list(...)
}


C_erlang3 <- function(c, r)
{
  b_result <- B_erlang(c, r)
  num <- c * b_result
  den <- c - (r * (1 - b_result))
  num / den    
}


# this saves one step of B_erlang, more efficient
C_erlang <- function(c, r)
{
  b_result <- B_erlang(c-1, r)
  num <- r * b_result
  den <- c - (r * (1 - b_result))
  num / den    
}



# recursive version, in R has problems of stack overflow when c is large
B_erlang2 <- function(c, u)
{

	f <- function(c, u)
	{
		if (c == 0) 1
		else ( 1 + ( f(c-1, u) * (c/u) ) )
	}
	1 / f(c, u)
}


# definition version
B_erlang3 <- function(c, u)
{
  n_fact <- 1
  u_power <- 1
  tot <- u_power / n_fact

  i <- 1
  while (i <= c)
  {
		n_fact <- i * n_fact
    u_power <- u_power * u
    tot <- tot + (u_power / n_fact)
    i <- i + 1
  }
  
  (u_power / n_fact) / tot
	
}


B_erlang <- function(c, u)
{

  tot <- 1
  aux <- 1 / u
  i <- 1

  while (i <= c)
  {
    tot <- 1 + (tot * aux)
    aux <- aux + (1 / u)
    i <- i + 1
  }

  1/tot
	
}



############################################################
## Model M/M/C
############################################################
NewInput.MMC <- function(lambda=0, mu=0, c=1, n=0)
{
  res <- list(lambda = lambda, mu = mu, c = c, n = n)
  class(res) <- "i_MMC"
  res
}

CheckInput.i_MMC <- function(x, ...)
{
	MMC_r_c_warning <- "( lambda/(mu*c) ) has to be less than one!!"
  MMC_c_warning <- "c has to be at least one!!"
  MMC_mu_positive <- "mu must be greater than zero"
  MMC_lambda_zpositive <- "lambda must be equal or greater than zero"
  MMC_class <- "the class of the object x has to be M/M/C (i_MMC)" 
  MMC_n_zpositive <- "the number of clients must be equal or greater than zero"
  MMC_anomalous <- "Some value of lambda, mu, c or n is anomalous. Check the values."


  if (class(x) != "i_MMC")
   	stop(MMC_class)

  if (is.anomalous(x$lambda) || is.anomalous(x$mu) ||
      is.anomalous(x$c) || is.anomalous(x$n)
  )
    stop(MMC_anomalous)    

  r <- x$lambda / x$mu  

	if (x$c < 1)
    stop(MMC_c_warning)

	if (x$lambda < 0)
		stop(MMC_lambda_zpositive)

	if (x$mu <= 0)
		stop(MMC_mu_positive)
  
  if (r >= x$c)
  {
    ro <- r/x$c
    print(paste("Throughput is: ", x$mu * x$c, sep=""))
    print(paste("Utilization exceeds 100% use!!: ", ro * 100, sep=""))
    stop(MMC_r_c_warning)
  }

  if (x$n < 0)
		stop(MMC_n_zpositive)

}


MMC_InitPn <- function(x)
{
    	
    r <- x$lambda / x$mu
    ro <- r / x$c
    one_minus_ro <- 1 - ro
    
    prod <- 1
  	acum <- prod  	
  	pn <- numeric()

		i <- 1
    pn[i] <- prod

		while ( i <= (x$c - 1) )
  	{
    	prod <- prod * r/i
    	acum <- acum + prod
      pn[i+1] <- prod
    	i <- i + 1
  	}

  	prod <- prod * r/x$c
    pn[x$c+1] <- prod
    
    p0 <- 1 / (acum + (prod / one_minus_ro))
      
    if (x$n > x$c)
    {
     	for (j in (x$c+1):x$n)
			{
				prod <- prod * r/x$c
        pn[j+1] <- prod
			}
    }    

    # Now, calculate the complete probabilities
    pn <- p0 * pn
    list(Pn = pn, Acum = p0 * acum)
   	  
}


QueueingModel.i_MMC <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MMC(x, ...)

  r <- x$lambda / x$mu
  RO <- r / x$c
  one_minus_ro <- 1 - RO
  inverse_lambda <- 1 / x$lambda
  Throughput <- x$lambda

  aux <- MMC_InitPn(x)
  Pn <- aux$Pn
  Lq <- ( 1 - aux$Acum ) * (RO / one_minus_ro)
  
  Wq <- Lq * inverse_lambda
  L <- Lq + r  
  W <- L * inverse_lambda
  WWq <- 1 / (x$c * one_minus_ro * x$mu)  

  FW <- function(t)
  {
    1 - ( C_erlang(x$c, r) * exp( (-1) * (1 - RO) * x$c * x$mu * t ) )  
  }

  FWq <- function(t)
  {
      
    if (r == (x$c - 1))
    {
      res <- 1 - ( 1 + C_erlang(x$c, r) * x$mu * t * exp(-x$mu * t) )
    }
    else
    {
      aux1 <- ( r - x$c + 1 - C_erlang(x$c, r) ) * exp(-x$mu * t)
      aux2 <- C_erlang(x$c, r) * exp( (-1) * (1 - RO) * x$c * x$mu * t )
      aux <- ( x$c - 1 - r ) * (aux1 + aux2)
      res <- 1 + aux
    }
    res
  }

  res <- list(
    Inputs = x, RO = RO, Lq = Lq, Wq = Wq, Throughput = Throughput, L = L, W = W, WWq = WWq,
    Pn = Pn, FW = FW, FWq = FWq
  )

  class(res) <- "o_MMC"
  res
}

Inputs.o_MMC <- function(x, ...){ x$Inputs }
RO.o_MMC <- function(x, ...){ x$RO }
Lq.o_MMC <- function(x, ...){ x$Lq }
Wq.o_MMC <- function(x, ...){ x$Wq }
L.o_MMC <- function(x, ...){ x$L }
W.o_MMC <- function(x, ...){ x$W }
WWq.o_MMC <- function(x, ...){ x$WWq } 
Pn.o_MMC <- function(x, ...){ x$Pn }
Throughput.o_MMC <- function(x, ...) { x$Throughput }

summary.o_MMC <- function(object, ...)
{ 
  Ls <- object$L - object$Lq

  print("The inputs of the model M/M/c are:")
  print(object$Inputs)
  print("", quote=FALSE)
  print("The outputs of the model M/M/c are:")
  print("", quote=FALSE)
  print(paste("The probability (p0, p1, ..., pk) of the clients in the system are:"))
  print(object$Pn)
  print("", quote=FALSE)
  print(paste("The traffic intensity is: ", Ls))
  print(paste("The server use is: ", object$RO))
  print(paste("The mean number of clients in the system is: ", object$L))
  print(paste("The mean number of clients in the queue is: ", object$Lq))
  print(paste("The mean number of clients in the server is: ", Ls))
  print(paste("The mean time spend in the system is: ", object$W))
  print(paste("The mean time spend in the queue is: ", object$Wq))
  print(paste("The mean time spend in the server is: ", object$W - object$Wq))
  #print(paste("The mean number of clients in queue when there is queue is: ", object$LLq))
  print(paste("The mean time spend in the queue when there is queue is: ", object$WWq))
  print(paste("The throughput is: ", object$Throughput))
}


############################################################
############################################################
## MODEL M/M/1
############################################################
############################################################
NewInput.MM1 <- function(lambda=0, mu=0, n=0)
{
  res <- list(lambda = lambda, mu = mu, n = n)
  class(res) <- "i_MM1"
  res
}

CheckInput.i_MM1 <- function(x, ...)
{
 MM1_ro_warning <- "ro is greater or equal to one!!"
 MM1_mu_positive <- "mu must be greater than zero"
 MM1_lambda_zpositive <- "lambda must be equal or greater than zero"
 MM1_n_zpositive <- "the number of clients must be equal or greater than zero"
 MM1_class <- "the class of the object x has to be M/M/1 (i_MM1)"
 MM1_anomalous <- "Some value of lambda, mu or n is anomalous. Check the values." 

 if (class(x) != "i_MM1")
  stop(MM1_class)

 if (is.anomalous(x$lambda) || is.anomalous(x$mu) || is.anomalous(x$n))
    stop(MM1_anomalous)

 if (x$mu <= 0)
 	stop(MM1_mu_positive)

 if (x$lambda < 0)
	stop(MM1_lambda_zpositive)

 if (x$n < 0)
  stop(MM1_n_zpositive)

 ro <- x$lambda / x$mu
	if (ro >= 1)
 {
  print(paste("Throughput is equal to: ", x$mu, sep=""))
  print(paste("Utilization is equal to: ", ro * 100, sep=""))
 	stop(MM1_ro_warning)
 }
}



QueueingModel.i_MM1 <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MM1(x, ...)

  # variables to improve the eficiency of the computing
  aux <- (x$mu - x$lambda)
  aux1 <- (x$mu * aux)
  
  RO <- x$lambda / x$mu
  Lq <- (x$lambda^2) / aux1
  Wq <- x$lambda / aux1
  L <- x$lambda / aux
  W <- 1 / aux
  LLq <- x$mu / aux
  Throughput <- x$lambda

  Pn <- numeric()
  
  #i <- 0
  #//while (i <= x$n)
  #{
  #  Pn[i+1] <- dgeom(i, 1-(x$lambda/x$mu))
  #  i <- i + 1
  #}

  Pn <- sapply(seq(0, x$n, 1), function(i){dgeom(i, 1-RO)})

  # The distribution functions
  FWq <- function(t) { 1 - (RO * exp(-t/W)) }
  FW <- function(t) { 1 - exp(-t/W) }


  res <- list(
    Inputs = x, RO = RO, Lq = Lq, Wq = Wq, Throughput = Throughput, L = L, W = W, WWq = Wq, LLq = LLq,
    Pn = Pn, FW = FW, FWq = FWq 
  )

  class(res) <- "o_MM1"
  res
} 

RO.o_MM1 <- function(x, ...){ x$RO }
Pn.o_MM1 <- function(x, ...){ x$Pn }
Lq.o_MM1 <- function(x, ...){ x$Lq }
Wq.o_MM1 <- function(x, ...){ x$Wq }
L.o_MM1 <- function(x, ...){ x$L }
W.o_MM1 <- function(x, ...){ x$W }
WWq.o_MM1 <- function(x, ...){ x$WWq }
LLq.o_MM1 <- function(x, ...){ x$LLq }
Inputs.o_MM1 <- function(x, ...){ x$Inputs }
Throughput.o_MM1 <- function(x, ...) { x$Throughput }

summary.o_MM1 <- function(object, ...)
{ 
  Ls <- object$L - object$Lq

  print("The inputs of the model M/M/1 are:")
  print(object$Inputs)
  print("", quote=FALSE)
  print("The outputs of the model M/M/1 are:")
  print("", quote=FALSE)
  print(paste("The probability (p0, p1, ..., pn) of the clients in the system are:"))
  print(object$Pn)
  print("", quote=FALSE)
  print(paste("The traffic intensity is: ", object$RO))
  print(paste("The server use is: ", object$RO))
  print(paste("The mean number of clients in the system is: ", object$L))
  print(paste("The mean number of clients in the queue is: ", object$Lq))
  print(paste("The mean number of clients in the server is: ", Ls))
  print(paste("The mean time spend in the system is: ", object$W))
  print(paste("The mean time spend in the queue is: ", object$Wq))
  print(paste("The mean time spend in the server is: ", object$W - object$Wq))
  print(paste("The mean number of clients in queue when there is queue is: ", object$LLq))
  print(paste("The mean time spend in the queue when there is queue is: ", object$WWq))
  print(paste("The throughput is: ", object$Throughput))
}


###############################################################
###############################################################
## MODEL M/M/1/K/K - Finite Poblation.                       ##
###############################################################
###############################################################
NewInput.MM1KK <- function(lambda=0, mu=0, k=1)
{
  res <- list(lambda = lambda, mu = mu, k = k)
  class(res) <- "i_MM1KK"
  res
}

CheckInput.i_MM1KK <- function(x, ...)
{

  MM1KK_lambda_zpositive <- "lambda must be equal or greater than zero"
  MM1KK_mu_positive <- "mu must be greater than zero"
  MM1KK_k_positive <- "k must be at least one"
  MM1KK_class <- "The class of the object x has to be M/M/1/K/K (i_MM1KK)"
  MM1KK_anomalous <- "Some value of lambda, mu or k is anomalous. Check the values."

 if (class(x) != "i_MM1KK")
  stop(MM1KK_class)

 if (is.anomalous(x$lambda) || is.anomalous(x$mu) || is.anomalous(x$k))
    stop(MM1KK_anomalous)

 if (x$lambda < 0)
	stop(MM1KK_lambda_zpositive)

 if (x$mu <= 0)
 	stop(MM1KK_mu_positive)

 if (x$k < 1)
 	stop(MM1KK_k_positive)
}


MM1KK_InitPn <- function(x)
{
  pn <- numeric()

  z <- x$mu / x$lambda
	u <- x$lambda / x$mu
  
  pn[1] <- B_erlang(x$k, z)
  
  totu <- 1
	totk <- 1

  i <- 2
  while (i <= x$k)
  {
    totu <- totu * u
		totk <- totk * (x$k - i + 1)
    pn[i] <- pn[1] * totu * totk 
		i <- i + 1
  }	

  pn
}


QueueingModel.i_MM1KK <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MM1KK(x, ...)
  z <- x$mu / x$lambda

  Pn <- MM1KK_InitPn(x)

  RO <- 1 - Pn[1]
  Throughput <- x$mu * RO
    
  L <- x$k - (Throughput / x$lambda)
  W <- (x$k / Throughput) - ( 1 / x$lambda)
  Wq <- W - (1 / x$mu)
  Lq <- Throughput * Wq
  WWq <- Wq / RO
  WWs <- (x$k / RO) - z
  SP <- 1 + z

  FW <- function(t){
    Qn <- function(n){ Pn[n] * (x$k - n) / (x$k - L) }
    dist <- function(n) { ppois(n, x$mu * t) }
    aux <- function(i) { Qn(i) * dist(i) }
    1 - sum(sapply(seq(1, x$k, 1), aux))
  }

  FWq <- function(t){
    Qn <- function(n){ Pn[n] * (x$k - n) / (x$k - L) }
    dist <- function(i) { ppois(i, x$mu * t) }
    aux <- function(i) { Qn(i+1) * dist(i) }
    1 - sum(sapply(seq(1, x$k-1, 1), aux))
  }


  # The result
  res <- list(
    Inputs=x, RO = RO, Lq = Lq, Wq = Wq, Throughput = Throughput, L = L, W = W, WWq = WWq, WWs = WWs, SP = SP,
    Pn = Pn, FW = FW, FWq = FWq
  )
  
  class(res) <- "o_MM1KK"
  res

} 

Inputs.o_MM1KK <- function(x, ...) { x$Inputs }
RO.o_MM1KK <- function(x, ...) { x$RO }
Lq.o_MM1KK <- function(x, ...) { x$Lq }
Wq.o_MM1KK <- function(x, ...) { x$Wq }
L.o_MM1KK <- function(x, ...) { x$L }
W.o_MM1KK <- function(x, ...) { x$W }
WWq.o_MM1KK <- function(x, ...) { x$WWq }
WWs.o_MM1KK <- function(x, ...) { x$WWs }
SP.o_MM1KK <- function(x, ...) { x$SP }
Pn.o_MM1KK <- function(x, ...) { x$RO }
Throughput.o_MM1KK <- function(x, ...) { x$Throughput }

summary.o_MM1KK <- function(object, ...)
{ 
  Ls <- object$L - object$Lq
  print("The inputs of the model M/M/1/K/K are:")
  print(object$Inputs)
  print("", quote=FALSE)
  print("The outputs of the model M/M/1/K/K are:")
  print("", quote=FALSE)
  print(paste("The probability (p0, p1, ..., pk) of the clients in the system are:"))
  print(object$Pn)
  print("", quote=FALSE)
  print(paste("The traffic intensity is: ", object$RO))
  print(paste("The server use is: ", object$RO))
  print(paste("The mean number of clients in the system is: ", object$L))
  print(paste("The mean number of clients in the queue is: ", object$Lq))
  print(paste("The mean number of clients in the server is: ", Ls))
  print(paste("The mean time spend in the system is: ", object$W))
  print(paste("The mean time spend in the queue is: ", object$Wq))
  print(paste("The mean time spend in the server is: ", object$W - object$Wq))
  #print(paste("The mean number of clients in queue when there is queue is: ", object$LLq))
  print(paste("The mean time spend in the queue when there is queue is: ", object$WWq))
  print(paste("The throughput is: ", object$Throughput))
  print(paste("The normalized average response time is: ", object$WWs))
  print(paste("The saturation point is: ", object$SP))  
}


###############################################################
###############################################################
## MODEL M/M/c/K/K - Finite Plobation, c servers        		 ##
###############################################################
###############################################################
NewInput.MMCKK <- function(lambda=0, mu=0, c=1, k=1)
{
  res <- list(lambda = lambda, mu = mu, c = c, k = k)
  class(res) <- "i_MMCKK"
  res
}


CheckInput.i_MMCKK <- function(x, ...)
{
 MMCKK_lambda_zpositive <- "lambda must be equal or greater than zero"
 MMCKK_mu_positive <- "mu must be greater than zero"
 MMCKK_k_one <- "k must be greater or equal than one"
 MMCKK_c_one <- "c must be greater or equal than one"
 MMCKK_k_c <- "k must be equal or greater than the number of servers c"
 MMCKK_class <- "The class of the object x has to be M/M/c/K/K (i_MMCKK)"
 MMCKK_anomalous <- "Some value of lambda, mu, c or k is anomalous. Check the values."

 if (class(x) != "i_MMCKK")
  stop(MMCKK_class)

 if (is.anomalous(x$lambda) || is.anomalous(x$mu) ||
      is.anomalous(x$c) || is.anomalous(x$k)
  )
    stop(MMCKK_anomalous)

 if (x$lambda < 0)
	stop(MMCKK_lambda_zpositive)

 if (x$mu <= 0)
 	stop(MMCKK_mu_positive)

 if (x$c < 1)
 	stop(MMCKK_c_one)

 if (x$k < 1)
 	stop(MMCKK_k_one)

 if (x$k < x$c)
	stop(MMCKK_k_c)
}


bak_MMCKK_InitPn <- function(x)
{
		pn <- c(0:x$k)
		fn <- c(0:x$k)

		u <- x$lambda / x$mu
		totu <- 1
		totfact <- 1
		factn <- 1
		factc <- 0
		totaux <- 1
		potc <- 1
    sumpn <- 0
		
		pn[1] <- totu * totfact
    fn[1] <- totfact
    sumpn <- sumpn + pn[1]

		i <- 1
		while (i <= x$k)
		{
			totu <- totu * u
		  factn <- factn * i
			# Factorial calculus
		  if (i <= x$k/2)
			{
				totfact <- totfact * (x$k - i + 1) / i
				fn[i+1] <- totfact
      }
			else
				totfact <- fn[x$k - i + 1]
				
			if (i == x$c) factc <- factn

			if (i > x$c)
			{
				potc <- potc * x$c
				totaux <- factn / (factc * potc)
		    pn[i+1] <- totfact * totu * totaux
        sumpn <- sumpn + pn[i+1]
			}
			else
      {
        pn[i+1] <- totfact * totu
        sumpn <- sumpn + pn[i+1]
      }
			i <- i + 1
		}
    pn/sumpn
}


MMCKK_InitPn <- function(x)
{
		pn <- rep(0, x$k+1)
		fn <- rep(0, x$k+1)

		u <- x$lambda / x$mu
		totu <- 1
		totfact <- 1
		factn <- 1
		factc <- 0
		totaux <- 1
		potc <- 1
    sumpn <- 0

    #variables to control the overflow
    qt_over <- 1e100
		
		pn[1] <- totu * totfact
    fn[1] <- totfact
    sumpn <- sumpn + pn[1]

		i <- 1
		while (i <= x$k)
		{

      # overflow control
      if (totu * u >= Inf)
      {
        totu <- (totu / qt_over) * u
        print(paste("totu: ", totu))
        pn <- pn / qt_over
        sumpn <- sumpn / qt_over
      }
      else
     	  totu <- totu * u

      # overflow control
      if (factn * i >= Inf)
      {
        factn <- (factn / qt_over) * i
        pn <- pn / qt_over
        sumpn <- sumpn / qt_over
      }
      else
        factn <- factn * i

			# Factorial calculus
		  if (i <= x$k/2)
			{
        # overflow control
        if (totfact * ((x$k - i + 1) / i) >= Inf)
        {
          totfact <- (totfact / qt_over) * ((x$k - i + 1) / i)
          print(paste("totfact: ", totfact))
          fn <- fn / qt_over
          pn <- pn / qt_over
          sumpn <- sumpn / qt_over
        }
        else
        {
				  totfact <- totfact * ((x$k - i + 1) / i)
        }
        fn[i+1] <- totfact
      }
			else
				totfact <- fn[x$k - i + 1]
				
			if (i == x$c) factc <- factn

			if (i > x$c)
			{
        if (potc * x$c >= Inf)
        {
          potc <- (potc / qt_over) * x$c
          pn <- pn * qt_over
          sumpn <- sumpn * qt_over
          print(paste("potc: ", potc))
        }
        else
				  potc <- potc * x$c

				totaux <- factn / (factc * potc)
        print(paste("totaux: ", totaux))
		    pn[i+1] <- totfact * totu * totaux
        if (pn[i+1] >= Inf | is.na(pn[i+1]))
          print(paste("Error en el paso a, : ", i))
        sumpn <- sumpn + pn[i+1]
        print(paste("sumpn: ", sumpn))
			}
			else
      {
        pn[i+1] <- totfact * totu
        if (pn[i+1] >= Inf | is.na(pn[i+1]))
          print(paste("Error en el paso a, : ", i))
        sumpn <- sumpn + pn[i+1]
        print(paste("sumpn: ", sumpn))
      }
			i <- i + 1
		}
    pn/sumpn
}



QueueingModel.i_MMCKK <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MMCKK(x, ...)

  Pn <- MMCKK_InitPn(x)

  k_per_pk <- c(0:x$k) * Pn[1:(x$k+1)]
  sum_pn_0_c_minus_1 <- sum(Pn[1:x$c])

  L <- sum(k_per_pk)

  Lq <- L - x$c - sum(k_per_pk[1:x$c]) + (x$c * sum_pn_0_c_minus_1)

  Throughput <- x$lambda * (x$k - L)
  
  W <- L / Throughput

  RO <-  Throughput / (x$c * x$mu)

  Wq <- Lq / Throughput

  WWq <- Wq / (1-sum_pn_0_c_minus_1)
  
  FW <- function(t){
    Qn <- function(n){ Pn[n] * (x$k - n) / (x$k - L) }
    dist <- function(n) { ppois(n, x$c * x$mu * t) }
    aux <- function(i) { Qn(i) * dist(i) }
    1 - sum(sapply(seq(1, x$k, 1), aux))
  }

  FWq <- function(t){
    Qn <- function(n){ Pn[n] * (x$k - n) / (x$k - L) }
    dist <- function(i) { ppois(i, x$c * x$mu * t) }
    aux <- function(i) { Qn(i+x$c) * dist(i) }
    1 - sum(sapply(seq(1, x$k-x$c, 1), aux))
  }

  # The result
  res <- list(
    Inputs=x, RO = RO, Lq = Lq, Wq = Wq, Throughput = Throughput, L = L, W = W, WWq = WWq,
    Pn = Pn, FW = FW, FWq = FWq
  )
  
  class(res) <- "o_MMCKK"
  res

} 

Inputs.o_MMCKK <- function(x, ...) { x$Inputs }
L.o_MMCKK <- function(x, ...) { x$L }
Lq.o_MMCKK <- function(x, ...) { x$Lq }
Throughput.o_MMCKK <- function(x, ...) { x$Throughput }
W.o_MMCKK <- function(x, ...) { x$W }
RO.o_MMCKK <- function(x, ...) { x$RO }
Wq.o_MMCKK <- function(x, ...) { x$Wq }
WWq.o_MMCKK <- function(x, ...) { x$WWq }
Pn.o_MMCKK <- function(x, ...) { x$Pn }

summary.o_MMCKK <- function(object, ...)
{ 
  Ls <- object$L - object$Lq

  print("The inputs of the model M/M/c/K/K are:")
  print(object$Inputs)
  print("", quote=FALSE)
  print("The outputs of the model M/M/c/K/K are:")
  print("", quote=FALSE)
  print(paste("The probability (p0, p1, ..., pk) of the clients in the system are:"))
  print(object$Pn)
  print("", quote=FALSE)
  print(paste("The traffic intensity is: ", Ls))
  print(paste("The server use is: ", object$RO))
  print(paste("The mean number of clients in the system is: ", object$L))
  print(paste("The mean number of clients in the queue is: ", object$Lq))
  print(paste("The mean number of clients in the server is: ", Ls))
  print(paste("The mean time spend in the system is: ", object$W))
  print(paste("The mean time spend in the queue is: ", object$Wq))
  print(paste("The mean time spend in the server is: ", object$W - object$Wq))
  #print(paste("The mean number of clients in queue when there is queue is: ", object$LLq))
  print(paste("The mean time spend in the queue when there is queue is: ", object$WWq))
  print(paste("The throughput is: ", object$Throughput))
}


###############################################################
###############################################################
## MODEL M/M/c/K/m - Finite Plobation, c servers, system     ##
## capacity lesser or equal than the poblation        		 	 ##
###############################################################
###############################################################
NewInput.MMCKM <- function(lambda=0, mu=0, c=1, k=1, m=1)
{
  res <- list(lambda = lambda, mu = mu, c = c, k = k, m = m)
  class(res) <- "i_MMCKM"
  res
}


CheckInput.i_MMCKM <- function(x, ...)
{

 MMCKM_lambda_zpositive <- "lambda must be equal or greater than zero"
 MMCKM_mu_positive <- "mu must be greater than zero"
 MMCKM_k_one <- "k must be greater than one"
 MMCKM_m_one <- "m must be greater than one"
 MMCKM_c_one <- "c must be greater than one"
 MMCKM_k_c <- "k must be equal or greater than the number of servers c"
 MMCKM_m_k <- "k must be equal or lesser than the poblation"
 MMCKM_class <- "The class of the object x has to be M/M/c/K/m (i_MMCKM)"
 MMCKM_anomalous <- "Some value of lambda, mu, c, k or m is anomalous. Check the values."

 if (class(x) != "i_MMCKM")
  stop(MMCKM_class)

 if (is.anomalous(x$lambda) || is.anomalous(x$mu) ||
      is.anomalous(x$c) || is.anomalous(x$k) || is.anomalous(x$m)
  )
    stop(MMCKM_anomalous)

 if (x$lambda < 0)
	stop(MMCKM_lambda_zpositive)

 if (x$mu <= 0)
 	stop(MMCKM_mu_positive)

 if (x$c < 1)
 	stop(MMCKM_c_one)

 if (x$k < 1)
 	stop(MMCKM_k_one)

 if (x$m < 1)
 	stop(MMCKM_m_one)

 if (x$k < x$c)
	stop(MMCKM_k_c)
 
 if (x$m < x$k)
	stop(MMCKM_m_k)
}

MMCKM_InitPn <- function(x)
{
		pn <- c(0:x$k)
		fn <- c(0:x$k)

		u <- x$lambda / x$mu
		totu <- 1
		totfact <- 1
		factn <- 1
		factc <- 0
		totaux <- 1
		potc <- 1
		
		pn[1] <- totu * totfact
    fn[1] <- totfact
    sum <- pn[1]

		i <- 1
		while (i <= x$k)
		{
			totu <- totu * u
		  factn <- factn * i
			# Factorial calculus
		  if (i <= x$m/2)
			{
				totfact <- totfact * (x$m - i + 1) / i
				fn[i+1] <- totfact
      }
			else
				totfact <- fn[x$m - i + 1]
		  	
			if (i == x$c) factc <- factn

			if (i > x$c)
			{
				potc <- potc * x$c
				totaux <- factn / (factc * potc)
		    pn[i+1] <- totfact * totu * totaux
        sum <- sum + pn[i+1]
			}
			else
      {
        pn[i+1] <- totfact * totu
        sum <- sum + pn[i+1]
      }		

			i <- i + 1
		}
    pn/sum
}


QueueingModel.i_MMCKM <- function(x, ...)
{
 CheckInput.i_MMCKM(x, ...)
 Pn <- MMCKM_InitPn(x)

 i_per_pn_i <- (0:x$k) * Pn[1:(x$k+1)]
 sum_pn_0_c_minus_1 <- sum(Pn[1:x$c])

 L <- sum(i_per_pn_i)
 Throughput <- x$lambda * (x$m - L)
 Lq <- L - x$c - sum(i_per_pn_i[1:x$c]) + (x$c * sum_pn_0_c_minus_1)
 W <- L / Throughput
 Wq <- Lq / Throughput 
 RO <- Throughput / (x$c * x$mu)
 WWq <- Wq / (1-sum_pn_0_c_minus_1) 

# The result
  res <- list(Inputs=x, RO = RO, Lq = Lq, Wq = Wq, Throughput = Throughput, L = L, W = W, WWq = WWq, Pn = Pn)

  class(res) <- "o_MMCKM"
  res
} 

Inputs.o_MMCKM <- function(x, ...){  x$Inputs }
L.o_MMCKM <- function(x, ...) { x$L }
Lq.o_MMCKM <- function(x, ...) { x$Lq }
Throughput.o_MMCKM <- function(x, ...) { x$Throughput }
W.o_MMCKM <- function(x, ...) { x$W }
RO.o_MMCKM <- function(x, ...) { x$RO }
Wq.o_MMCKM <- function(x, ...) { x$Wq }
WWq.o_MMCKM <- function(x, ...) { x$WWq }
Pn.o_MMCKM <- function(x, ...) { x$Pn }
summary.o_MMCKM <- function(object, ...)
{ 
  Ls <- object$L - object$Lq

  print("The inputs of the model M/M/c/K/m are:")
  print(object$Inputs)
  print("", quote=FALSE)
  print("The outputs of the model M/M/c/K/m are:")
  print("", quote=FALSE)
  print(paste("The probability (p0, p1, ..., pk) of the clients in the system are:"))
  print(object$Pn)
  print("", quote=FALSE)
  print(paste("The traffic intensity is: ", Ls))
  print(paste("The server use is: ", object$RO))
  print(paste("The mean number of clients in the system is: ", object$L))
  print(paste("The mean number of clients in the queue is: ", object$Lq))
  print(paste("The mean number of clients in the server is: ", Ls))
  print(paste("The mean time spend in the system is: ", object$W))
  print(paste("The mean time spend in the queue is: ", object$Wq))
  print(paste("The mean time spend in the server is: ", object$W - object$Wq))
  #print(paste("The mean number of clients in queue when there is queue is: ", object$LLq))
  print(paste("The mean time spend in the queue when there is queue is: ", object$WWq))
  print(paste("The throughput is: ", object$Throughput))
}


############################################################
############################################################
## MODEL M/M/Infinite
############################################################
############################################################
NewInput.MMInf <- function(lambda=0, mu=0, n=0)
{
  res <- list(lambda = lambda, mu = mu, n = n)
  class(res) <- "i_MMInf"
  res
}


CheckInput.i_MMInf <- function(x, ...)
{
  MMInf_mu_positive <- "mu must be greater than zero"
  MMInf_lambda_zpositive <- "lambda must be equal or greater than zero"
  MMInf_class <- "The class of the object x has to be M/M/Inf (i_MMInf)"
  MMInf_n_zpositive <- "the number of clients must be equal or greater than zero"
  MMInf_anomalous <- "Some value of lambda, mu, or n is anomalous. Check the values."

  if (class(x) != "i_MMInf")
   	stop(MMInf_class)	

  if (is.anomalous(x$lambda) || is.anomalous(x$mu) || is.anomalous(x$n))
    stop(MMInf_anomalous)

  if (x$mu <= 0)
 		stop(MMInf_mu_positive)

 	if (x$lambda < 0)
		stop(MMInf_lambda_zpositive)

  if (x$n < 0)
		stop(MMInf_n_zpositive)
}

QueueingModel.i_MMInf <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MMInf(x, ...)

  # Calculate the output parameters of the model 
  W <- 1 / x$mu
  L <- x$lambda * W
  
  Throughput <- x$lambda

  # we're going to calculate the probability distribution
  
  Pn <- sapply(0:x$n, dpois, L)

  FW <- function(t){ exp(x$mu) }
  FWq <- function(t){ 0 }

  # The result
  res <- list(
    Inputs=x, RO = 0, Lq = 0, Wq = 0, Throughput = Throughput, L = L, W = W, LLq = 0, WWq = 0,
    Pn = Pn, FW = FW, FWq = FWq
  )

  class(res) <- "o_MMInf"
  res

} 

Inputs.o_MMInf <- function(x, ...) { x$Inputs }
L.o_MMInf <- function(x, ...) { x$L }
W.o_MMInf <- function(x, ...) { x$W }
RO.o_MMInf <- function(x, ...) { x$RO }
Lq.o_MMInf <- function(x, ...) { x$Lq }
Wq.o_MMInf <- function(x, ...) { x$Wq }
WWq.o_MMInf <- function(x, ...) { x$WWq }
LLq.o_MMInf <- function(x, ...) { x$LLq }
Pn.o_MMInf <- function(x, ...) { x$Pn }
Throughput.o_MMInf <- function(x, ...) { x$Throughput }

summary.o_MMInf <- function(object, ...)
{ 
  Ls <- object$L - object$Lq

  print("The inputs of the model M/M/Infinite are:")
  print(object$Inputs)
  print("", quote=FALSE)
  print("The outputs of the model M/M/Infinite are:")
  print("", quote=FALSE)
  print(paste("The probability (p0, p1, ..., pn) of the clients in the system are:"))
  print(object$Pn)
  print("", quote=FALSE)
  print(paste("The traffic intensity is: ", Ls))
  print(paste("The server use is: ", object$RO))
  print(paste("The mean number of clients in the system is: ", object$L))
  print(paste("The mean number of clients in the queue is: ", object$Lq))
  print(paste("The mean number of clients in the server is: ", Ls))
  print(paste("The mean time spend in the system is: ", object$W))
  print(paste("The mean time spend in the queue is: ", object$Wq))
  print(paste("The mean time spend in the server is: ", object$W - object$Wq))
  #print(paste("The mean number of clients in queue when there is queue is: ", object$LLq))
  print(paste("The mean time spend in the queue when there is queue is: ", object$WWq))
  print(paste("The throughput is: ", object$Throughput))
}


############################################################
############################################################
## MODEL M/M/1/K - Capacity limited of the system         ##
############################################################
############################################################
NewInput.MM1K <- function(lambda=0, mu=0, k=1)
{
  res <- list(lambda = lambda, mu = mu, k = k)
  class(res) <- "i_MM1K"
  res
}

CheckInput.i_MM1K <- function(x, ...)
{
  MM1K_mu_positive <- "mu must be greater than zero"
  MM1K_lambda_zpositive <- "lambda must be equal or greater than zero"
  MM1K_k_one <- "k must be equal or greater than one"
  MM1K_class <- "the class of the object x has to be M/M/1/K (i_MM1K)"
  MM1K_anomalous <- "Some value of lambda, mu, or k is anomalous. Check the values."

  if (class(x) != "i_MM1K")
   	stop(MM1K_class)	

  if (is.anomalous(x$lambda) || is.anomalous(x$mu) || is.anomalous(x$k))
    stop(MM1K_anomalous)

  if (x$mu <= 0)
 	  stop(MM1K_mu_positive)

  if (x$lambda < 0)
	  stop(MM1K_lambda_zpositive)

  if (x$k < 1)
	  stop(MM1K_k_one) 
}


MM1K_InitPn <- function(x)
{

  pn <- numeric()
  
  if (x$lambda == x$mu)
  {
    pn[1:(x$k+1)] <- 1 / (x$k + 1)    
  }
  else
  {
    pow <- function(e, b, k){k * (b^e)}
    u <- x$lambda / x$mu
    aux <- (1 - u) / (1 - (u^(x$k+1)))
    pn <- sapply(0:x$k, pow, u, aux)
  }
	pn
}


MM1K_L <- function(x)
{
 if (x$lambda == x$mu) ( x$k / 2 )
 else
 {
		u <- x$lambda / x$mu
    u_up_k <- u^(x$k)
    u_up_k_plus_1 <- u_up_k * u
    numerator <- x$lambda * (1 - ((x$k + 1) * u_up_k) + (x$k * u_up_k_plus_1))
    denominator <- (x$mu - x$lambda) * (1 - u_up_k_plus_1)
		numerator / denominator
 }
}

QueueingModel.i_MM1K <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MM1K(x, ...)

  Pn <- MM1K_InitPn(x)

  RO <- 1 - Pn[1]

  L <- MM1K_L(x)
  Lq <- L - RO
  Throughput <- x$lambda * (1 - Pn[x$k+1])
  W <- L / Throughput
  Wq <- Lq / Throughput
  WWq <- Wq / RO
  
  FW <- function(t){
    qn <- Pn[seq(1, x$k, 1)] / (1 - Pn[x$k+1])
    dist <- function(n) { ppois(n, x$mu * t) }
    aux <- function(i) { qn[i] * dist(i) }
    1 - sum(sapply(seq(1, x$k, 1), aux))
  }

  FWq <- function(t){
    qn <- Pn[seq(1, x$k, 1)] / (1 - Pn[x$k+1])
    dist <- function(i) { ppois(i, x$mu * t) }
    aux <- function(i) { qn[i+1] * dist(i) }
    1 - sum(sapply(seq(1, x$k-1, 1), aux))
  }


  # The result
  res <- list(
    Inputs=x, RO = RO, Lq = Lq, Wq = Wq, Throughput = Throughput, L = L, W = W, WWq = WWq,
    Pn = Pn, FW = FW, FWq = FWq)

  class(res) <- "o_MM1K"
  res
} 

Inputs.o_MM1K <- function(x, ...) { x$Inputs }
L.o_MM1K <- function(x, ...) { x$L }
W.o_MM1K <- function(x, ...) { x$W }
RO.o_MM1K <- function(x, ...) { x$RO }
Lq.o_MM1K <- function(x, ...) { x$Lq }
Wq.o_MM1K <- function(x, ...) { x$Wq }
WWq.o_MM1K <- function(x, ...) { x$WWq }
Pn.o_MM1K <- function(x, ...) { x$Pn }
Throughput.o_MM1K <- function(x, ...) { x$Throughput }

summary.o_MM1K <- function(object, ...)
{ 
  Ls <- object$L - object$Lq

  print("The inputs of the model M/M/1/K are:")
  print(object$Inputs)
  print("", quote=FALSE)
  print("The outputs of the model M/M/1/K are:")
  print("", quote=FALSE)
  print(paste("The probability (p0, p1, ..., pk) of the clients in the system are:"))
  print(object$Pn)
  print("", quote=FALSE)
  print(paste("The traffic intensity is: ", object$RO))
  print(paste("The server use is: ", object$RO))
  print(paste("The mean number of clients in the system is: ", object$L))
  print(paste("The mean number of clients in the queue is: ", object$Lq))
  print(paste("The mean number of clients in the server is: ", Ls))
  print(paste("The mean time spend in the system is: ", object$W))
  print(paste("The mean time spend in the queue is: ", object$Wq))
  print(paste("The mean time spend in the server is: ", object$W - object$Wq))
  #print(paste("The mean number of clients in queue when there is queue is: ", object$LLq))
  print(paste("The mean time spend in the queue when there is queue is: ", object$WWq))
  print(paste("The throughput is: ", object$Throughput))
}


###############################################################
###############################################################
## MODEL M/M/c/K - Capacity limited of the system, c servers.##
###############################################################
###############################################################
NewInput.MMCK <- function(lambda=0, mu=0, c=1, k=1)
{
  res <- list(lambda = lambda, mu = mu, c = c, k = k)
  class(res) <- "i_MMCK"
  res
}

CheckInput.i_MMCK <- function(x, ...)
{
 MMCK_c_one <- "c has to be at least one!!"
 MMCK_lambda_zpositive <- "lambda must be equal or greater than zero"
 MMCK_k_one <- "k has to be at least one!!"
 MMCK_k_c <- "k must be equal or greater than the number of servers c"
 MMCK_mu_positive <- "mu must be greater than zero"
 MMCK_class <- "the class of the object x has to be M/M/C/K (i_MMCK)"
 MMCK_anomalous <- "Some value of lambda, mu, c or k is anomalous. Check the values."

 if (class(x) != "i_MMCK")
   	stop(MMCK_class)

 if (is.anomalous(x$lambda) || is.anomalous(x$mu) ||
      is.anomalous(x$c) || is.anomalous(x$k)
  )
    stop(MMCK_anomalous)

 if (x$lambda < 0)
	stop(MMCK_lambda_zpositive)

 if (x$mu <= 0)
 	stop(MMCK_mu_positive)

 if (x$c < 1)
	stop(MMCK_c_one)

 if (x$k < 1)
	stop(MMCK_k_one)

 if (x$k < x$c)
	stop(MMCK_k_c)
}



MMCK_InitPn <- function(x)
{
  
  pn <- numeric()

  u <- x$lambda / x$mu
  ro <- u / x$c

  prod <- 1
	acum <- 1
  aux <- 1

  i <- 1
  pn[i] <- prod # in the final, to multiply by p0

	while (i <= x$c-1)
  {
		prod <- prod * u/i
   	acum <- acum + prod
    pn[i+1] <- prod
    i <- i + 1
  }  

  prod <- prod * ro # this is the case of i = c
  pn[x$c+1] <- prod

  if (ro == 1)
		p0 <- 1 / (acum + (prod * (x$k - x$c + 1)))
	else
		p0 <- 1 / (acum + (prod * ((1 - ro^(x$k - x$c + 1)) / (1 - ro))))

 # from c+1 to k
 i <- x$c + 1

 while (i <= x$k)
 {
   prod <- prod * u/x$c
   pn[i+1] <- prod
   i <- i + 1
 }

 p0 * pn
}


QueueingModel.i_MMCK <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MMCK(x, ...)


  Pn <- MMCK_InitPn(x)
 
 	aux <- x$lambda / (x$c * x$mu)
	queue_max_length <- x$k - x$c

  Lq <-
    if (aux == 1)
		  Pn[x$c+1] * queue_max_length * (queue_max_length + 1) / 2
	  else
	  {
      one_minus_aux <- 1 - aux
		  aux_up_queue_max_length <- aux^queue_max_length
		  aux_up_queue_max_length_plus_1 <- aux_up_queue_max_length * aux
		  tmp1 <- 1 - aux_up_queue_max_length_plus_1 - ((queue_max_length+1) * aux_up_queue_max_length * one_minus_aux)
		  tmp2 <- one_minus_aux^2
		  Pn[x$c+1] * aux * (tmp1 / tmp2)
	  }

  Throughput <- x$lambda * (1 - Pn[x$k+1])

  #@@@
  # To comment to teacher's book
  L <- Lq + (Throughput / x$mu)

  RO <- Throughput / (x$mu * x$c)
  W <- L / Throughput
  Wq <- W - (1/x$mu)
  WWq <- Wq / (1 - sum(Pn[1:x$c]))

  FW <- function(t){
    qn <- Pn[seq(1, x$k, 1)] / (1 - Pn[x$k+1])
    dist <- function(n) { ppois(n, x$c * x$mu * t) }
    aux <- function(i) { qn[i] * dist(i) }
    1 - sum(sapply(seq(1, x$k, 1), aux))
  }

  FWq <- function(t){
    qn <- Pn[seq(1, x$k, 1)] / (1 - Pn[x$k+1])
    dist <- function(i) { ppois(i, x$c * x$mu * t) }
    aux <- function(i) { qn[i+x$c] * dist(i) }
    1 - sum(sapply(seq(1, x$k-x$c, 1), aux))
  }


  # The result
  res <- list(
    Inputs = x, RO = RO, Lq = Lq, Wq = Wq, Throughput = Throughput, L = L, W = W, WWq = WWq,
    Pn = Pn, FW = FW, FWq = FWq
  )

  class(res) <- "o_MMCK"
  res

} 

Inputs.o_MMCK <- function(x, ...){ x$Inputs }
L.o_MMCK <- function(x, ...) { x$L }
W.o_MMCK <- function(x, ...) { x$W }
RO.o_MMCK <- function(x, ...) { x$RO }
Lq.o_MMCK <- function(x, ...) { x$Lq }
Wq.o_MMCK <- function(x, ...) { x$Wq }
WWq.o_MMCK <- function(x, ...) { x$WWq }
Pn.o_MMCK <- function(x, ...) { x$Pn }
Throughput.o_MMCK <- function(x, ...) { x$Throughput }

summary.o_MMCK <- function(object, ...)
{ 
  Ls <- object$L - object$Lq

  print("The inputs of the model M/M/c/K are:")
  print(object$Inputs)
  print("", quote=FALSE)
  print("The outputs of the model M/M/c/K are:")
  print("", quote=FALSE)
  print(paste("The probability (p0, p1, ..., pk) of the clients in the system are:"))
  print(object$Pn)
  print("", quote=FALSE)
  print(paste("The traffic intensity is: ", Ls))
  print(paste("The server use is: ", object$RO))
  print(paste("The mean number of clients in the system is: ", object$L))
  print(paste("The mean number of clients in the queue is: ", object$Lq))
  print(paste("The mean number of clients in the server is: ", Ls))
  print(paste("The mean time spend in the system is: ", object$W))
  print(paste("The mean time spend in the queue is: ", object$Wq))
  print(paste("The mean time spend in the server is: ", object$W - object$Wq))
  #print(paste("The mean number of clients in queue when there is queue is: ", object$LLq))
  print(paste("The mean time spend in the queue when there is queue is: ", object$WWq))
  print(paste("The throughput is: ", object$Throughput))
}


###############################################################
###############################################################
## MODEL M/M/c/c - Capacity limited of the system, c servers.##
## truncated model, Erlang-B function #########################
###############################################################
###############################################################
NewInput.MMCC <- function(lambda=0, mu=0, c=1)
{
  res <- list(lambda = lambda, mu = mu, c = c)
  class(res) <- "i_MMCC"
  res
}

CheckInput.i_MMCC <- function(x, ...)
{
 MMCC_c_one <- "c has to be at least one!!"
 MMCC_lambda_zpositive <- "lambda must be equal or greater than zero"
 MMCC_mu_positive <- "mu must be greater than zero"
 MMCC_class <- "the class of the object x has to be M/M/C/C (i_MMCC)"
 MMCC_anomalous <- "Some value of lambda, mu or c is anomalous. Check the values."

 if (class(x) != "i_MMCC")
   stop(MMCC_class)

 if (is.anomalous(x$lambda) || is.anomalous(x$mu) || is.anomalous(x$c))
    stop(MMCC_anomalous)

 if (x$lambda < 0)
	stop(MMCC_lambda_zpositive)

 if (x$mu <= 0)
 	stop(MMCC_mu_positive)

 if (x$c < 1)
	stop(MMCC_c_one)
}



MMCC_InitPn <- function(x)
{
  pn <- numeric()

  u <- x$lambda / x$mu
  ro <- u / x$c
  
  prod <- 1
  acum <- 1

  i <- 1
  pn[i] <- prod

  while (i <= x$c-1)
  {
    prod <- prod * u/i
    acum <- acum + prod
    pn[i+1] <- prod    
    i <- i + 1
  }

  prod <- prod * ro
  pn[i+1] <- prod # i has the value c

  p0 <- 1 / (acum + prod)
  pn <- p0 * pn
}


QueueingModel.i_MMCC <- function(x, ...)
{
  # Is everything fine??
  CheckInput.i_MMCC(x, ...)

  Pn <- MMCC_InitPn(x)
  
  Lq <- 0
  Wq <- 0
  WWq <- 0
  LLq <- 0

  aux <- x$lambda / x$mu
  one_minus_b_erlang <- 1 - B_erlang(x$c, aux)
  L <- aux * one_minus_b_erlang
  Throughput <- x$lambda * one_minus_b_erlang
  RO <- Throughput / (x$mu * x$c)
  W <- 1 / x$mu 

  FW <- function(t){
    exp(x$mu)
  }

  FWq <- function(t){0}

  # The result
  res <- list(
    Inputs = x, RO = RO, Lq = Lq, Wq = Wq, Throughput = Throughput, L = L, W = W, LLq = LLq, WWq = WWq,
    Pn = Pn, FW = FW, FWq = FWq
  )
  
  class(res) <- "o_MMCC"
  res

} 

Inputs.o_MMCC <- function(x, ...) { x$Inputs }
L.o_MMCC <- function(x, ...) { x$L }
W.o_MMCC <- function(x, ...) { x$W }
RO.o_MMCC <- function(x, ...) { x$RO }
Lq.o_MMCC <- function(x, ...) { x$Lq }
Wq.o_MMCC <- function(x, ...) { x$Wq }
WWq.o_MMCC <- function(x, ...) { x$WWq }
Pn.o_MMCC <- function(x, ...) { x$Pn }
Throughput.o_MMCC <- function(x, ...) { x$Throughput }

summary.o_MMCC <- function(object, ...)
{
  Ls <- object$L - object$Lq
   
  print("The inputs of the model M/M/c/c are:")
  print(object$Inputs)
  print("", quote=FALSE)
  print("The outputs of the model M/M/c/c are:")
  print("", quote=FALSE)
  print(paste("The probability (p0, p1, ..., pc) of the clients in the system are:"))
  print(object$Pn)
  print("", quote=FALSE)
  print(paste("The traffic intensity is: ", Ls))
  print(paste("The server use is: ", object$RO))
  print(paste("The mean number of clients in the system is: ", object$L))
  print(paste("The mean number of clients in the queue is: ", object$Lq))
  print(paste("The mean number of clients in the server is: ", Ls))
  print(paste("The mean time spend in the system is: ", object$W))
  print(paste("The mean time spend in the queue is: ", object$Wq))
  print(paste("The mean time spend in the server is: ", object$W - object$Wq))
  #print(paste("The mean number of clients in queue when there is queue is: ", object$LLq))
  print(paste("The mean time spend in the queue when there is queue is: ", object$WWq))
  print(paste("The throughput is: ", object$Throughput))
}


###############################################################
###############################################################
## Open Jackson networks
###############################################################
###############################################################

clambda <- function(x)
{
  res <- numeric()
  i <- 1
  while (i <= length(x))
  {
    res[i] <- x[[i]]$lambda
    i <- i + 1
  }
  cbind(res)
}


newNodes <- function(rawNodes, arrivals)
{
  res <- list()
  i <- 1
  while (i <= length(rawNodes))
  {
    rawNode = rawNodes[[i]] ;

    if (class(rawNode) == "i_MM1")
      res[[i]] <- NewInput.MM1(lambda = arrivals[i], mu = rawNode$mu, n = rawNode$n)
    else if (class(rawNode) == "i_MMC")
      res[[i]] <- NewInput.MMC(lambda = arrivals[i], mu = rawNode$mu, c = rawNode$c, n = rawNode$n)
    else
      stop(paste(paste("Node ", i), "is not of class i_MM1 or i_MMC!!"))
    i <- i + 1
  }
  res
}


doModel <- function(x, newNodes, tLambda)
{
  prob <- numeric()
  Li <- numeric()
  Wi <- numeric()
  ro <- numeric()
  Throughputi <- numeric()
  is_prob_a_matrix <- (class(x$prob) == "matrix")
  totalL <- 0
  
  i <- 1
  while (i <= length(newNodes))
  {
    aux <- QueueingModel(newNodes[[i]])
    prob[i] <- Pn(aux)
 
    if (is_prob_a_matrix)
      Wi[i] <- W(aux)
    else 
      Wi[i] <- W(aux) * x$prob[i]

    auxL <- L(aux)
    Li[i] <- auxL
    ro[i] <- RO(aux)
    Throughputi[i] <- Throughput(aux)
    totalL <- totalL + auxL
    i <- i + 1
  }

  w <- totalL/tLambda
  throughput <- sum(tLambda)

  res <-
    list(
      Inputs = x, RO = ro, Pn = prob, Li = Li, L = totalL, W = w, Throughputi = Throughputi,
      Throughput = throughput, Wi = Wi
    )

  class(res) <- "o_OJN"
  res

}



CheckInput.i_OJN <- function(x, ...)
{
 x_class_OJN <- "x has to be of class i_OJN (Open Jackson Network)" 
 x_anomalous <- "x has some anomalous value. Check the value(s)."
 row_distinct_col <- "x$prob has the number of rows distinct of the number of columns."
 row_distinct_nodes <- "x$prob has distinct number of rows that the number of nodes in x$nodes."
 visit_ratios_wrong <- "x$prob contains a different number of visit ratios than x$nodes."
 prob_zero <- "If neither a routing x$prob is given nor a visit ratio vector, x$prob should be 0"
 all_lambda_equals <- "if visit ratios are given, all nodes must have the same lambda (the sum of all external arrivals)"
 
 is_prob_a_matrix <- (class(x$prob) == "matrix")

 if (is.anomalous(x$prob) || is.anomalous(x$nodes))
    stop(x_anomalous)

 if (class(x) != "i_OJN")
   stop(x_class_OJN) 

 num_nodes <- length(x$nodes)

 if (is_prob_a_matrix)
 {
   if (nrow(x$prob) != ncol(x$prob))
     stop(row_distinct_col)

   if (nrow(x$prob) != num_nodes)
     stop(row_distinct_nodes)
 }
 else
 {
   if (length(x$prob) != num_nodes)
     stop(visit_ratios_wrong)
 }

 i <- 1
 while (i <= num_nodes)
 {
   n = x$nodes[[i]]

   if (class(n) != "i_MM1" && class(n) != "i_MMC")
     stop(paste(paste("Node ", i), "is not of class i_MM1 or i_MMC!!"))
   
   if (!is_prob_a_matrix && (x$nodes[[i]]$lambda != x$nodes[[1]]$lambda))
     stop(all_lambda_equals)

   i <- i + 1
 }

}


QueueingModel.i_OJN <- function(x, ...)
{
  CheckInput(x)

  if (class(x$prob) == "matrix")
  {
    vlambda <- -clambda(x$nodes)
    tProb <- t(x$prob) 
    sol <- solve(tProb - diag(nrow=nrow(tProb)), vlambda)
    newNd <- newNodes(x$nodes, sol)  
    model <- doModel(x, newNd, -sum(vlambda))
  }
  else
  {
    lambda <- x$nodes[[1]]$lambda
    arrivals <- x$prob * lambda
    newNd <- newNodes(x$nodes, arrivals)
    model <- doModel(x, newNd, lambda)
  }
  
  model

}


NewInput.OJN <- function(prob=NULL, ...)
{
  nds <- list(prob=prob, nodes=nodes(...))
  class(nds) <- "i_OJN"
  nds
}


Inputs.o_OJN <- function(x, ...) { x$Inputs }
L.o_OJN <- function(x, ...) { x$L }
W.o_OJN <- function(x, ...) { x$W }
RO.o_OJN <- function(x, ...) { x$RO }
Li.o_OJN <- function(x, ...) { x$Li }
Pn.o_OJN <- function(x, ...) { x$Pn }
Throughput.o_OJN <- function(x, ...) { x$Throughput }
Throughputi.o_OJN <- function(x, ...) { x$Throughputi }
Wi.o_OJN <- function(x, ...) { x$Wi }


summary.o_OJN <- function(object, ...)
{
   
  print("The inputs of the open Jackson network are:")
  print(object$Inputs)
  print("", quote=FALSE)
  print("The outputs of the open Jackson network are:")
  print("", quote=FALSE)

  print("--------- Per node ---------------------------------")
 
  i <- 1
  while (i <= length(object$RO))
  {
    print(paste("The use of node ", i, " is: ", object$RO[i]))
    print(paste("The mean number of clients in node ", i, " is: ", object$Li[i]))
    print(paste("The mean time spend in node ", i, " is: ", object$Wi[i]))
    print(paste("The probability (p0, p1, ..., pn) or visit ratio of node ", i, " is: "))
    print(object$Pn[[i]])
    print(paste("The throughput of node ", i, " is: ", object$Throughputi[i]))
    i <- i + 1
  }

  print("---------- Complete network -------------------------")
  print(paste("The mean number of clients in the network is: ", object$L))
  print(paste("The mean time spend in the network is: ", object$W))
  print(paste("The throughput of the network is: ", object$Throughput))
  print("", quote=FALSE)
}

#######################################################################################
## Closed Jackson Network
#######################################################################################

NewInput.CJN <- function(prob=NULL, n=0, z=0, ...)
{
  nds <- list(prob=prob, n=n, z=z, nodes=nodes(...))
  class(nds) <- "i_CJN"
  nds
}


CheckInput.i_CJN <- function(x, ...)
{
 x_class_CJN <- "x has to be of class i_CJN (Closed Jackson Network)" 
 x_anomalous <- "x has some anomalous value. Check the value(s)."
 row_distinct_col <- "x$prob (matrix class) has the number of rows distinct of the number of columns."
 row_distinct_nodes <- "x$prob (matrix class) has distinct number of rows that the number of nodes in x$nodes."
 visit_ratios_wrong <- "x$prob contains a different number of visit ratios than x$nodes."
 n_greater_zero <- "n has to be greater than zero"

 if (is.anomalous(x$prob) || is.anomalous(x$nodes) || is.anomalous(x$n) || is.anomalous(x$z))
    stop(x_anomalous)

 if (class(x) != "i_CJN")
   stop(x_class_CJN) 

 if (x$n <= 0)
   stop(n_greater_zero)

 num_nodes <- length(x$nodes)

 is_prob_a_matrix <- (class(x$prob) == "matrix")

 if (is_prob_a_matrix)
 {
   if (nrow(x$prob) != ncol(x$prob))
     stop(row_distinct_col)

   if (nrow(x$prob) != num_nodes)
     stop(row_distinct_nodes)
 }
 else
 {
   if (length(x$prob) != num_nodes)
     stop(visit_ratios_wrong)
 }

 i <- 1
 while (i <= num_nodes)
 {
   n = x$nodes[[i]]

   if (class(n) != "i_MM1" && class(n) != "i_MMC" && class(n) != "i_MMInf")
     stop(paste(paste("Node ", i), "is not of class i_MM1 or i_MMC or i_MMInf!!"))
     
   i <- i + 1
 }

}


QueueingModel.i_CJN <- function(x, ...)
{
  num_nodes <- length(x$nodes)

  if (class(x$prob) == "matrix")
  {
    ident <- diag(dim(x$prob)[1])
    const <- matrix(data=1, nrow=dim(x$prob)[1], ncol=1)
    all1 <- matrix(data=1, ncol=dim(x$prob)[2], nrow=dim(x$prob)[1])
    prob_est <- t(solve(t(x$prob + all1 - ident), const))
  }
  else
    prob_est <- x$prob
    
  #print(paste("prob_est: ", prob_est))
  
  # count the number of nodes of type i_MMC
  num_mmc <- 0

  # create the list to hold the prob
  mclass <- list()

  k <- 1
  while (k <= num_nodes)
  {
    if (class(x$nodes[[k]]) == "i_MMC" && x$nodes[[k]]$c > 1)
      mclass <- c(mclass, list(array(0, dim=c(x$nodes[[k]]$c, x$n))))
    k <- k + 1
  }

  # array initialization
  #print(paste("num_mmc: ", num_mmc))
  #print(paste("x$n: ", x$n))

  wi <- numeric()
  li <- rep(times=num_nodes, 0)

  alfa <- function(i, c)
  {
    if (i <= c)
      res <- i
    else
      res <- c

    res
  }

  CalcProb <- function
    (n, c, mu, prob_est_mmcnode, thro, probC)
  {
    if (n == 1)
    {
      probC[, n] <- 0
      probC[1, n] <- 1
    }
    else #n>=2
    {
      sum <- 0
      j <- c
      while (j>=2)
      {
        probC[j, n] <- 
          ( (prob_est_mmcnode * thro) / (mu * alfa(j-1, c)) ) * probC[j-1, n-1]
        #print(paste("iterating inside, prob cond de ", j, " | ", n, " es: ", probC[j, n]))
        sum <- sum + ( (c - (j-1)) * probC[j, n] ) 
        j <- j - 1
      }

      sum <- (1/c) * (sum + (prob_est_mmcnode * thro / mu)) 
      probC[1, n] <- 1 - sum
      #print(paste("inside, prob cond de 1", " | ", n, " es: ", probC[1, n]))
    }
    probC
  }

  #print("Se compila la funcion, se va a entrar en el bucle")

  i <- 1
  while (i <= x$n)
  {  
    # change before putting the x$z  
    #tmp <- 0
    tmp <- x$z
    k <- 1
    while (k <= num_nodes)
    {
      num_mmc <- 0
      if (class(x$nodes[[k]]) == "i_MMInf")
        wi[k] <-  x$nodes[[k]]$mu
      else
      {
        if (class(x$nodes[[k]]) == "i_MM1" ||
            (class(x$nodes[[k]]) == "i_MMC" && x$nodes[[k]]$c == 1)
           )
          {
            #print("entrando en node mm1")
            wi[k] <- (1 + li[k]) / x$nodes[[k]]$mu
          }
        else
        {
          num_mmc <- num_mmc + 1

          #calculate probabilities
          mclass[[num_mmc]] <- CalcProb(i, x$nodes[[k]]$c,
            x$nodes[[k]]$mu, prob_est[k], throughput, mclass[[num_mmc]]
          )

          aux1 <- 1/(x$nodes[[k]]$mu * x$nodes[[k]]$c)
          sum <- 1 + li[k]
  
          z <- 1
          while (z <= x$nodes[[k]]$c - 1)
          {
            sum <- sum + 
             ( (x$nodes[[k]]$c -(z-1) -1) * (mclass[[num_mmc]][z, i]) )
            z <- z + 1
          }
          wi[k] <- aux1 * sum
        }
      }
              
      tmp <- tmp + (prob_est[k] * wi[k])
      #print(paste("wi[", k, "]: ", wi[k]))
      #print(paste("tmp: ", tmp))

      k <- k + 1
    }

    throughput <- i / tmp
    #print(paste("throughput: ", throughput))

    k <- 1
    while (k <= num_nodes)
    {
      li[k] <- throughput * prob_est[k] * wi[k]
      #print(paste("li[", k, "]: ", li[k]))
      k <- k + 1
    }

    i <- i + 1
  }

  lambdai <- prob_est * throughput
  
  RO <- numeric()

  k <- 1
  while (k<=num_nodes)
  {
    if (class(x$nodes[[k]]) == "i_MMInf")
      RO[k] <- 0
    else if (class(x$nodes[[k]]) == "i_MM1")
      RO[k] <- lambdai[k] * (1/x$nodes[[k]]$mu)
    else #class i_MMC  
      RO[k] <- lambdai[k] * (1/(x$nodes[[k]]$mu * x$nodes[[k]]$c))
    k <- k + 1
  }

  W <- (x$n / throughput) - x$z
  L <- x$n - (throughput * x$z)

  #if (class(x$prob) != "matrix")
  #  wi <- wi * x$prob

  res <- list(Inputs=x, Lambdai = lambdai, RO=RO, L=L, W=W, Wi = wi, Li = li, Throughput = throughput)
  class(res) <- "o_CJN"
  res
}


Inputs.o_CJN <- function(x, ...) { x$Inputs }
L.o_CJN <- function(x, ...) { x$L }
W.o_CJN <- function(x, ...) { x$W }
RO.o_CJN <- function(x, ...) { x$RO }
Li.o_CJN <- function(x, ...) { x$Li }
Wi.o_CJN <- function(x, ...) { x$Wi }
Throughput.o_CJN <- function(x, ...) { x$Throughput }


summary.o_CJN <- function(object, ...)
{
   
  print("The inputs of the closed Jackson network are:")
  print(object$Inputs)
  print("", quote=FALSE)
  print("The outputs of the closed Jackson network are:")
  print("", quote=FALSE)
  
  print("--------- Per node ---------------------------------")
 
  i <- 1
  while (i <= length(object$RO))
  {
    print(paste("The use of node ", i, " is: ", object$RO[i]))
    print(paste("The mean number of clients in node ", i, " is: ", object$Li[i]))
    print(paste("The response time in node ", i, " is: ", object$Wi[i]))
    print(paste("The throughput of node ", i, " is: ", object$Throughput[i]))
    i <- i + 1
  }

  print("---------- Complete network -------------------------")
  print(paste("The mean number of clients in the network is: ", object$L))
  print(paste("The mean time spend in the network is: ", object$W))  
  print("", quote=FALSE)
}




