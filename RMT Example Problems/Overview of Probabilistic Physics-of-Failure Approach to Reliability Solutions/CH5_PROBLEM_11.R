# CHAPTER 5 PROBLEM 11
# Reuel Smith
# =================================================
# Compute standard constant for computing the time to reach pit transition point
C.STANDARD <- ((2*96514*7.87*(1/10000^2)*pi)/(55.85*0.8))*(250^2)
35000/8.314

Mean.Hoop.Stress <- (7.5*(0.01+0.5*0.25))/0.01

# Crack Propagtion Time
t_2 <- (24/15)*(1/(1e-12*(1.12^3)*(202.5^3)*(pi^1.5)*(-0.5)))*((0.01^-0.5) - (0.00025^-0.5))

(55.85*0.8)/(2*96514*7.87*(1/10000^2)*pi)

(250^2)
# Now Run an MCMC to find the first crack to breach the pipe
T <- 500              # Set Temperature to something relatively high (K)
# Set the pit growth standard deviation to 0.52 and randomly draw actual size of pit to recalculate growth time
SD.Pitsize <- 1/sqrt((1/0.8^2) + (1/0.6905921^2))

# Set surface area of the Pipe in cm^2
A.PIPE <- pi*50*8e5

# RUN MCMC FOR MINIMUM TIME TO FAILURE
N.SAMPLE <- 10000
for(i in 1:N.SAMPLE){
  N.PITS <- ceiling(A.PIPE*rnorm(1,202.79, 28.78))     # RANDOM DRAW PIT NUMBER
  r.NEW <- rlnorm(N.PITS,log(250),SD.Pitsize)          # Sample true pit radius of whole set i

  # Split radii that are beyond the threshold and radii that are before
  r.NEW.GREATER <- r.NEW[which(r.NEW > 250)]
  r.NEW.LESS <- r.NEW[which(r.NEW <= 250)]

}
