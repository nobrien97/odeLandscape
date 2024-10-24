# Generate a table of molecular trait combinations
setwd("/mnt/c/GitHub/odeLandscape/ODESolver/tests")
# Randomly sample 10000 genotypes
aZ <- exp(rnorm(10000))
bZ <- exp(rnorm(10000))
KZ <- exp(rnorm(10000))
KXZ <- exp(rnorm(10000))
base <- exp(rnorm(10000))
n <- exp(rnorm(10000))
XMult <- exp(rnorm(10000))

samples_NARPAR <- data.frame(
    aZ = aZ,
    bZ = bZ,
    KZ = KZ,
    KXZ = KXZ,
    base = base,
    n = n,
    XMult = XMult
)

write.table(samples_NARPAR, "./samples_PARNAR.csv", sep = ",", row.names = FALSE,
    col.names = FALSE, quote = FALSE)

aY <- exp(rnorm(10000))
bY <- exp(rnorm(10000))
KY <- exp(rnorm(10000))

samples_FFL <- data.frame(
    aY = aY,
    bY = bY,
    KY = KY,
    aZ = aZ,
    bZ = bZ,
    KXZ = KXZ,
    base = base,
    n = n,
    XMult = XMult
)

write.table(samples_FFL, "./samples_FFL.csv", sep = ",", row.names = FALSE,
    col.names = FALSE, quote = FALSE)

aX <- exp(rnorm(10000))
KZX <- exp(rnorm(10000))

samples_FFBH <- data.frame(
    aX = aX,
    KZX = KZX,
    aY = aY,
    bY = bY,
    KY = KY,
    aZ = aZ,
    bZ = bZ,
    KXZ = KXZ,
    base = base,
    n = n,
    XMult = XMult
)

write.table(samples_FFBH, "./samples_FFBH.csv", sep = ",", row.names = FALSE,
    col.names = FALSE, quote = FALSE)

# calculate samples for each motif
system("ODELandscaper -i ./samples_PARNAR.csv -o ./out_NAR.csv -s NAR -t 8 -p 2 -w 0.05")
system("ODELandscaper -i ./samples_PARNAR.csv -o ./out_PAR.csv -s PAR -t 8 -p 2 -w 0.05")
system("ODELandscaper -i ./samples_FFL.csv -o ./out_FFLC1.csv -s FFLC1 -t 8 -p 2 -w 0.05")
system("ODELandscaper -i ./samples_FFL.csv -o ./out_FFLI1.csv -s FFLI1 -t 8 -p 2 -w 0.05")
system("ODELandscaper -i ./samples_FFBH.csv -o ./out_FFLBH.csv -s FFBH -t 8 -p 2 -w 0.05")

library(tidyverse)
library(deSolve)
library(DescTools)

# ODE system for feedback autoregulation:
ODEs_NAR <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- (t > Xstart && t <= Xstop)
    dZ <- base * X + bZ * (t > Xstart && t <= Xstop) * ((X^Hilln)/(KXZ^Hilln + X^Hilln)) * ((KZ^Hilln)/(KZ^Hilln + Z^Hilln)) - aZ*Z
    return(list(c(dZ)))
  })
}

# ODE system for positive feedback autoregulation:
# When Z = 0, this is a stable point
# Will need to give a basal expression level of Z to activate the circuit
ODEs_PAR <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    X <- XMult * (t > Xstart && t <= Xstop)
    
    dZ <- base * X + bZ * (X^Hilln/(KXZ^Hilln + X^Hilln)) * ((Z^Hilln)/((KZ^Hilln)+(Z^Hilln))) - aZ*Z
    return(list(c(dZ)))
  })
}

# ODE system for C1 FFL:
ODEs_C1_FFL <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {

    X <- XMult * (t > Xstart && t <= Xstop)
    dY <- base * X + bY * X^Hilln/(KY^Hilln + X^Hilln) - aY*Y
    dZ <- base * X + bZ * ((X * Y)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dY, dZ)))
  })
}

# I1 FFL
ODEs_I1_FFL <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {

    X <- XMult * (t > Xstart && t <= Xstop)
    
    dY <- base * X + bY * X^Hilln/(KY^Hilln + X^Hilln) - aY*Y
    dZ <- base * X + bZ *  ((X * KY)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dY, dZ)))
  })
}

# ODE system for Feed forward/back hybrid:
ODEs_FFBH <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    
    # X is step function augmented by Z product
    # baseline X given by environmental cue in Xstart -> Xstop
    # change in X at Xstart is 1, change in X at Xstop is -1
    # Z adds a bit more signal to that
    
    # Manually set X
    X <- XMult * (t >= Xstart && t <= Xstop)
    X <- X + XH
    
    # Hill function component of X, XH
    dXH <- ( Z^Hilln / (KZX^Hilln + Z^Hilln) ) - aX*XH
    
    # Update X
    X <- X + dXH
    
    dY <- base * X + bY * X^Hilln/( KY^Hilln + X^Hilln ) - aY*Y
    dZ <- base * X + bZ *  ((X * Y)^Hilln)/((KXZ^Hilln + X^Hilln) * (KY^Hilln + Y^Hilln)) - aZ*Z
    
    return(list(c(dXH, dY, dZ)))
  })
}


# Calculate using deSolve and compare
for (i in 1:nrow(samples_NARPAR)) {
  combos <- samples_NARPAR[i,]
  combos$Xstart <- 1
  combos$Xstop <- 6
      
    #Values for the NAR
    state <- c(Z = 0)
    times <- seq(0, 10, by = 0.1)
    params <- c(Xstart = combos$Xstart, Xstop = combos$Xstop, aZ = combos$aZ, 
                bZ = combos$bZ, KZ = combos$KZ, KXZ = combos$KXZ, base = combos$base,
                Hilln = combos$n, XMult = combos$XMult)
    
    solution <- ode(state, times, ODEs_NAR, params) %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(X = ifelse(time >= params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
      select(time, X, Z)
    
    out <- cbind(i, AUC(solution$time, solution$Z, absolutearea = T))
    colnames(out) = c("modelindex", "Z")
    
    write.table(out, file = "out_AUC_NAR_R.tsv", append = T, sep = "\t", row.names = F,
                col.names = F)


    # Repeat for PAR
    state <- c(Z = 0)
    solution <- ode(state, times, ODEs_PAR, params) %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(X = ifelse(time >= params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
      select(time, X, Z)
    
    out <- cbind(i, AUC(solution$time, solution$Z, absolutearea = T))
    colnames(out) = c("modelindex", "Z")
    
    write.table(out, file = "out_AUC_PAR_R.tsv", append = T, sep = "\t", row.names = F,
                col.names = F)
}

for (i in 1:nrow(samples_FFL)) {
  combos <- samples_FFL[i,]
  combos$Xstart <- 1
  combos$Xstop <- 6
      
    #Values for the C1 FFL
    state <- c(Y = 0, Z = 0)
    times <- seq(0, 10, by = 0.1)
    params <- c(Xstart = combos$Xstart, Xstop = combos$Xstop, aY = combos$aY, 
                bY = combos$bY, KY = combos$KY, KXZ = combos$KXZ, aZ = combos$aZ,
                bZ = combos$bZ, base = combos$base,
                Hilln = combos$n, XMult = combos$XMult)
    
    solution <- ode(state, times, ODEs_C1_FFL, params) %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(X = ifelse(time >= params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
      select(time, X, Y, Z)
    
    out <- cbind(i, AUC(solution$time, solution$Z, absolutearea = T))
    colnames(out) = c("modelindex", "Z")
    
    write.table(out, file = "out_AUC_FFLC1_R.tsv", append = T, sep = "\t", row.names = F,
                col.names = F)


    # Repeat for I1 FFL
    state <- c(Y = 0, Z = 0)    
    solution <- ode(state, times, ODEs_I1_FFL, params) %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(X = ifelse(time >= params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
      select(time, X, Y, Z)
    
    out <- cbind(i, AUC(solution$time, solution$Z, absolutearea = T))
    colnames(out) = c("modelindex", "Z")
    
    write.table(out, file = "out_AUC_FFLI1_R.tsv", append = T, sep = "\t", row.names = F,
                col.names = F)
}

# Now for the FFBH
for (i in 1:nrow(samples_FFBH)) {
  combos <- samples_FFBH[i,]
  combos$Xstart <- 1
  combos$Xstop <- 6
      
    #Values for the C1 FFL
    state <- c(Y = 0, Z = 0)
    times <- seq(0, 10, by = 0.1)
    params <- c(Xstart = combos$Xstart, Xstop = combos$Xstop, aX = combos$aX, KZX = combos$KZX,
                bY = combos$bY, KY = combos$KY, KXZ = combos$KXZ, aZ = combos$aZ,
                bZ = combos$bZ, base = combos$base,
                Hilln = combos$n, XMult = combos$XMult)

    solution <- ode(state, times, ODEs_FFBH, params) %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(X = ifelse(time >= params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
      select(time, X, Y, Z)
    
    out <- cbind(i, AUC(solution$time, solution$Z, absolutearea = T))
    colnames(out) = c("modelindex", "Z")
    
    write.table(out, file = "out_AUC_FFBH_R.tsv", append = T, sep = "\t", row.names = F,
                col.names = F)
}



# Combine data frames
out_NAR_R <- read.table("out_AUC_NAR_R.tsv", sep = "\t", header = F)
out_PAR_R <- read.table("out_AUC_PAR_R.tsv", sep = "\t", header = F)
out_FFLC1_R <- read.table("out_AUC_PAR_R.tsv", sep = "\t", header = F)
out_FFLI1_R <- read.table("out_AUC_PAR_R.tsv", sep = "\t", header = F)

out_R <- read.table("out_AUC_R.tsv", sep = "\t", header = F)

names(out_slim) <- c("modelindex", "A", "B")
names(out_R) <- c("modelindex", "A", "B")

out_slim[names(lhc_samples)[2:4]] <- NA
for (index in unique(out_slim$modelindex)) {
  out_slim[out_slim$modelindex == index, names(lhc_samples)] <- lhc_samples[index, ]
}

out_R[names(lhc_samples)[2:4]] <- NA
for (index in unique(out_R$modelindex)) {
  out_R[out_R$modelindex == index, names(lhc_samples)] <- lhc_samples[index, ]
}

df_combined <- rbind(out_slim, out_R)
df_combined$modelindex[df_combined$modelindex > 512] <- 1:512


# Scale data

df_combined$AScaled <- scale(df_combined$A)
df_combined$BScaled <- scale(df_combined$B)

# Plot distances?
df_combined$distA <- NA
df_combined$distB <- NA

for (idx in unique(df_combined$modelindex)) {
  df_combined[idx,]$distA <- c(dist(df_combined$A[df_combined$modelindex == idx]))
  df_combined[idx,]$distB <- c(dist(df_combined$B[df_combined$modelindex == idx]))
  df_combined[idx,]$scaleddistA <- c(dist(df_combined$AScaled[df_combined$modelindex == idx]))
  df_combined[idx,]$scaleddistB <- c(dist(df_combined$BScaled[df_combined$modelindex == idx]))
}

# Calculate maximum distances

max(df_combined$distA[!is.na(df_combined$distA)])
max(df_combined$distB[!is.na(df_combined$distB)])

write.table(df_combined, file = "out_AUC.tsv", sep = "\t", row.names = F)
