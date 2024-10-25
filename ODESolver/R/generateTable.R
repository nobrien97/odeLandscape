library(tidyverse)
library(deSolve)
library(DescTools)
library(paletteer)
library(ggh4x)
library(ggpmisc)

# Generate a table of molecular trait combinations
setwd("/mnt/c/GitHub/odeLandscape/ODESolver/tests")

NUM_SAMPLES <- 1000

# Randomly sample NUM_SAMPLES genotypes
aZ <- exp(rnorm(NUM_SAMPLES))
bZ <- exp(rnorm(NUM_SAMPLES))
KZ <- exp(rnorm(NUM_SAMPLES))
KXZ <- exp(rnorm(NUM_SAMPLES))
base <- exp(rnorm(NUM_SAMPLES))
n <- exp(rnorm(NUM_SAMPLES))
XMult <- exp(rnorm(NUM_SAMPLES))

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

aY <- exp(rnorm(NUM_SAMPLES))
bY <- exp(rnorm(NUM_SAMPLES))
KY <- exp(rnorm(NUM_SAMPLES))

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

aX <- exp(rnorm(NUM_SAMPLES))
KZX <- exp(rnorm(NUM_SAMPLES))

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

# ODE system for feedback autoregulation:
ODEs_NAR <- function(t, state, parameters) {
  with (as.list(c(state, parameters)), {
    # step function leads to numerical issues in lsoda:
    #dZ <- bZ * (t > Xstart && t <= Xstop & Z<1) - aZ*Z
    # use Hill function instead:
    X <- XMult * (t > Xstart && t <= Xstop)
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
    state <- c(XH = 0, Y = 0, Z = 0)
    times <- seq(0, 10, by = 0.1)
    params <- c(Xstart = combos$Xstart, Xstop = combos$Xstop, aX = combos$aX, KZX = combos$KZX,
                aY = combos$aY, bY = combos$bY, KY = combos$KY, KXZ = combos$KXZ, aZ = combos$aZ,
                bZ = combos$bZ, base = combos$base,
                Hilln = combos$n, XMult = combos$XMult)

    solution <- ode(state, times, ODEs_FFBH, params) %>%
      as.data.frame() %>%
      as_tibble() %>%
      mutate(X = ifelse(time >= params["Xstart"] & time <= params["Xstop"], 1, 0)) %>%
      select(time, X, XH, Y, Z)
    
    out <- cbind(i, AUC(solution$time, solution$Z, absolutearea = T))
    colnames(out) = c("modelindex", "Z")
    
    write.table(out, file = "out_AUC_FFBH_R.tsv", append = T, sep = "\t", row.names = F,
                col.names = F)
}


# Combine data frames

# Load in R
out_NAR_R <- read.table("out_AUC_NAR_R.tsv", sep = "\t", header = F)
out_PAR_R <- read.table("out_AUC_PAR_R.tsv", sep = "\t", header = F)
out_FFLC1_R <- read.table("out_AUC_FFLC1_R.tsv", sep = "\t", header = F)
out_FFLI1_R <- read.table("out_AUC_FFLI1_R.tsv", sep = "\t", header = F)
out_FFBH_R <- read.table("out_AUC_FFBH_R.tsv", sep = "\t", header = F)

names(out_NAR_R) <- c("modelindex", "Z")
names(out_PAR_R) <- c("modelindex", "Z")
names(out_FFLC1_R) <- c("modelindex", "Z")
names(out_FFLI1_R) <- c("modelindex", "Z")
names(out_FFBH_R) <- c("modelindex", "Z")

out_NAR_R$motif <- "NAR"
out_PAR_R$motif <- "PAR"
out_FFLC1_R$motif <- "FFLC1"
out_FFLI1_R$motif <- "FFLI1"
out_FFBH_R$motif <- "FFBH"

out_NAR_R$calcMode <- "R"
out_PAR_R$calcMode <- "R"
out_FFLC1_R$calcMode <- "R"
out_FFLI1_R$calcMode <- "R"
out_FFBH_R$calcMode <- "R"

# Adjust modelindex so they are unique
out_PAR_R$modelindex = out_PAR_R$modelindex + 100
out_FFLC1_R$modelindex = out_FFLC1_R$modelindex + 200
out_FFLI1_R$modelindex = out_FFLI1_R$modelindex + 300
out_FFBH_R$modelindex = out_FFBH_R$modelindex + 400


out_NAR_C <- read.table("out_NAR.csv", sep = ",", header = F)
out_PAR_C <- read.table("out_PAR.csv", sep = ",", header = F)
out_FFLC1_C <- read.table("out_FFLC1.csv", sep = ",", header = F)
out_FFLI1_C <- read.table("out_FFLI1.csv", sep = ",", header = F)
out_FFBH_C <- read.table("out_FFLBH.csv", sep = ",", header = F)

out_NAR_C <- as.data.frame(out_NAR_C[,2])
out_PAR_C <- as.data.frame(out_PAR_C[,2])
out_FFLC1_C <- as.data.frame(out_FFLC1_C[,2])
out_FFLI1_C <- as.data.frame(out_FFLI1_C[,2])
out_FFBH_C <- as.data.frame(out_FFBH_C[,2])

names(out_NAR_C) <- c("Z")
names(out_PAR_C) <- c("Z")
names(out_FFLC1_C) <- c("Z")
names(out_FFLI1_C) <- c("Z")
names(out_FFBH_C) <- c("Z")

out_NAR_C$modelindex <- 1:nrow(out_NAR_C)
out_PAR_C$modelindex <- 1:nrow(out_PAR_C) + 100
out_FFLC1_C$modelindex <- 1:nrow(out_FFLC1_C) + 200
out_FFLI1_C$modelindex <- 1:nrow(out_FFLI1_C) + 300
out_FFBH_C$modelindex <- 1:nrow(out_FFBH_C) + 400

out_NAR_C$motif <- "NAR"
out_PAR_C$motif <- "PAR"
out_FFLC1_C$motif <- "FFLC1"
out_FFLI1_C$motif <- "FFLI1"
out_FFBH_C$motif <- "FFBH"

out_NAR_C$calcMode <- "C"
out_PAR_C$calcMode <- "C"
out_FFLC1_C$calcMode <- "C"
out_FFLI1_C$calcMode <- "C"
out_FFBH_C$calcMode <- "C"

# Combine data frames
df_combined <- rbind(out_NAR_C, out_PAR_C, out_FFLC1_C, out_FFLI1_C, out_FFBH_C,
                     out_NAR_R, out_PAR_R, out_FFLC1_R, out_FFLI1_R, out_FFBH_R)
# Scale data
df_combined$ZScaled <- scale(df_combined$Z)

df_combined <- pivot_wider(df_combined, names_from = calcMode, values_from = Z)

df_combined <- df_combined %>%
  group_by(modelindex, motif) %>%
  fill(R, C) %>%
  fill(R, C, .direction = "up") %>%
  distinct()

ggplot(df_combined,
  aes(x = C, y = R, colour = motif)) +
  facet_nested(.~"Motif" + motif) +
  geom_point() +
  stat_poly_line(method = "lm") +
  stat_poly_eq(label.y = 0.95, rr.digits = 5) +          
  scale_colour_paletteer_d("MetBrewer::Johnson") +
  labs(x = "Z AUC Ascent", y = "Z AUC deSolve", colour = "Motif") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "bottom")
    
ggsave("AUC_comparison.png", device = png, width = 10, height = 3)
