# Generate a table of molecular trait combinations
aZ <- seq(0.0, 30.0, by = 1)
bZ <- seq(0.0, 30.0, by = 1)
KZ <- seq(0.0, 30.0, by = 1)
KXZ <- seq(0.0, 30.0, by = 1)

samples <- expand.grid(aZ, bZ, KZ, KXZ)

write.table(samples, "./samples.csv", sep = ",", row.names = FALSE,
    col.names = FALSE, quote = FALSE)
