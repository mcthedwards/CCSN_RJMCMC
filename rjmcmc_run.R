# Source the functions and libraries
source("rjmcmc_functions.R")
source("rjmcmc_sampler.R")

# Note: We are using a simplified version of the model
#       Doing this in the time domain with standard white noise
#       Standardising signals and PCs
#       Will eventually do this in the frequency domain with Whittle likelihood

# Working directory
path = paste0(getwd(), "/")

# Read in PC design matrix
PC = read.table(paste0(path, "PCs.txt"), header = TRUE)

# Read in test signal
inj = read.table(paste0(path, "test_set.txt"), 
                 header = TRUE)
wave = 1                             # Which test signal to inject?
signal = inj[, wave]                 # Signal
signal = signal / sqrt(var(signal))  # Standardise
n = length(signal)                   # Extract sample size (n = 2^14)
noise = rnorm(n, sd = 1)             # Generate standard white noise
yt = signal + noise                  # Add noise to signal

# Plot data
time = seq(0, 1, length = n)
plot(time, yt, type = "l", col = "grey80",
     xlab = "Time (Seconds)", ylab = "GW Strain")
lines(time, signal, col = "navyblue")

# Run MCMC
mcmc = rjmcmc(yt, PC, 40000, 20000, 10)

# Plot signal reconstruction
res = (n / 2 + 1 - n / 32):(n / 2 + 1 + n / 32 - 1)  # Part to zoom in on
plot(time[res], yt[res], type = "l", col = "grey80",
     xlab = "Time (Seconds)", ylab = "GW Strain")
lines(time[res], signal[res], col = "navyblue")
lines(time[res], mcmc$recon.mean[res], lty = 2, col = 2)
lines(time[res], mcmc$recon.05[res], lty = 3, col = 2)
lines(time[res], mcmc$recon.95[res], lty = 3, col = 2)

# Plot proportion 
index = 1:ncol(PC)
plot.new()
plot.window(xlim = range(index), ylim = c(0, 1))
plot(index, mcmc$PC.prop, pch = 20,
     xlab = "PC Number", ylab = "Pr(Inclusion | Data)",
     col = "navyblue")
segments(index, 0, index, mcmc$PC.prop, col = "navyblue")


