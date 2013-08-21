require("adaptive.sampling.gp")


#x  <- as.matrix(expand.grid(x = 1:20/4, y = 1:20/4))

start <- Sys.time ()
nx = 30	# dimensionality of x
nxp = 10
x = as.matrix(expand.grid(x = seq(1, 5, length=nx), y = seq(1, 5, length=nx)))
gp <- new.gp(x, 0.5, kernel.exponential(0.9, 1))

xp <- matrix(0,2, nrow=2)
xp = as.matrix(expand.grid(seq(1, 5, length=nxp), seq(1, 5, length=nxp)))

#xp[1,] <- c(2,2)
#xp[2,] <- c(4,4)
#xp[3,] <- c(1,2)
#xp[4,] <- c(3,3)

#yp <- c(0.2, 0.8)
yp <- rep(0.2, nrow(xp))
yp[floor(nrow(xp)/2+1):nrow(xp)] = 0.8
yp[88:90]=0.1
#yp[3] = 0.8
# measurement noise
ep <- rep(0.1, nrow(xp))

gp <- posterior(gp, xp, yp, ep)

print(Sys.time () - start)

plot(gp)

ex = new.experiment(kernelf=kernel.exponential(0.9, 1), type="gaussian")
for(i in 1:length(yp))
{
	add.measurement(ex, yp[i], xp[i,])
}

gp1 = posterior(ex, x)
