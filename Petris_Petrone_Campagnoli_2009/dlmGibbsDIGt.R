
#################################################################
###                                                           ###
###       Unequal variances - t-distributed innovations       ###
###                                                           ###
#################################################################

dlmGibbsDIGt <- function(y, mod, A_y, B_y, A_theta = A_y, B_theta = A_y,
                         nuRange = c(1 : 10, seq(20, 100, by = 10)),
                         alpha_y = 1 / length(nuRange),
                         alpha_theta,
                         n.sample = 1,
                         thin = 0, ind, save.states = FALSE,
                         progressBar = interactive())
#################################################################
#################################################################
###                                                           ###
### y      : data (univariate time series)                    ###
### mod    : model skeleton                                   ###    
### A_y    : upper bound for a_y (prior)                      ###
### B_y    : upper bound for b_y (prior)                      ###
### A_theta: upper bounds for the components of a_theta,      ###
###            recycled if needed (prior)                     ###
### B_theta: upper bounds for the components of b_theta,      ###
###            recycled if needed (prior)                     ###
### nuRange: set of possible values for the degree-of-freedom ###
###            parameter                                      ###
### alpha_y: parameter of the Dirichlet distribution of pi_y  ###
###            (prior)                                        ### 
### alpha_theta: vector parameters of the Dirichlet           ###
###            distributions of the pi_theta's, stored as     ###
###            columns of a matrix; recycled if vector of     ###
###            length 1 or length(nuRange) (prior)            ###
### n.sample: number of MCMC to output                        ###
### thin  : number of MCMC iteration to skip between two      ###        
###           consecutived stored draws                       ###
### ind   : vector of integers specifying the position of the ###
###           unknown, nonconstant variances. If not          ###
###           specified, all are considered unknown           ###
### save.states: logical; should the simulated states be      ###
###           included in the output?                         ###
### progressBar: logical                                      ###
###                                                           ###
#################################################################
#################################################################
{
    m <- NCOL(y)
    nobs <- NROW(y)
    mod$JW <- matrix(0, nrow = ncol(mod$FF), ncol = ncol(mod$FF))
    r <- ncol(mod$FF)
    if ( hasArg(ind) ) {
        ind <- sort(unique(as.integer(ind)))
        s <- 1 : r
        perm <- s[c(ind, s[ !(s %in% ind)])]
        FF(mod) <- mod$FF[, perm, drop = FALSE]
        GG(mod) <- mod$GG[perm, perm, drop = FALSE]
        mod$W <- mod$W[perm, perm, drop = FALSE]
        p <- length(ind)
    }
    else {
        perm <- 1 : r
        p <- r
    }
    diag(mod$JW)[ 1 : p ] <- 1 : p
    mod$JV <- matrix(p + 1)
    mod$X <- matrix(0, nrow = nobs, ncol = p + 1)
    K <- length(nuRange)
    if ( is.numeric(thin) && (thin <- as.integer(thin)) >= 0 ) 
    {
        every <- thin + 1
        mcmc <- n.sample * every
    }
    else
        stop("\"thin\" must be a nonnegative integer")
    if (!all(c(length(A_theta), length(B_theta)) %in% c(1,p)))
        warning("Unexpected length of \"A_theta\" and/or \"B_theta\"")
    A_theta <- rep(A_theta, length.out = p)
    B_theta <- rep(B_theta, length.out = p)
    ablim_theta <- cbind(A_theta, B_theta)
    ablim_y <- c(A_y, B_y)
    if ( !(length(alpha_y) %in% c(1,K)) )
        warning("Unexpected length of \"alpha_y\"")
    alpha_y <- rep(alpha_y, length.out = K)
    if ( hasArg(alpha_theta) ) {
        if ( is.matrix(alpha_theta) ) {
            if ( nrow(alpha_theta) != K || ncol(alpha_theta) != p )
                stop("Wrong dimension of \"alpha_theta\"")
        } else {
            if ( !(length(alpha_theta) %in% c(1,K)) )
                warning("Unexpected length of \"alpha_theta\"")
            alpha_theta <- matrix(rep(alpha_theta, length.out = K*p), nr = K)
        }
    }
    else 
        alpha_theta <- matrix(rep(alpha_y, length.out = K*p), nr = K)
    theta <- matrix(0, nobs + 1, nrow(mod$W))
    ## initialize
    omega_y <- rep(1, nobs)
    nu_y <- rep(100, nobs)
    lambda_y <- 1
    ab_y <- 0.5 * ablim_y
    pi_y <- alpha_y / sum(alpha_y)
    ##
    omega_theta <- matrix(1, nrow = nobs, ncol = p)
    nu_theta <- matrix(100, nrow = nobs, ncol = p)
    lambda_theta <- rep(1, p)
    ab_theta <- 0.5 * ablim_theta
    pi_theta <- alpha_theta / rep(colSums(alpha_theta), each = K) 

    ## memory allocation
    if ( save.states ) 
        gibbsTheta <- array(0, dim = c(dim(theta), n.sample))
    gibbsOmega_y <- matrix(0, nrow = n.sample, ncol = nobs)
    gibbsNu_y <- matrix(0, nrow = n.sample, ncol = nobs)
    gibbsLambda_y <- numeric(n.sample)
    gibbsAB_y <- matrix(0, nrow = n.sample, ncol = 2)
    gibbsPi_y <- matrix(0, nrow = n.sample, ncol = length(nuRange))
    ##
    gibbsOmega_theta <- array(0, dim = c(nobs, p, n.sample))
    gibbsNu_theta <- array(0, dim = c(nobs, p, n.sample))
    gibbsLambda_theta <- matrix(0, nrow = n.sample, ncol = p)
    gibbsAB_theta <- array(0, dim = c(p, 2, n.sample))
    gibbsPi_theta <- array(0, dim = c(K, p, n.sample))

    ## log target and support for use with `arms'
    ldens.ab <- function(x, lambda, ...) {
        rate <- x[1] / x[2]
        dgamma(lambda, rate * x[1], rate, log = TRUE)
    }
    ind.ab <- function(x, ablim, ...) {
        all( x > 1e-6 & x < ablim )
    }

    ## draw from dirichlet distribution, from package gtools
    rdirichlet <- function (n, alpha) 
    {
        l <- length(alpha)
        x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
        sm <- x %*% rep(1, l)
        x/as.vector(sm)
    }

    it.save <- 0
    if ( progressBar )
        pb <- txtProgressBar(0, mcmc, style = 3)
    for (it in 1 : mcmc)
    {
        if ( progressBar )
            setTxtProgressBar(pb, it)
        mod$X[] <- 1 / cbind( omega_theta * rep(lambda_theta, each = nobs),
                             omega_y * lambda_y)

        ## generate states - FFBS
        modFilt <- dlmFilter(y, mod, simplify=TRUE)
        theta[] <- dlmBSample(modFilt)

        ## generate omega_y
        Sq.y.center <- (y - tcrossprod(theta[-1,],mod$FF))^2
        omega_y <- rgamma(nobs, 0.5 * (1 + nu_y),
                                           0.5 * (lambda_y * Sq.y.center + nu_y))
            
        ## generate nu_y
        for (i in 1:nobs) {
            probs <- dgamma(omega_y[i], 0.5 * nuRange, 0.5 * nuRange) * pi_y
            nu_y[i] <- sample(nuRange, size = 1, prob = probs)
        }

        ## generate pi_y
        nuTable <- table(factor(nu_y, levels=nuRange))
        pi_y <- rdirichlet(1, nuTable + alpha_y)
        
        ## generate lambda_y
        SSy <- crossprod(Sq.y.center, omega_y)
        u <- ab_y[1] / ab_y[2]
        lambda_y <- rgamma(1, u * ab_y[1] + 0.5 * nobs,
                           u + 0.5 * SSy)
        
        ## generate a_y and b_y
        ab_y <- arms(ab_y, myldens = ldens.ab, indFunc = ind.ab, n = 1,
                     lambda = lambda_y, ablim = ablim_y)
        
        ## same story for the 'theta' parameters
        ## omega_theta_i
        Sq.theta.center <- (theta[-1,1:p] - tcrossprod(theta[-(nobs + 1), ],
                                                       mod$GG)[,1:p])^2
        omega_theta[] <- rgamma(nobs * p, 0.5 * (1 + nu_theta),
                                0.5 * (rep(lambda_theta, each=nobs) * Sq.theta.center +
                                       nu_theta))
        ## nu_theta_i
        for (j in 1 : nobs)
            for (i in 1 : p)
            {
                probs <- dgamma(omega_theta[j,i], 0.5 * nuRange, 0.5 * nuRange) *
                    pi_theta[, i]
                nu_theta[j, i] <- sample(nuRange, size = 1, prob = probs)
            }

        ## pi_theta_i
        for (i in 1 : p)
        {
            nuTable <- table(factor(nu_theta[, i], levels = nuRange))
            pi_theta[, i] <- rdirichlet(1, nuTable + alpha_theta[, i])
        }

        ## lambda_theta_i
        SStheta <- colSums(Sq.theta.center * omega_theta)
        u <- ab_theta[, 1] / ab_theta[, 2]
        lambda_theta <- rgamma(p, u * ab_theta + 0.5 * nobs,
                               u + 0.5 * SStheta)

        ## a_theta_i & b_theta_i
        for (i in 1 : p)
            ab_theta[i, ] <- arms(ab_theta[i, ], myldens = ldens.ab,
                                 indFunc = ind.ab, n = 1,
                                 lambda = lambda_theta[i], ablim = ablim_theta[i, ])
        
        ## save
        if ( !(it %% every) )
        {
            it.save <- it.save + 1
            if ( save.states )
                gibbsTheta[,,it.save] <- theta
            gibbsOmega_y[it.save,] <- omega_y
            gibbsNu_y[it.save,] <- nu_y
            gibbsPi_y[it.save,] <- pi_y
            gibbsLambda_y[it.save] <- lambda_y
            gibbsAB_y[it.save,] <- ab_y
            ##
            gibbsOmega_theta[,,it.save] <- omega_theta
            gibbsNu_theta[,,it.save] <- nu_theta
            gibbsPi_theta[,,it.save] <- pi_theta
            gibbsLambda_theta[it.save,] <- lambda_theta
            gibbsAB_theta[,,it.save] <- ab_theta
       }
    }

    if ( progressBar )
        close(pb)
    if ( save.states )
        return(list(omega_y = gibbsOmega_y, nu_y = gibbsNu_y,
                    pi_y = gibbsPi_y, lambda_y = gibbsLambda_y, ab_y = gibbsAB_y,
                    omega_theta = gibbsOmega_theta, nu_theta = gibbsNu_theta,
                    pi_theta = gibbsPi_theta, lambda_theta = gibbsLambda_theta,
                    ab_theta = gibbsAB_theta,
                    theta = gibbsTheta[, order(perm), , drop = FALSE]))
    else
        return(list(omega_y = gibbsOmega_y, nu_y = gibbsNu_y,
                    pi_y = gibbsPi_y, lambda_y = gibbsLambda_y, ab_y = gibbsAB_y,
                    omega_theta = gibbsOmega_theta, nu_theta = gibbsNu_theta,
                    pi_theta = gibbsPi_theta, lambda_theta = gibbsLambda_theta,
                    ab_theta = gibbsAB_theta))

#### end dlmGibbsDIGt ####
}
