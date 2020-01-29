## aux functions

# zeros |----------------------------------------------------------------------
ones <- function(n = NULL, k = NULL) {
    if (!is.null(n) & !is.null(k)) {
        return(matrix(1, nrow = n, ncol = k))
    } else if(any(!is.null(n), !is.null(k))) {
        return(rep(1, max(n, k)))
    }
}

zeros <- function(n, k) {
    if (!is.null(n) & !is.null(k)) {
        return(matrix(0, nrow = n, ncol = k))
    } else if(any(!is.null(n), !is.null(k))) {
        return(rep(0, max(n, k)))
    }
}

eye <- function(d) {
    return(diag(d))
}


# mlag2 |----------------------------------------------------------------------
mlag2 <- function(X, p) {
    size = dim(X)
    Traw = size[1]
    N = size[2]
    Xlag = zeros(n = Traw, k = N * p)
    for (ii in 1:p) {
        Xlag[((p + 1):Traw), (N * (ii - 1) + 1):(N * ii)] <-
            X[(p + 1 - ii):(Traw - ii), 1:N]
    }
    return(Xlag)
}

# create_RHS |-----------------------------------------------------------------
create_RHS <- function(YY, M, p, t) {
    # K is the number of elements in the state vector
    K <-  M + p * (M ^ 2)
    # Create x_t matrix.
    # first find the zeros in matrix x_t
    x_t <-  zeros((t - p) * M, K)
    for (i in 1:(t - p)) {
        ztemp = eye(M)
        for (j in 1:p) {
            xtemp = YY[i, ((j - 1) * M + 1):(j * M)]
            xtemp = kronecker(eye(M), xtemp)
            ztemp = cbind(ztemp, xtemp)
        }
        x_t[((i - 1) * M + 1):(i * M),] = ztemp
    }
    return(list(x_t, K))
}
