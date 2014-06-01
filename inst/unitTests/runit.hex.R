### --- Test setup ---
library("gtm")

## tolerance
mytol = .Machine$double.eps^0.5

## hexagon function; find a use for this in gtm.stp2
gtm.hex.old <-function (xDim, yDim) {
    if ((xDim < 2) || (yDim < 2) || (yDim != floor(yDim)) || 
        (xDim != floor(xDim))) {
      stop("Invalid grid dimensions")
    }
    r1 = 0:(xDim - 1)
    r2 = (yDim - 1):0
    fx = function(x, y) return(x)
    fry = function(x, y) return(y * sqrt(3)/2)
    X = matrix(outer(r1, r2, fx), ncol = yDim)
    Y = outer(r1, r2, fry)
    i = 2
    while (i <= yDim) {
        X[, i] = X[, i] + 0.5
        i = i + 2
    }
    grid = cbind(matrix(X, ncol = 1), matrix(Y, ncol = 1))
    maxVal = max(grid)
    grid = grid * (2/maxVal)
    maxXY = apply(grid, 2, max)
    grid[, 1] = grid[, 1] - maxXY[1]/2
    grid[, 2] = grid[, 2] - maxXY[2]/2
    grid
}

old = gtm.hex.old(10,10)
new = gtm.hex(10,10)

## check for regressions if you neaten up the function
test.hex.diff <- function() {  
  checkEqualsNumeric(old, new, tolerance= mytol)
}
