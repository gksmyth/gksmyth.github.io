neldmead <- function(fobj, startvec, print.level=0, tol=1e-6, ...)
#  Function minimization by the Nelder-Mead simplex method. This is an
#  inefficient minimizer but is useful for obtaining a decent starting
#  point and scale information for one of the efficient methods.
#  Ref: Press (Numerical Methods in C) p. 305.

#  Posted by Bill Clark to the S-News news group in March 1997.

#  Minor modifications by
#  Gordon Smyth, Walter and Eliza Hall Institute, smyth@wehi.edu.au
#  2 Oct 1999. Last modified 24 April 2001.

#  Corrections to the algorithm by 
#  David Clifford, University of Chicago, clifford@galton.uchicago.edu
#  9 Aug 2002.
{
# Standard tuning parameters.
alpha <- 1.0 # alpha >= 1
beta  <- 0.5 # 0 < beta < 1
gamma <- 2.0 # gamma > alpha
# Build initial simplex by scaling single elements of starting vector by
# a factor exp(+/-lambda).
if (print.level==2) cat('Initializing simplex...\n')
lambda <- 2
npar <- length(startvec)
nverts <- npar + 1
verts <- matrix(startvec,nverts,npar,byrow=T)
fvals <- rep(0,nverts)
for (iv in 1:npar)
  repeat
    {
      lamb <- lambda
      {
        verts[iv,iv] <- exp(lamb) * startvec[iv]
        if (verts[iv,iv]==0) verts[iv,iv] <- 1
        fvals[iv] <- fobj(verts[iv,],...)
        if (is.finite(fvals[iv])) break
        verts[iv,iv] <- exp(-lamb) * verts[iv,iv]
        fvals[iv] <- fobj(verts[iv],...)
        if (is.finite(fvals[iv])) break
        lamb <- lamb/2
      }
    }
verts[nverts,] <- startvec
fvals[nverts] <- fobj(startvec,...)
scatter.init <- apply(verts,2,stdev)
centroid <- startvec
# Start iteration.
iter <- 0
repeat
  {
    # Check for convergence.
    if (print.level==2) {
      cat("Press <return> to continue, q<return> to quit\n")
      out <- readline()
      if(out=="q") break
    }
    iter <- iter + 1
    scatter <- apply(verts,2,stdev)
    conv.test <- max(scatter/scatter.init)
    fmax <- max(fvals)
    iv.fmax <- min((1:nverts)[fvals==fmax])
    fnext <- max(fvals[(1:nverts) != iv.fmax])
    iv.fnext <- min((1:nverts)[fvals==fnext])
    fmin <- min(fvals)
    iv.fmin <- min((1:nverts)[fvals==fmin])
    if (conv.test < tol) break
    if (print.level > 0)
      {
        fcen <- fobj(centroid,...)
        cat('Iteration',iter,'... convergence test =',signif(conv.test,3),
            'fval =',signif(fcen,8),'\n')
      }
    # Make a move. First try projecting across centroid from high point.
    # The centroid referred to here is the centroid of the face
    # across from the point where the function is largest.
    centroid <- if (npar==1) verts[iv.fmin] else
    apply(verts[(1:nverts) != iv.fmax,],2,mean)
    proj.vert <- centroid + alpha * (centroid - verts[iv.fmax,])
    proj.fval <- fobj(proj.vert,...)
    if (print.level == 2)
      {
        cat('Vertices:\n')
        for (iv in 1:nverts) cat(iv,signif(verts[iv,]),'\n')
        cat('Function values:\n',signif(fvals),'\n')
        cat('Projected vertex:\n',signif(proj.vert),'\n')
        cat('Projected function value:',signif(proj.fval),'\n')
      }
    # Make sure function values are different.
    if (all(fvals==fmax))
      {
        cat(paste('nldrmd error: Function values same at all vertices.',
                  'Change startvec.\n'))
        return(NA)
      }
        # Decide what to do based on projected function value.
    if (proj.fval >= fmin & proj.fval <= fnext)
      # Middling result. Take it.
      {
        if (print.level==2) cat('Taking projection.\n')
        verts[iv.fmax,] <- proj.vert
        fvals[iv.fmax] <- proj.fval
        next
      }
    if (proj.fval < fmin)
      # Better than the minimum. Perform an expansion.
      {
        if (print.level==2) cat('Performing expansion.\n')
        #### HERE IS THE PROBLEM>>> its fixed, ... use proj.vert instead
        exp.vert <- centroid - gamma * (centroid - proj.vert)
        exp.fval <- fobj(exp.vert,...)
        if(print.level==2) {
          cat('Expansion vertex:\n',signif(exp.vert),'\n')
          cat('Expansion function value:',signif(exp.fval),'\n')
        }
        if(exp.fval < fmin) {
          verts[iv.fmax,] <- exp.vert
          fvals[iv.fmax] <- exp.fval
        } else {
          verts[iv.fmax,] <- proj.vert
          fvals[iv.fmax] <- proj.fval
        }
        next
      }
    if (proj.fval > fnext)
      # Worse than second highest. Try a one-dimensional contraction.
      {
        if(proj.fval < fvals[iv.fmax]) {
          verts[iv.fmax,] <- proj.vert
          fvals[iv.fmax] <- proj.fval
        }
        cont.vert <- (1-beta) * centroid + (beta) * verts[iv.fmax,]
        cont.fval <- fobj(cont.vert,...)
        if (cont.fval < fmax)
          {
            if (print.level==2) cat('Performing 1D contraction.\n')
            if(print.level==2) {
              cat('Contraction vertex:\n',signif(cont.vert),'\n')
              cat('Contraction function value:',signif(cont.fval),'\n')
            }
            verts[iv.fmax,] <- cont.vert
            fvals[iv.fmax] <- cont.fval
            next
          }
        # Didn't work. Shrink toward best point.
        if (print.level==2) cat('Performing wholesale contraction.\n')
        for (iv in 1:nverts)
          {
            verts[iv,] <- beta * verts[iv.fmin,] + (1 - beta) * verts[iv,]
            fvals[iv] <- fobj(verts[iv,],...)
          }
      }
  # Done with update.
  }
# Convergence test succeeded.
return(x=verts[iv.fmin,],fmin=min(fvals),iter=iter-1)
}
