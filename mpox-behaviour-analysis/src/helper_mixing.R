# authors: Jesse Knight

# constants
tol    = 1e-12
tri.5  = upper.tri(diag(5),TRUE)
diag.5 = c(1,7,13,19,25)
dd.15  = abs(outer(1:15,1:15,`-`))
norm.1 = function(x){ x/sum(x) }
norm.2 = function(x){ x/rowSums(x) }

apply.mix.odds = function(M0,OR){
  M0 = M0 + tol    # avoid NaN issues
  m1 = rowSums(M0) # target row margin
  m2 = colSums(M0) # target col margin
  M = M0 * exp(OR) # apply odds of mixing
  for (i in 1:100){ # iterative proportional fitting: recover target margins
    r1 = m1/rowSums(M); M = sweep(M,1,r1,`*`)
    r2 = m2/colSums(M); M = sweep(M,2,r2,`*`)
    if (all(abs(r1-1) < tol & abs(r2-1) < tol)){ break } } # close enough
  return(M-tol) }
gen.mix.h.odds = function(or){ # odds of mixing by hiv in age:hiv matrix
  OR = matrix(0,10,10)
  OR[1:5,1:5] = OR[6:10,6:10] = or
  return(OR) }
gen.mix.a.odds = function(ors){ # odds of mixing by hiv in age:hiv matrix
  OR = matrix(0,5,5)
  OR[tri.5] = ors                  # copy to upper triangle
  OR = OR + t(OR) - diag(diag(OR)) # copy to lower triangle
  OR = OR = rbind(cbind(OR,OR),cbind(OR,OR))
  return(OR) }
gen.mix.s.odds = function(or,sd=2.5){ # odds of mixing by activity
  OR = or*dnorm(dd.15,sd=sd)*sd
  return(OR) }
gen.mix.s = function(xi,xp,or){ # mixing by activity
  S = outer(xi,xp) / (.5*sum(xi)+.5*sum(xp))
  S = apply.mix.odds(S,gen.mix.s.odds(or)) }
gen.mix.ah = function(x,ors){ # mixing by age:hiv
  AH0 = outer(x,x)/sum(x) # random mixing, population scale
  OR = gen.mix.h.odds(ors[1]) + gen.mix.a.odds(ors[2:16]) # odds of mixing
  AH = apply.mix.odds(AH0,OR) }
aggr.mix.ah.a = function(AH){ # age:hiv -> age (population scale)
  A = AH[1: 5,1: 5] +
      AH[6:10,1: 5] +
      AH[1: 5,6:10] +
      AH[6:10,6:10] }
aggr.mix.ah.h = function(AH){ # age:hiv -> hiv (population scale)
  H = matrix(2,2,data=c(
    sum(AH[1: 5,1: 5]),
    sum(AH[6:10,1: 5]),
    sum(AH[1: 5,6:10]),
    sum(AH[6:10,6:10]))) }
reshape.mix.ah = function(AH){
  AH4 = array(norm.2(AH),c(5,2,5,2)) # ah:hiv (10,10) -> (5,2,5,2) & make probability
}
gen.mix.error = function(ors,x,A,H,fix){
  if (fix){ ors[1] = fix } # only optimize ors[1] for mtl, else fix it
  AH = gen.mix.ah(x,ors)
  err = mean(abs(AH/t(AH)-1)) + # symmetry error
        mean(abs(A-norm.2(aggr.mix.ah.a(AH)))) + # age error
        ifelse(fix,0,mean(abs(H-norm.2(aggr.mix.ah.h(AH))))) # hiv error
  err = min(err,1e3,na.rm=TRUE) } # fix NaN issues
fit.mix.ah = function(x,A,H,fix=0){
  opt = optim(rep(0,16),gen.mix.error,
    x=x,A=A,H=H,fix=fix, # inputs
    method='L-BFGS-B',lower=-6,upper=+6, # max = +/- 6; exp(6) ~= 400x more/less vs random
    control=list(maxit=1e4))
  print(sprintf('%s: err = %.5f @ %d iter',y,opt$value,opt$counts[1]))
  return(opt$par) } # return odds *vector*
