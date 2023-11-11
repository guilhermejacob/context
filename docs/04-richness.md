# Richness



Richness (or Affluence) measures provide another approach for understanding 
income concentration. Unlike inequality measures, that are sensitive to the 
changes over the entire income distribution, richness measures are restricted 
to the top. 

In that sense, they work like "inverted" poverty measures: while the poverty 
focus axiom states that poverty measures should be insensitive to changes 
in the incomes of the non-poor, the richness focus 
axiom states that richness measures should remain unaltered by changes in the 
incomes of the non-rich.

Also like poverty measures, richness measures also rely on a "richness threshold":
those above that threshold are regarded as rich, while those below are regarded 
as non-rich.^[Not necessarily, but including the poor] 
Like for poverty measurement, these richness thresholds can be somewhat arbitrary; 
@medeiros2006 provide a nice review of proposed richness thresholds and also propose
one possible approach with practical and theoretical grounds.

@peichl2010 presented an axiomatic study of richness measures. They classify 
richness measures in two types depending on the effect of income transfers among
the rich. 
A *concave* measure follows the concave transfer axiom (T1): the measure should 
increase when a progressive transfer occurs among two rich individuals.
On the other hand, a *convex* measure follows the convex transfer axiom (T2): 
the measure should *decrease* when a progressive transfer occurs among two rich 
individuals.

The reasoning behind the concave transfer axiom is that such progressive transfer would make the income distribution among the rich more internally homogeneous. By clustering the incomes and assets of the rich away from the rest of society, the distribution of income becomes more polarized. On the other hand, the convex transfer axiom means that a progressive transfer among rich individuals reduces inequality among them, so the convex inequality measure should register a reduction. In a sense, concave measures are related to the income polarization approach, meaning that increases in income polarization should increase the value of the richness measure, while convex measures are related to the income inequality approach, meaning that a reduction in inequality among the rich should also result in a reduction of a convex richness measure.


## Richness Measures (svyrich)


```r
✔️ focuses on the top -- i.e., the "rich"
✔️ can have an inequality or polarization interpretation
✔️ interesting complementary approach to income inequality or polarization
❌ hard to interpret as parameters change
❌ convex richness measures are severely affected by outliers, unreliable for inference
❌ requires a richness line defition
❌ hardly ever used
```

@peichl2010 also presented a general class of richness measures, proposing three
particular sub-classes: the (concave) Chakravarty class, the concave-FGT class and the convex-FGT class, defined at the population level as:

$$
\begin{aligned}
R^{Cha}_\gamma &= \frac{1}{N} \sum_{k \in U} \bigg[( 1 - \bigg( \frac{z_r}{y_k}\bigg)^\gamma \bigg] \delta( y_k \geq z_r ) , \quad \gamma > 0 \\
R^{FGT,T1}_\gamma &= \frac{1}{N} \sum_{k \in U} \bigg( \frac{y_k - z_r}{y_k} \bigg)^\gamma \delta( y_k \geq z_r ) , \quad \gamma \in [0,1) \\
R^{FGT,T2}_\gamma &= \frac{1}{N} \sum_{k \in U} \bigg( \frac{y_k - z_r}{z_r} \bigg)^\gamma \delta( y_k \geq z_r ) , \quad \gamma > 1 \\
\end{aligned}
$$
\noindent where $z_r$ is the richness threshold and $\gamma$ is a sensitivity parameter. 

To estimate these measures, @brz2014 proposed the estimators

$$
\begin{aligned}
\widehat{R}^{Cha}_\gamma &= \frac{1}{\widehat{N}} \sum_{k \in s} w_k \bigg[( 1 - \bigg( \frac{z_r}{y_k}\bigg)^\gamma \bigg] \delta( y_k \geq z_r ) , \quad \gamma > 0 \\
\widehat{R}^{FGT,T1}_\gamma &= \frac{1}{\widehat{N}} \sum_{k \in s} w_k \bigg( \frac{y_k - z_r}{y_k} \bigg)^\gamma \delta( y_k \geq z_r ) , \quad \gamma \in [0,1) \\
\widehat{R}^{FGT,T2}_\gamma &= \frac{1}{\widehat{N}} \sum_{k \in s} w_k \bigg( \frac{y_k - z_r}{z_r} \bigg)^\gamma \delta( y_k \geq z_r ) , \quad \gamma > 1 \\
\end{aligned}
$$
\noindent where $w_k$ is the sampling (or calibration) weight.

In order to estimate the variance of these estimators, @brz2014 derived influence functions under the @deville1999 approach. These functions are the ones used in the `svyrich` function.

@brz2014 also studied the reliability of the inference based on these estimators 
using a model-based Monte Carlo simulation, which are also valid for (design-based) 
simple random sampling with replacement. His results showed that inference for convex 
measures are unreliable for convex measures, the estimator generally produces bad results.  It is highly sensitive to outliers, and confidence intervals are invalid (actual coverage is much smaller than 95%).  The vignette below shows for a design-based approach with the same results.

---

### Monte Carlo Simulation

@brz2014 presented results using a model-based Monte Carlo --- i.e., he simulated 
several samples from the model and compared the behavior of the estimators with 
the model-based parameters. 
In the simulation below, we take a design-based approach: we take several samples 
from a fixed finite population using a particular sampling design, compute the estimator 
on each and compare them to the finite population parameter (not the model-based 
parameter!).

For the sake of similarity, we start by simulating a large finite population ($N = 10^5$) 
using the same distribution from @brz2014:


```r
# load libraries
library(survey)
library(convey)
library(sampling)
library(VGAM)

# set random seed
set.seed(2023)

# superpopulation parameters
scale.x0 <- 1
shape.theta <- 2
cha.beta <- 2
fgt.alpha <- 1.5
n.pop <- as.integer(10 ^ 5)

# generate finite population
pop.df <-
  data.frame(y1 = rparetoI(n.pop  ,  scale.x0 , shape.theta))
```

Then, we compute the finite population parameters using the simulated population:


```r
# richness measures: finite population parameters
cha.scores <-
  function(y , b , rho)
    ifelse(y > rho , (1 - (rho / y) ^ b) , 0)
fgtt2.scores <-
  function(y , g , rho)
    ifelse(y > rho , (y / rho - 1) ^ a , 0)
median.fp <- quantile(pop.df$y1 , .50)
rho.fp <- 3 * median.fp
rHC.fp <- mean(pop.df$y1 > rho.fp)
rCha.fp <- mean(cha.scores(pop.df$y1  , cha.beta , rho.fp))
rFGTT2.fp <- mean(cha.scores(pop.df$y1  , fgt.alpha , rho.fp))
```

For our sampling design, we select $n = 1000$ units using multinomial sampling, 
with the variable `x.aux` as the size variable for the inclusion probabilities:


```r
# define sample size
n.sample <- 1000L

# inclusion probability
pop.df$x.aux <- plogis( pop.df$y1 ) / 1.1
pop.df$pi1 <- sampling::inclusionprobabilities( pop.df$x.aux , n.sample )
```

We run the procedure N and store the estimate objects using the code below:


```r
# define the number of simulation runs
mc.reps <- 5000L

# simulation runs
rep.list <- lapply(seq_len(mc.reps) , function(this.iter) {
  # multinomial sampling
  this.sample <- sampling::UPmultinomial(pop.df$pi1)
  this.sample <- rep(1:n.pop , this.sample)
  sample.df <- pop.df[this.sample ,]
  sample.df$weights <- 1 / sample.df$pi1
  des.obj <-
    svydesign(
      ids = ~ 1 ,
      weights = ~ weights ,
      data = sample.df ,
      nest = FALSE
    )
  
  # run estimation
  des.obj <- convey_prep(des.obj)
  rCha.hat <-
    svyrich(
      ~ y1 ,
      des.obj ,
      type_measure = "Cha" ,
      g = cha.beta ,
      type_thresh = "relq" ,
      percent = 3
    )
  suppressWarnings(
    rHC.hat <-
      svyrich(
        ~ y1 ,
        des.obj ,
        type_measure = "FGTT1" ,
        g = 0 ,
        type_thresh = "relq"  ,
        percent = 3
      )
  )
  suppressWarnings(
    rFGTT2.hat <-
      svyrich(
        ~ y1 ,
        des.obj ,
        type_measure = "FGTT2" ,
        g = fgt.alpha ,
        type_thresh = "relq"  ,
        percent = 3
      )
  )
  est.list <- list(rHC.hat , rCha.hat , rFGTT2.hat)
  est.list
  
})
```

To study the behavior of the estimators, we estimate their expected values, empirical 
variance (for the main parameter) and mean squared error. To study the validity 
of the normal approximation, we also estimate the percent coverage rate of the nominal 
95% confidence interval. This is done using the function below:


```r
sim.compile <- function(ll ,
                        pv ,
                        level = .95 ,
                        na.rm = FALSE) {
  # collect estimates
  mhat.vec <- sapply(ll , coef)
  vhat.vec <- sapply(ll , vcov)
  
  # estimate expected value
  mhat.exp <- mean(mhat.vec , na.rm = na.rm)
  vhat.exp <- mean(vhat.vec , na.rm = na.rm)
  
  # calculate empirical variance
  mhat.empvar <- var(mhat.vec , na.rm = na.rm)
  
  # estimate squared bias
  mhat.bias2 <- (mhat.exp - pv) ^ 2
  vhat.bias2 <- (vhat.exp - mhat.empvar) ^ 2
  
  # estimate mse
  mhat.mse <- mhat.bias2 + mhat.empvar
  
  # estimate coverage rate
  ci.hats <- t(sapply(ll , confint))
  ci.check <-
    matrix(as.logical(NA) , nrow = nrow(ci.hats) , ncol = 3)
  ci.check[, 1] <- ci.hats[, 1] <= pv
  ci.check[, 2] <- ci.hats[, 2] >= pv
  ci.check[, 3] <- apply(ci.check[, 1:2] , 1 , all)
  pcr.emp <- mean(ci.check[, 3] , na.rm = na.rm)
  
  # setup final table
  data.frame(
    "mc.reps" = length(ll) ,
    "theta" = pv ,
    "theta.hat" = mhat.exp ,
    "theta.bias2" = mhat.bias2 ,
    "theta.empvar" = mhat.empvar ,
    "theta.hat.mse" = mhat.mse ,
    "theta.varhat" = vhat.exp ,
    "pcr" = pcr.emp
  )
}
```

For the Headcount Richness Ratio (computed using the concave FGT measure), we have:


```r
rhc.list <- lapply(rep.list , `[[` , 1)
sim.compile(rhc.list, rHC.fp)
```

```
##   mc.reps   theta  theta.hat  theta.bias2 theta.empvar theta.hat.mse
## 1    5000 0.05469 0.05452618 2.683797e-08 4.112616e-05    4.1153e-05
##   theta.varhat    pcr
## 1 4.336252e-05 0.9524
```

```r
stopifnot(round(
  sim.compile(rhc.list, rHC.fp)["theta.bias2"] / sim.compile(rhc.list, rHC.fp)["theta.hat.mse"] ,
  4
) == 0.0007)

stopifnot(round(
  sim.compile(rhc.list, rHC.fp)["theta.varhat"] / sim.compile(rhc.list, rHC.fp)["theta.empvar"] ,
  2
) == 1.05)

stopifnot(round(sim.compile(rhc.list, rHC.fp)["pcr"] , 4) == 0.9524)
```

Under this approach, the squared bias accounts for approx. 0.07% of the MSE, indicating 
that the MSE of the estimator is reasonably approximated by its variance. Additionally, 
the ratio between the expected value of the variance estimator and the empirical 
variance is approx. 1.05, indicating that the variance estimates are expected to 
be a (slighlty conservative, but) good approximation of the empirical variance. 
Finally, the estimated percent coverage rate of 95.24% is close to the nominal 
level of 95%, indicating that the confidence intervals are approximately valid.

For the Chakravarty measure, we have:


```r
rcha.list <- lapply(rep.list , `[[` , 2)
sim.compile(rcha.list, rCha.fp)
```

```
##   mc.reps      theta  theta.hat  theta.bias2 theta.empvar theta.hat.mse
## 1    5000 0.02743915 0.02737177 4.538966e-09 1.428639e-05  1.429093e-05
##   theta.varhat    pcr
## 1 1.402416e-05 0.9428
```

```r
stopifnot(round(
  sim.compile(rcha.list, rCha.fp)["theta.bias2"] / sim.compile(rcha.list, rCha.fp)["theta.hat.mse"] ,
  4
) == 0.0003)

stopifnot(round(
  sim.compile(rcha.list, rCha.fp)["theta.varhat"] / sim.compile(rcha.list, rCha.fp)["theta.empvar"] ,
  2
) == 0.98)

stopifnot(round(sim.compile(rcha.list, rCha.fp)["pcr"] , 4) == 0.9428)
```

Under this approach, the squared bias is approx 0.03% of the MSE, indicating that
the MSE of the estimator is reasonably approximated by its variance. Additionally, 
the ratio between the expected value of the variance estimator and the empirical 
variance is approx. 0.98, indicating that the variance estimates are expected to 
be a good approximation of the empirical variance. Finally, the estimated percent 
coverage rate of 94.28% is close to the nominal level of 95%, indicating that the 
confidence intervals are approximately valid.

For the convex FGT measure, we have:


```r
rfgtt2.list <- lapply(rep.list , `[[` , 3)
sim.compile(rfgtt2.list, rFGTT2.fp)
```

```
##   mc.reps      theta theta.hat theta.bias2 theta.empvar theta.hat.mse
## 1    5000 0.02349189 0.1045973  0.00657808   0.01865219    0.02523027
##   theta.varhat    pcr
## 1   0.01857144 0.5212
```

```r
stopifnot(round(
  sim.compile(rfgtt2.list, rFGTT2.fp)["theta.bias2"] / sim.compile(rfgtt2.list, rFGTT2.fp)["theta.hat.mse"] ,
  2
) == 0.26)

stopifnot(round(
  sim.compile(rfgtt2.list, rFGTT2.fp)["theta.varhat"] / sim.compile(rfgtt2.list, rFGTT2.fp)["theta.empvar"] ,
  2
) == 1)

stopifnot(round(sim.compile(rfgtt2.list, rFGTT2.fp)["pcr"] , 4) == 0.5212)
```

Under this approach, the squared bias is approx 26% of the MSE, indicating that 
the *bias is substantial*. The ratio between the expected value of the variance 
estimator and the empirical variance is approx. 1.00, indicating that the variance 
estimates are expected to be a good approximation of the empirical variance (but 
not of the MSE!). Finally, the estimated percent coverage rate of 52.12% is far 
from the nominal level of 95%, indicating that the confidence intervals are invalid. 
This comes from the fact that the estimator is very sensitive to the extreme values 
and is the reason why @brz2014 does not recommend using convex richness measures.
