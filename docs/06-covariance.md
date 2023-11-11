---
output:
  pdf_document: default
  html_document: default
---


# Covariance Matrix

The current `convey` software does not account for the covariance.  The example below describes this calculation and presents what is, in most cases, a limited impact on the final error terms.  We acknowledge that, under ideal circumstances, the covariance would be properly accounted for in the software.  Modifying our software to account for the covariance remains a future goal.

Influence functions and "resampling" replicates can be used to improve the inference 
about differences and changes between estimates. Put simply, the variance of the difference between 
two estimates $\widehat{T}_1$ and $\widehat{T}_2$ is given by

$$
Var \big( \widehat{T}_1 - \widehat{T}_2 \big) = 
    Var \big( \widehat{T}_1 \big) + Var \big( \widehat{T}_2 \big) - 2 Cov \big( \widehat{T}_1 , \widehat{T}_2 \big)
$$

\noindent where: $Var \big( \widehat{T}_1 \big)$ and $Var \big( \widehat{T}_2 \big)$ 
are the variances of the estimators $\widehat{T}_1$ and $\widehat{T}_2$; 
$Cov \big( \widehat{T}_1 , \widehat{T}_2 \big)$ is the covariance of these estimators.
Usually, the covariance term is ignored. But, specially for changes, the covariance
tends to be positive, which results in a conservative variance estimator for the difference.

Then, how can we estimate the covariance? The `survey` package already provides an 
approach for that using the `svyby`function. Based on the `?survey::svyby` examples, 
we have: 


```r
# load the survey library
library(survey)
```

```
## Carregando pacotes exigidos: grid
```

```
## Carregando pacotes exigidos: Matrix
```

```
## Carregando pacotes exigidos: survival
```

```
## 
## Attaching package: 'survey'
```

```
## The following object is masked from 'package:graphics':
## 
##     dotchart
```

```r
#  load data set
data( api )

# declare sampling design
dclus1 <- svydesign( id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc )

# estimate means
mns <-svyby(~api99, ~stype, dclus1, svymean,covmat=TRUE)

# collect variacnce-covariance matrix of estimates
( m <- vcov( mns ) )
```

```
##          E         H         M
## E 520.5973  573.0404  684.6562
## H 573.0404 1744.2317  747.1989
## M 684.6562  747.1989 1060.1954
```

```r
# compute variance terms
var.naive <- sum( diag( m[c(1,3),c(1,3)] ) )
cov.term <- sum( diag( m[ c(1,3),c(3,1)] ) )

# "naive" SE of the difference
sqrt( var.naive )
```

```
## [1] 39.75918
```

```r
# SE of the difference
sqrt( var.naive - cov.term )
```

```
## [1] 14.54236
```

```r
#... or using svycontrast
svycontrast( mns , c(E = 1, M = -1) )
```

```
##          contrast     SE
## contrast -0.80833 14.542
```

Notice that, because the covariance terms are positive, the (actual) variance of 
the difference is smaller than the "naive" variance.

A similar idea can be implemented with other estimators, such as inequality and poverty measures.
However, this has not yet been implemented for linearization/influence function methods.
In the next section, we show an example with the Gini index.

## Example Calculation

We start by showing what has been implemented --- the resampling based approach:


```r
# load the convey package
library(convey)

# load the survey library
library(survey)

# load the laeken library
library(laeken)

# load the synthetic EU statistics on income & living conditions
data(eusilc)

# make all column names lowercase
names(eusilc) <- tolower(names(eusilc))

# construct a survey.design
# using our recommended setup
des_eusilc <-
  svydesign(
    ids = ~ rb030 ,
    strata = ~ db040 ,
    weights = ~ rb050 ,
    data = eusilc
  )

# create bootstrap-based design object
des_eusilc_rep <- as.svrepdesign( des_eusilc , type = "bootstrap" )

# immediately run the convey_prep function on it
des_eusilc <- convey_prep( des_eusilc )
des_eusilc_rep <- convey_prep(des_eusilc_rep)

# estimate gini by gender, with and without covariance matrices
ginis.rep <- svyby( ~ py010n , ~rb090 , des_eusilc_rep , svygini , na.rm = TRUE , covmat = TRUE )
ginis.rep.nocov <- svyby( ~ py010n , ~rb090 , des_eusilc_rep , svygini , na.rm = TRUE , covmat = FALSE )

# variance of the gini difference: naive
svycontrast( ginis.rep.nocov , quote( `male` - `female` ) )
```

```
## Warning in vcov.svyby(stat): Only diagonal elements of vcov() available
```

```
##             nlcon     SE
## contrast -0.15006 0.0073
```

```r
# variance of the gini difference: accounting for covariance
svycontrast( ginis.rep , quote( `male` - `female` ) )
```

```
##             nlcon     SE
## contrast -0.15006 0.0069
```

While not implemented yet, a linearization-based approach can be applied using:


```r
# estimate gini by gender, with and without covariance matrices
ginis.lin.nocov <- svyby( ~ py010n , ~rb090 , des_eusilc , svygini , na.rm = TRUE , covmat = FALSE )
gini.female <- svygini( ~ py010n , subset( des_eusilc , rb090 == "female" ), na.rm = TRUE , influence = TRUE )
gini.male   <- svygini( ~ py010n , subset( des_eusilc , rb090 == "male" ), na.rm = TRUE , influence = TRUE )

# variance of the gini difference: naive
linmat <- rep( 0 , nrow( des_eusilc$variables ) )
linmat[ attr( gini.female , "index" ) ] <- attr( gini.female , "influence" )
linmat[ attr( gini.male , "index" ) ] <- attr( gini.male , "influence" )
linmat <- linmat * des_eusilc$prob
lintot <- svyby( linmat , ~rb090 , des_eusilc , svytotal , covmat = TRUE )

# collect covariance matrix
( m <- vcov( lintot ) )
```

```
##                male       female
## male   2.954413e-05 2.764584e-10
## female 2.764584e-10 2.292063e-05
```

```r
# SE of the gini difference: naive
SE( svycontrast( ginis.lin.nocov , quote( `male` - `female` ) ) )
```

```
## Warning in vcov.svyby(stat): Only diagonal elements of vcov() available
```

```
##    contrast 
## 0.007243256
```

```r
# SE of the gini difference: accounting for covariance
SE( svycontrast( lintot , quote( `male` - `female` ) ) )
```

```
##    contrast 
## 0.007243218
```
