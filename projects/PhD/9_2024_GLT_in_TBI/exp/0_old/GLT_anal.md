GLT distribution analysis
=========================

## References
- https://stats.stackexchange.com/questions/435644/is-there-a-method-to-look-for-significant-difference-between-two-linear-regressi


## LM
### Ctrl
Call:
lm(formula = int ~ cluster, data = filter(df, group == "cont"))

Residuals:
      Min        1Q    Median        3Q       Max 
-13634446  -6812003  -3690274   2935023  47113880 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept) -40181010   14864044  -2.703   0.0146 *  
cluster        125825      23554   5.342 4.46e-05 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 13610000 on 18 degrees of freedom
Multiple R-squared:  0.6132,	Adjusted R-squared:  0.5917 
F-statistic: 28.54 on 1 and 18 DF,  p-value: 4.456e-05

### Cef
Cef group
Call:
lm(formula = int ~ cluster, data = filter(df, group == "cef"))

Residuals:
     Min       1Q   Median       3Q      Max 
-8449937 -4134655 -2104998  1313614 20958263 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1912566    2358078   0.811    0.423    
cluster        15528       2469   6.290  4.7e-07 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 6539000 on 32 degrees of freedom
Multiple R-squared:  0.5529,	Adjusted R-squared:  0.5389 
F-statistic: 39.57 on 1 and 32 DF,  p-value: 4.698e-07

### TBI
TBI group
Call:
lm(formula = int ~ cluster, data = filter(df, group == "TBI"))

Residuals:
     Min       1Q   Median       3Q      Max 
-9489379 -1643052  -388204   805481 15075286 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) -1314142    1409040  -0.933    0.358    
cluster        18166       1448  12.549 4.13e-14 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 4236000 on 33 degrees of freedom
Multiple R-squared:  0.8268,	Adjusted R-squared:  0.8215 
F-statistic: 157.5 on 1 and 33 DF,  p-value: 4.131e-14

## LM comparison by Z-transform
*p*=0.35