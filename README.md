# MSbox

Common MS tools for analytical chemistry

# Installation 

```r
## stable version
install.packages('MSbox')
```

# Functions

1. Check element isotopes

examples:

```r
E_iso('C') # element symbol, case insensitive
E_iso('Carbon') # element full name, case insensitive
E_iso('carBon') # element full name, case insensitive
```

2. Calculate monoisitopic mass

example:

```r
M_mass('C7H6O1') # case insensitive
```

3. Calculate exact m/z values

example:

```r
mz('C7H7O4', z = -1) # case insensitive
```

4. Calculate the mass accuracy of measured m/z

examples:

```r
ppm(155.03383, 155.03388) # with m/z value
ppm(155.03383, mz('C7H7O4', z = -1)) # with ion formula, case insensitive
```

5. Check if an m/z value originates from possible contaminant

examples

```r
contam(33.0335, ppm = 10, mode = '+')
contam(44.998, ppm = 10, mode = '-')
```

6. Calculate common adduct ions in pos or neg ion mode

examples

```r
adduct('C1H4',mode = '-') # case insensitive
adduct('C1H4',mode = '+') # case insensitive
```

7. Calculate tanimoto similarities of different compounds

example

```r
x <- data.frame(Samp1=c(0,0,0,1,1,1,0,0,1), Samp2=c(1,1,1,1,1,1,0,0,1))
tanimoto(x)
```
