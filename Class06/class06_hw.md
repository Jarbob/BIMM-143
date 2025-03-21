# class06_hw
Mason Lew (PID: A17533139)

> Q6 How would you generalize the original code above to work with any
> set of input protein structures?

To do this question I will make a function using some of the tools found
in the `bio3d` library

``` r
#import tools from the bio3d library
library(bio3d)

#Creates the function named 'ploot', with argument 'input'
Ploot <- function(input) {
  
  #uses read.pdb on the input and stores the protein information in 'S0"
  S0 <- read.pdb(input)
  
  #Creates a smaller PDB object with a subset of atoms 
  S0.chainA <- trim.pdb(S0, chain="A", elety="CA")
  
  #Stores specific column of values within S0.b to be graphed
  S0.b <- S0.chainA$atom$b
  
  #Plots the data, outputs a graph comparing residue (x) and Bfactor (Y)
  plotb3(S0.b, sse=S0.chainA, typ="l", ylab="Bfactor")
}

#Calls the function with `Ploot()`, 
#The argument is the protein accession number for your protein of interest
#Protein 4AKE has been used as an example
Ploot("4AKE")
```

      Note: Accessing on-line PDB file

![](class06_hw_files/figure-commonmark/unnamed-chunk-1-1.png)
