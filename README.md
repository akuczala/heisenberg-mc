# heisenberg-mc
Monte Carlo stochastic series expansion (MC-SSE) simulation of the [Heisenberg](https://en.wikipedia.org/wiki/Heisenberg_model_(quantum)) chain or 2d lattice. Returns an estimate of the expected energy at the specified temperature. See [here](http://physics.bu.edu/~sandvik/programs/ssebasic/ssebasic.html) for algorithm details.

## Input file format

```
<dim> <itSteps> <states> <prec> <nSites> <nTemps>
<space separated list of system sizes>
<space separated list of temperatures>
```

where

- dim = Dimension of the spin system (1 = chain, 2 = lattice)
- itSteps = number of monte carlo steps
- states = Number of states in the SSE sequence
- prec = Convergence precision. Algorithm stops if the estimated energy changes by less than prec between iteration
- nSites = Number of system sizes to run algorithm with
- nTemps = Number of temperatures to run algorithm with

An example input file is included.
