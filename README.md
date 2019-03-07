# Poisson-rates-whole-biotas

This repository provides R functions below implement different forms of the likelihood function used by [Bacon et al. 2015 Proc. Nat. Acad. Sci.](https://www.pnas.org/content/112/19/6110). This function models a series of dated events (e.g. divergence times, dated migration events) as arising from a Poisson process and mesure the most likely (Poisson) rate at which these events happened. This model is typically useful for studying the rate at which some events (e.g. divergence, dispersal) happened across an entire biota (i.e. in widely different clades) rather than within a single clade for which we have a phylogeny. 

The derivation of this likelihood function, as well as some important discussion of its limitations, is described in the Supplementary Information of the paper by Bacon et al. (section 1.6) and is due to the authors of this paper. I'm only responsible for this R implementation of this function.
