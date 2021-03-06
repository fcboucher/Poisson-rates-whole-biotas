# Poisson-rates-whole-biotas

This repository provides R functions to implement different forms of the likelihood function used by [Bacon et al. 2015 Proc. Nat. Acad. Sci.](https://www.pnas.org/content/112/19/6110). This function models a series of dated events (e.g. divergence times, dated migration events, etc.) as arising from a Poisson process and measure the most likely (Poisson) rate at which these events happened. This model is typically useful for studying the rate at which some events (e.g. divergence, dispersal, etc.) happened across an entire biota (i.e. in widely different clades) rather than within a single clade for which we have a well sampled phylogeny. 

The derivation of this likelihood function, as well as some important discussion of its limitations, is described in the Supplementary Information of the paper by Bacon et al. (section 1.6) and is due to the authors of this paper. I'm only responsible for this R implementation of this function, which was used in the following study: [Husson et al. 2020 Journal of Biogeography](https://onlinelibrary.wiley.com/doi/full/10.1111/jbi.13762).

The repository contains two R scripts. The first one (Poisson_rate_...) has the likelihood functions written, it is heavily commented so that you understand what is going on. The other one (Example_script.r) contains an example script to show you how to use these functions.
