# DomainDLRS
DomainDLRS is a hierarchical, generative probabilistic model containing three levels corresponding to species, genes, and domains, respectively. It simultanously reconstructs the domain and gene evlulitonay histories in the light of a given species tree. It is implemented in Java. It takes a dated species tree together with a multiple sequence alignment for each domain family as input and outputs an estimated posterior distribution over reconciled gene and domain trees.

# Usage:
```
#! /bin/bash
java -jar domainDLRS.jar -p <ParameterFile>
```

# Requirements:
1. Efficient Java Matrix Library (EJML) <http://www.ejml.org>
2. Commons Math: The Apache Commons Mathematics Library <http://commons.apache.org/proper/commons-math/>
3. Java Evolutionary Biology Library <https://github.com/Biomatters/JEBL>
4. MersenneTwisterRNG Random Number Generator from Uncommons Maths <https://maths.uncommons.org>

# Citation:
Sayyed Auwn Muhammad, Bengt Sennblad, and Jens Lagergren. Species tree aware simultaneous reconstruction of gene and domain evolution., 2018.
