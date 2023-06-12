# PAKDD_Graph-Mining - README

## Abstract of the paper

Mining patterns in a dynamic attributed graph has received
more and more attention recently. However, it is a complex task because
both graph topology and attributes values of each vertex can change
over time. In this work, we focus on the discovery of frequent sequential
subgraph evolutions (FSSE) in such a graph. These FSSE patterns rep-
resent frequent evolutions of general sets of connected vertices’ attribute
values. A novel algorithm, named FSSEMiner, is proposed to mine FSSE
patterns. This algorithm is based on a new strategy (graph addition) to
guarantee mining efficiency. Experiments performed on benchmark and
real-world datasets show the interest of our approach and its scalability.

## Summary of the experimental evaluation

We run two experiments with the ffollowing three objectives : 

+ to evaluate the efficiency of our proposed algorithm to mine the datasets by comparing its performances in terms of number attributes, attributes values and size of size-1 patterns
+ to evaluate the scalability of our proposed algorithm by comparing its  performances on different numbers of temporal graphs
+ to evaluate the accuracy and pertinence of the results

## Technical environnement

To execute this code, you need : 

+ A c++ envrionnement using g++
+ Openmp installed 

To compile the code, run the following command : 

```bash
g++  -g "/your/path/to/Source.cpp" -fopenmp -o "/your/path/to/Source.out"
```

Then execute *Source.out*.
