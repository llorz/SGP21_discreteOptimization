
# Discrete Optimization for Shape Matching
This is an example code for the paper "Discrete Optimization for Shape Matching" 
by [Jing Ren](https://ren-jing.com/), [Simone Melzi](https://sites.google.com/site/melzismn/), [Peter Wonka](http://peterwonka.net/), and [Maks Ovsjanikov](http://www.lix.polytechnique.fr/~maks/).


In this paper we propose a *discrete solver* that can be used to optimize different functional map based energies. Specifically, for a given functional map based energy E(C), our solver can minimize the input energy with the hard constraint that the functional map C is proper, i.e., corresponding to a pointwise map. 
<p align="center">
  <img align="center"  src="/figs/overview.png", width=800>
</p>



## Example 01: Minimize the area-preserving and the angle-preserving energy
<p align="center">
  <img align="center"  src="/figs/eg_box.png", width=800>
</p>

- In this example we use our solver to minimize the **area-preserving** energy E1 (othogonality) and **angle-preserving** energy E2 (conformality)
- We test on two shape pairs (S1, S2) and (S1, S3), where S2 has the same surface area as S1, and S3 is conformal to S1.
- We compare our discrete solver (D) to the standard continuous solver (C). 
- We can see that, using our discrete solver on minimizing the area-preserving energy on shape pair (S1, S2), and minimizing the angle-preserving energy on shape pair (S1, S3) leads to more desirable pointwise maps than the other setting. 
- The script ```eg1_box.m``` reproduce this result. 


## Example 02: Minimize descriptor-preserving energy
<p align="center">
  <img align="center"  src="/figs/eg_sphere.png", width=800>
</p>

- In this example, we use our solver to minimize the **descriptor-preserving** energy
- We show three examples. For each example, we minimize the functional map energy that enforce the correspondences between f and g. 
- We can see that final pointwise maps (last column) are consistent with the input corresponding descriptors. 
- The script ```eg2_sphere.m``` reproduce this result. 

## Example 03: Alternative to the multiplicative operator
<p align="center">
  <img align="center"  src="/figs/eg_faust.png", width=800>
</p>

- In this example, we compare two settings:
  - (1) [[NO17]](http://www.lix.polytechnique.fr/~maks/papers/fundescEG17.pdf): use a standard continuous solver to minimize the functional map energy that consists of the descriptor-preserving energy, Laplacian commutativity energy, and the **multiplicative operator commutativity** energy that are induced from the input descriptors and are proposed to make functional maps more proper (i.e., have an underlying pointwise map)
  - (2) Our setting: we used our discrete solver to minimize the same energy with the same parameters as (1) without the multiplicative term. 
- In this case, we can compare our discrete solver to the multiplicative term proposed in [NO17].
- The script ```eg3_faust.m``` reproduce this result. 

## Example 04: Minimize the Laplacian Commutativity
<p align="center">
  <img align="center"  src="/figs/eg_smal.png", width=800>
</p>

- In this example, we test our discrete solver on a set of SMAL animal shape pairs on minimizing the Laplacian Commutativity energy from random initialization. 
- We compare to the standard continuous solver. 
- The script ```eg4_smal.m``` reproduce this result. 

## Comments
- Please let us know (jing.ren@kaust.edu.sa) if you have any question regarding the algorithms/paper ʕ•ﻌ•ʔ or you find any bugs in the code ԅ(¯﹃¯ԅ)
- This work is licensed under a [Creative Commons Attribution-NonCommercial 4.0 International License](http://creativecommons.org/licenses/by-nc/4.0/). For any commercial uses or derivatives, please contact us (jing.ren@kaust.edu.sa, melzismn@gmail.com, peter.wonka@kaust.edu.sa, maks@lix.polytechnique.fr).
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)


