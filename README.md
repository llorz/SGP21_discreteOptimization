
# Discrete Optimization for Shape Matching
This is an example code for the paper "Discrete Optimization for Shape Matching" by Jing Ren, Simone Melzi, Peter Wonka, and Maks Ovsjanikov.

In this paper we propose a *discrete solver* that can be used to optimize different functional map based energies. Specifically, for a given functional map based energy $\alpha$

<p align="center">
  <img align="center"  src="/figs/overview.png", width=800>
</p>



## Example 01: Minimize the area-preserving energy and the angle-preserving energy
<p align="center">
  <img align="center"  src="/figs/eg_box.png", width=800>
</p>

- The script ```eg1_box.m``` reproduce this result. 


## Example 02: Minimize descriptor-preserving energy
<p align="center">
  <img align="center"  src="/figs/eg_sphere.png", width=800>
</p>

- The script ```eg2_sphere.m``` reproduce this result. 

## Example 03:
<p align="center">
  <img align="center"  src="/figs/eg_faust.png", width=800>
</p>

- The script ```eg3_faust.m``` reproduce this result. 

## Example 04:
<p align="center">
  <img align="center"  src="/figs/eg_smal.png", width=800>
</p>

- The script ```eg4_smal.m``` reproduce this result. 

## Comments
- The script ```eg2_shapePair.m``` shows how to find multiple high-quality maps between a shape pair.
- Please let us know (jing.ren@kaust.edu.sa) if you have any question regarding the algorithms/paper ʕ•ﻌ•ʔ or you find any bugs in the code ԅ(¯﹃¯ԅ)


[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

This work is licensed under a [Creative Commons Attribution-NonCommercial 4.0 International License](http://creativecommons.org/licenses/by-nc/4.0/). For any commercial uses or derivatives, please contact us (jing.ren@kaust.edu.sa, peter.wonka@kaust.edu.sa, maks@lix.polytechnique.fr).
