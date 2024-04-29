# Fabric Tessellation: Realizing Freeform Surfaces by Smocking
<p align="center">
<img src="./figs/teaser.jpg" width="800" /> 
</p>

This repository contains the code for the paper "Fabric Tessellation: Realizing Freeform Surfaces by Smocking", accetped to SIGGRAPH 2024, by [Aviv Segall](https://segaviv.github.io), [Jing Ren](https://ren-jing.com/), [Amir Vaxman](https://avaxman.github.io/) and [Olga Sorkine-Hornung](https://igl.ethz.ch/people/sorkine). 

Smocking is renowned for producing intricate geometric textures with voluminous pleats. However, it has been mostly used to realize flat shapes or manually designed, limited classes of curved surfaces (e.g. see [Canadian smocking](https://github.com/llorz/SmockingDesign) and [Italian smocking](https://github.com/nifzhou/ItalianSmocking)). 
In this project, we aim at solving the inverse design problem for smocking, called fabric tessellation.
We present a novel method for realizing freeform surfaces with pieces of flat fabric, where curvature is created by stitching together points on the fabric using the smocking technique. 

You can find more details at (coming soon): [[project page]]() | [[paper]]() | [[suppl. video]]() | [[slides]]()

## Methodology
<p align="center">
<img src="./figs/algorithm.jpg" width="1000" /> 
</p>

### Problem Formulation
- ***Input***: a 3D triangle mesh (with arbitrary topology and discretization) and desired smocking pleat type (such as the "arrow" pattern specified in the figure above)
- ***Output***: a modified 2D smocking pattern that, after fabrication where the specified points are sewn together, results in a textile that approximates the target shape and exhibits visually pleasing pleats (as shown in the teaser figure).

### Algorithm
Our method combines the computation of directional fields with continuous optimization of a Tangram graph in the plane, which together allow us to realize surfaces of arbitrary topology and curvature with smocking patterns of diverse symmetries. More specifically:
1. extract a so-called **Tangram graph** from the input smocking pattern : see Sec.4 and Appendix B for more detailed definitions and discussions regarding Tangram and its construction.
2. compute the **closed configuration** of the Tangram: the closing process of Tangram explains how the 3D pleats are formed during the smocking process (please see Fig.6 or [[web demo]]() for illustrations).
3. the closed Tangram is pulled-back to the input 3D shape to fetch the expected lengths (curvatures): this can also be interpreted as reparameterizing or tiling the input 3D surface using the closed Tangram structure.
4. optimize the Tangram in 2D in its **open configuration** to realize the expected lengths (i.e, once the optimized Tangram is closed, it will reproduce the input shape)
5. extract the 2D smocking pattern from Tangram since they have one-to-one correspondence

## Implementation
### main functions
blablabla
### how to use
blablabla 


## Comments
### Acknowledgements
The authors would like to thank the anonymous reviewers for their valuable feedback. 
This work was supported in part by the ERC Consolidator Grant No. 101003104 (MYCLOTH).
Special thanks to ***Ningfeng Zhou*** for her assitance in fabricating the heart and cloud shapes, and to ***all members of IGL*** for their insightful discussions and kind support. 

### Contact
Please let us know (aviv.segall, jing.ren, @inf.ethz.ch) if you have any question regarding the algorithms/paper or you find any bugs in the code ε-(´∀｀; )
This work is licensed under a [Creative Commons Attribution-NonCommercial 4.0 International License](http://creativecommons.org/licenses/by-nc/4.0/). 
For any commercial uses or derivatives, please contact us (aviv.segall, jing.ren, sorkine, @inf.ethz.ch, avaxman, @inf.ed.ac.uk). [![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
