Description
---------------------

Topological similarity estimation for 3D models (face-vertex meshes) using multiresolutional Reeb graphs (MRG).

The MRG method was proposed by Hilaga et al. in [1]. 
This implementation was used to obtain experimental results that were reported in [2] and [3]. 

### Implementation Notes

 - Programming language: **Java**
 - The code was developed using J2SE 1.4, so it does not use Java generics (introduced in J2SE 5.0)
 - **3D models must be specified in a specially-formatted VRML files** (see below)
   - Face-vertex meshes with triangle and quad faces are supported
   - 16 sample models are included in this package

Package Contents
---------------------

### Source Code

The source code is located in [`src/`](src/) directory. 
There are **two** main java classes:

* [`ExtractReebGraph`](src/ExtractReebGraph.java) constructs MRGs for 3D models and saves them into text files
 - MRG construction is described in Section 4 of [1]
* [`CompareReebGraph`](src/CompareReebGraph.java) implements matching algorithm for a pairwise comparison of MRGs
 - Matching algorithm is described in Section 5 of [1]

### Sample 3D Models 

The following 16 CAD models in VRML format can be found in [`models/`](models/) directory:

<a href="https://raw.github.com/dbespalov/reeb_graph/master/figs/sample_models.pdf"><img  width="300px" target="_blank" src="https://raw.github.com/dbespalov/reeb_graph/master/figs/sample_models.png"/></a>


VRML parser in [`ExtractReebGraph`](src/ExtractReebGraph.java) can **only** handle specially-formatted face-vertex meshes:

1. Vertex lists are assumed to be enclosed by strings `point [\n` and `]\n`
  * Coordinates for each point must appear on a separate line in the VRML file (e.g., `1.5 3.2 0.2,\n`)
2. Face lists are assumed to be enclosed by strings `coordIndex [\n` and `]\n`
  * Each face must appear on a separate line (e.g., `0, 2, 1, -1,\n` encodes a triangle face)
  

Usage
---------------------

License
---------------------
GNU General Public License


References
---------------------

1. Topology Matching for Fully Automatic Similarity Estimation of 3D Shapes.  
   Masaki Hilaga, Yoshihisa Shinagawa, Taku Kohmura, and Tosiyasu L. Kunii.  
   *SIGGRAPH*, 2001. [doi:10.1145/383259.383282](http://dx.doi.org/10.1145/383259.383282)                                                             

2. Reeb graph based shape retrieval for CAD.  
   Dmitriy Bespalov, William C. Regli, and Ali Shokoufandeh.  
   ASME IDETC, 2003. [doi:10.1115/DETC2003/CIE-48194](http://dx.doi.org/10.1115/DETC2003/CIE-48194)

3. Benchmarking search techniques for CAD.  
   Dmitriy Bespalov, Cheuk Yiu Ip, William C. Regli, and Joshua Shaffer.  
   *ACM SPM*, 2005. [doi:10.1145/1060244.1060275](http://dx.doi.org/10.1145/1060244.1060275)
