Description
---------------------

Topological similarity estimation for 3D models (face-vertex meshes) using multiresolutional Reeb graphs (MRG).

The MRG method was proposed by Hilaga et al. [1]. This implementation of MRG was used to obtain experimental results that were reported in [2] and [3]. 

### Implementation Notes

 - Programming language: **Java**
 - The code was developed using J2SE 1.4, so it does not use Java generics (introduced in J2SE 5.0)
 - **3D models must be specified in a specially-formatted VRML files** (see [below](3Dmodels) )
   - Face-vertex meshes with triangle and quad faces are supported
   - 16 sample models are included in this package

### References

1. Topology Matching for Fully Automatic Similarity Estimation of 3D Shapes.  
   Masaki Hilaga, Yoshihisa Shinagawa, Taku Kohmura, and Tosiyasu L. Kunii.  
   *SIGGRAPH*, 2001. [doi:10.1145/383259.383282](http://dx.doi.org/10.1145/383259.383282)                                                             

3. Reeb graph based shape retrieval for CAD.  
   Dmitriy Bespalov, William C. Regli, and Ali Shokoufandeh.  
   ASME IDETC, 2003. [doi:10.1115/DETC2003/CIE-48194](http://dx.doi.org/10.1115/DETC2003/CIE-48194)

2. Benchmarking search techniques for CAD.  
   Dmitriy Bespalov, Cheuk Yiu Ip, William C. Regli, and Joshua Shaffer.  
   *ACM SPM*, 2005. [doi:10.1145/1060244.1060275](http://dx.doi.org/10.1145/1060244.1060275)


Package Content
---------------------

### Source Code

[`src/`](src/)

[`src/ExtractReebGraph.java`](src/ExtractReebGraph.java)
[`src/CompareReebGraph.java`](src/CompareReebGraph.java)

### Sample 3D Models [3Dmodels]


Sample Usage
---------------------


License
---------------------
GNU General Public License
