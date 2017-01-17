Description
===========
RigidMotionsMapleTools is a set of Maple modules for which the main objective is to generate unique
neighborhood motion maps (NMM, for short). These neighborhood motion maps are the images of an image
patch which is a finite set of integer points. More information about the alpha version of the
algorithm were published in: [Pluta, Moroz, Kenmochi, Romon,Quadric Arrangement in Classifying Rigid
Motions of a 3D Digital Image, CASC16, Springer
2016](http://link.springer.com/chapter/10.1007%2F978-3-319-45641-6_27). You can obtain a preprint
from [HAL](https://hal.archives-ouvertes.fr/hal-01334257/document) for free.

Quick Install
=============

1. Download or clone the repository.
2. [Edit (or create) Maple initialization
   file](https://www.maplesoft.com/support/help/Maple/view.aspx?path=worksheet/reference/initialization)
   and add to it paths to the files (the order matters): 
   
   1. RealAlgebraicNumber.mpl
   2. RigidMotionsParameterSpaceEvent.mpl
   3. RigidMotionsParameterSpaceRegister.mpl
   4. RigidMotionsParameterSpaceCommon.mpl
   5. RigidMotionsParameterSpaceDecompositionRecursive.mpl
   6. RigidMotionsParameterSpaceDecomposition.mpl
   7. RigidMotionsRecoverNMM.mpl
3. Additionally, FGb library by Jean-Charles Faugère can be installed to improve generation of
   univariate polynomials. For more information see:
   [http://www-polsys.lip6.fr/~jcf/FGb/index.html](http://www-polsys.lip6.fr/~jcf/FGb/index.html).
   **We do recommend to install FGb library in version 1.61.**


Examples
================

To compute rotational sample points for 6-neighborhood (for documentation of the parameters see the
source file):
 
```
with(RigidMotionsParameterSpaceDecompostion);
LaunchComputeSamplePoints([a,b,c], "N1.db", "N1", [-1, 0, 1]);
```

Resume computations after a crash (for documentation of the parameters see the
source file):

```
with(RigidMotionsParameterSpaceDecompostion);
LaunchResumeComputations([a,b,c], "N1.db", "N1");
```

Compute neighborhood motion maps (for documentation of the parameters see the
source file):

```
with(RigidMotionsRecoverNMM);
LaunchComputeNMM([a,b,c], "N1", [-1,0,1], "N1.db");
```

Additional information
================

The computations are memory and time expensive. The implementation uses the Maple Grid framework to
distribute (locally) the computation on given number of nodes. While, the problem for 6-neighbor 
can be solved on a desktop machine in a relatively short time, computations for bigger image patch
can take weeks even for a machine with dozens of cores and hundreds of gigabytes of memory. Note that,
if the FGb library is installed then computations are relatively slower but less events are generated.

