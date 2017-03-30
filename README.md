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
2. ```cd <path to the cloned repository/scripts>```
3. ```./Install.sh -d="<path to the cloned repository>"```

Install
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

Quick Uninstall
=============
1. ```cd <path to the cloned repository/scripts>```
2. ```./Install.sh -u```

Examples
================

To compute a sorted list of events for 6-neighborhood (for documentation of the parameters see the
source file):
```
with(RigidMotionsParameterSpaceDecompostion);
LaunchComputeEvents([a,b,c], "N1.db", "N1", [-1, 0, 1]);
```

To run or resume computations (for documentation of the parameters see the source file):
```
with(RigidMotionsParameterSpaceDecompostion);
LaunchComputeSamplePoints([a,b,c], "N1.db", "N1", [-1, 0, 1]);
```

To compute neighborhood motion maps (for documentation of the parameters see the
source file):
```
with(RigidMotionsRecoverNMM);
LaunchComputeNMM([a,b,c], "N1", [-1,0,1], "N1.db");
```

Submitting jobs to a cluster
================

The directory ```scripts/cluster``` contains several scripts meant to easy the process of running
parts of the code in a cluster. The only part of the code which does not work in a cluster is
```LaunchComputeEvents```.

The main script ```Build.sh``` creates a self-executable bash script which contains all the code and
input data needed to run computation in a cluster. To build a self-executable script run:

``` ./Build.sh -d=<database.db> -n=<no. nodes> -o=<output.shx> -m=<scpript.mpl>```.

The input parameters are:

- -d -- path to a database which has to contain at least computed events
- -n -- number of nodes in the cluster
- -o -- where to write the output self-executable script
- -m -- a maple script which contains a set of instructions to be run on each node in the cluster
  (see an exemple of such a file: ```scripts/cluster/MapleScriptExample.mpl```)

To submit computations to the cluster you have to execute the self-executable script obtained from
```scripts/cluster/Build.sh``` and provide as an argument a path to a folder accessible via the same
path by each node in the cluster. The self-executable script will create a folder
```<SHARED_DIR>/selfextract.XXXXXX/DB```, where XXXXXX stands for some random string, where the
output is stored in a number (equal to the number of nodes provided to ```Build.sh``` via -n
parameter) of database files. For the moment the files have to merged by hand but a script to
facilitated this operation is planned.


Additional information
================

The computations are memory and time expensive. The implementation uses the Maple Grid framework to
distribute (locally) the computation on given number of nodes. While, the problem for 6-neighbor can
be solved on a desktop machine in a relatively short time, computations for bigger image patch can
take weeks even for a machine with dozens of cores and hundreds of gigabytes of memory. Note that,
if the FGb library is installed then computations are relatively slower but less events are
generated.

