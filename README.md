Introduction
============

This OpenFOAM-combo-pre15x repository was created as a complement to the unofficial OpenFOAM-combo repository, which tracks the git history of 1.5.x and up: https://github.com/wyldckat/OpenFOAM-combo

The objective is to have a single repository to track the changes that occurred since OpenFOAM 1.1 up to 1.5, without the git history for 1.5.x. This is because a git repository was not officially created for these versions and because the original _tarballs_ have code from third-parties.

Therefore, this repository will attempt to not include any source code that is unrelated to OpenFOAM itself, such as MICO (MICO Is COrba) and LAM/MPI, as well as the Doxygen documentation, in an effort to have a more easily browseable source code history of OpenFOAM itself.

*NOTE:* If you find any files that should be removed from this repository, because of source code that doesn't belong to OpenFOAM's own source code, please report this on the bug tracker of this project, so that the git tree can be properly pruned and rebuilt, namely here: https://github.com/wyldckat/OpenFOAM-combo-pre15x/issues

*Disclaimer*: This offering is not approved or endorsed by OpenCFD Limited, the producer of the OpenFOAM software and owner of the OPENFOAM®  and OpenCFD® trade marks.

This [git repository](https://github.com/wyldckat/OpenFOAM-combo-pre15x) was brought to you by Bruno Santos (wyldckat@github working at [blueCAPE Lda](http://www.bluecape.com.pt)).


How to use
==========

In order to only browse online, have a look into the [`combo-pre15x` branch](https://github.com/wyldckat/OpenFOAM-combo-pre15x/tree/combo-pre15x).

In order to have your own local copy, run these commands:
```
git clone https://github.com/wyldckat/OpenFOAM-combo-pre15x.git
cd OpenFOAM-combo-pre15x
git checkout combo-pre15x
```

Now, locally, if you want to checkout how `simpleFoam` has evolved over the versions, you can use `gitk` like this:
```
gitk applications/solvers/incompressible/simpleFoam/
```


Creating a repository similar to this one
=========================================

The original source code tarballs are available in the official project at SourceProject.Net: http://sourceforge.net/projects/foam/files/

The tarballs used to create this repository:
 * OpenFOAM-1.1.General.gtgz
 * OpenFOAM-1.2.General.gtgz
 * OpenFOAM-1.3.General.gtgz
 * OpenFOAM-1.4.General.gtgz
 * OpenFOAM-1.4.1.General.gtgz
 * OpenFOAM-1.5.General.gtgz


 In a nutshell, run these commands:
```
mkdir OpenFOAM-combo-pre15x
git init

tar -xf ../OpenFOAM-1.1.General.gtgz
mv OpenFOAM-1.1/* .
mv OpenFOAM-1.1/.[a-zA-Z]* .
rmdir OpenFOAM-1.1
find . -name lnInclude | xargs rm -r
rm -rf applications/utilities/mesh/conversion/ccm24ToFoam/libccmio*
rm applications/utilities/mesh/conversion/ccm24ToFoam/ccm24ToFoam/MeshCCMIO.* 
rm applications/utilities/mesh/conversion/ccm24ToFoam/ccm24ToFoam/Mesh.* 
rm applications/utilities/mesh/conversion/ccm24ToFoam/ccm24ToFoam/MeshNGEOM*
rm applications/utilities/mesh/conversion/ccm24ToFoam/ccm24ToFoam/polyfile.*
rm -rf src/lam-7.1.1
rm -rf src/malloc
rm -rf src/mico-2.3.11
rm -rf src/zlib-1.2.1
rm -rf applications/utilities/parallelProcessing/decompositionMethods/metisDecomp/metis-4.0
rm applications/utilities/surface/surfaceCoarsen/bunnylod/rab*
rm -rf doc/Doxygen/html
rm -rf ./applications/utilities/mesh/manipulation/patchTool/lib
rm -rf ./applications/utilities/preProcessing/FoamX/lib
find . -name "wmkdep" | xargs rm
find . -name "dirToString" | xargs rm

git add * .bashrc .cshrc .timeStamp .OpenFOAM-1.1
git commit -m "Source code for OpenFOAM 1.1, without third-party source code and generated data ('lnInclude', 'html' and '*.jar')."


rm -rf * .bashrc .cshrc .timeStamp .OpenFOAM-1.1

tar -xf ../OpenFOAM-1.2.General.gtgz
mv OpenFOAM-1.2/* .
mv OpenFOAM-1.2/.[a-zA-Z]* .
rmdir OpenFOAM-1.2
find . -name lnInclude | xargs rm -r
rm -rf applications/utilities/mesh/conversion/ccm24ToFoam/libccmio*
rm applications/utilities/mesh/conversion/ccm24ToFoam/ccm24ToFoam/MeshCCMIO* 
rm applications/utilities/mesh/conversion/ccm24ToFoam/ccm24ToFoam/Mesh.* 
rm -rf src/lam-7.1.1
rm -rf src/malloc
rm -rf src/mico-2.3.11
rm -rf src/zlib-1.2.1
rm -rf applications/utilities/parallelProcessing/decompositionMethods/metisDecomp/metis-4.0
rm applications/utilities/surface/surfaceCoarsen/bunnylod/rab*
rm -rf doc/Doxygen/html
rm -rf ./applications/utilities/mesh/manipulation/patchTool/lib
rm -rf ./applications/utilities/preProcessing/FoamX/lib
rm -rf applications/utilities/mesh/manipulation/setSet/readline-5.0
find . -name "wmkdep" | xargs rm
find . -name "dirToString" | xargs rm

git add -u
git add * .bashrc .cshrc .timeStamp .OpenFOAM-1.2
git commit -m "Source code for OpenFOAM 1.2, without third-party source code and generated data ('lnInclude', 'html' and '*.jar')."


rm -rf * .bashrc .cshrc .timeStamp .OpenFOAM-1.2

tar -xf ../OpenFOAM-1.3.General.gtgz
mv OpenFOAM-1.3/* .
mv OpenFOAM-1.3/.[a-zA-Z]* .
rmdir OpenFOAM-1.3
find . -name lnInclude | xargs rm -r
rm -rf applications/utilities/mesh/conversion/ccm24ToFoam/libccmio*
rm applications/utilities/mesh/conversion/ccm24ToFoam/ccm24ToFoam/MeshCCMIO* 
rm applications/utilities/mesh/conversion/ccm24ToFoam/ccm24ToFoam/Mesh.* 
rm -rf src/lam-7.1.1
rm -rf src/openmpi-1.0.2a7
rm -rf src/malloc
rm -rf src/mico-2.3.11
rm -rf src/zlib-1.2.1
rm -rf applications/utilities/parallelProcessing/decompositionMethods/metisDecomp/metis-4.0
rm applications/utilities/surface/surfaceCoarsen/bunnylod/rab*
rm -rf doc/Doxygen/html
rm -rf ./applications/utilities/mesh/manipulation/patchTool/lib
rm -rf ./applications/utilities/preProcessing/FoamX/lib
rm -rf applications/utilities/mesh/manipulation/setSet/readline-5.0
find . -name "wmkdep" | xargs rm
find . -name "dirToString" | xargs rm

git add -u
git add * .bashrc .cshrc .timeStamp .OpenFOAM-1.3
git commit -m "Source code for OpenFOAM 1.3, without third-party source code and generated data ('lnInclude', 'html' and '*.jar')."


rm -rf * .bashrc .cshrc .timeStamp .OpenFOAM-1.3

tar -xf ../OpenFOAM-1.4.General.gtgz
mv OpenFOAM-1.4/* .
mv OpenFOAM-1.4/.[a-zA-Z]* .
rmdir OpenFOAM-1.4
find . -name lnInclude | xargs rm -r
rm -rf applications/utilities/mesh/conversion/ccm26ToFoam/libccmio*
rm applications/utilities/mesh/conversion/ccm26ToFoam/ccm26ToFoam/MeshCCMIO* 
rm applications/utilities/mesh/conversion/ccm26ToFoam/ccm26ToFoam/Mesh.* 
rm -rf src/lam-7.1.2
rm -rf src/openmpi-1.2b3
rm -rf src/malloc
rm -rf src/mico-2.3.12
rm -rf src/zlib-1.2.1
rm applications/utilities/surface/surfaceCoarsen/bunnylod/rab*
rm -rf doc/Doxygen/html
rm -rf ./applications/utilities/mesh/manipulation/patchTool/lib
rm -rf ./applications/utilities/preProcessing/FoamX/lib
rm applications/test/Field/OpenFOAM.out
rm -rf applications/utilities/mesh/manipulation/setSet/readline-5.0
rm -rf applications/utilities/parallelProcessing/decompositionMethods/parMetisDecomp/ParMetis-3.1
rm -rf ./src/MGridGenGamgAgglomeration/ParMGridGen-1.0
find . -name "wmkdep" | xargs rm
find . -name "dirToString" | xargs rm

git add -u
git add * .bashrc .cshrc .timeStamp .OpenFOAM-1.4
git commit -m "Source code for OpenFOAM 1.4, without third-party source code and generated data ('lnInclude', 'html' and '*.jar')."


rm -rf * .bashrc .cshrc .timeStamp .OpenFOAM-1.4

tar -xf ../OpenFOAM-1.4.1.General.gtgz
mv OpenFOAM-1.4.1/* .
mv OpenFOAM-1.4.1/.[a-zA-Z]* .
rmdir OpenFOAM-1.4.1
find . -name lnInclude | xargs rm -r
rm -rf applications/utilities/mesh/conversion/ccm26ToFoam/libccmio*
rm -rf src/lam-7.1.2
rm -rf src/openmpi-1.2.3
rm -rf src/malloc
rm -rf src/mico-2.3.12
rm -rf src/zlib-1.2.1
rm -rf applications/utilities/parallelProcessing/decompositionMethods/metis-5.0pre2
rm -rf src/MGridGenGamgAgglomeration/ParMGridGen-1.0
rm applications/utilities/surface/surfaceCoarsen/bunnylod/rab*
rm -rf doc/Doxygen/html
rm -rf ./applications/utilities/mesh/manipulation/patchTool/lib
rm -rf ./applications/utilities/preProcessing/FoamX/lib
rm -rf applications/utilities/parallelProcessing/decompositionMethods/parMetisDecomp/ParMetis-3.1
rm -rf applications/utilities/parallelProcessing/decompositionMethods/metis-4.0
find . -name "wmkdep" | xargs rm
find . -name "dirToString" | xargs rm

git add -u
git add * .bashrc .cshrc .timeStamp .OpenFOAM-1.4.1
git commit -m "Source code for OpenFOAM 1.4.1, without third-party source code and generated data ('lnInclude', 'html' and '*.jar')."


rm -rf * .bashrc .cshrc .timeStamp .OpenFOAM-1.4.1

tar -xf ../OpenFOAM-1.5.General.gtgz
mv OpenFOAM-1.5/* .
mv OpenFOAM-1.5/.[a-zA-Z]* .
rmdir OpenFOAM-1.5
find . -name lnInclude | xargs rm -r
rm applications/utilities/surface/surfaceCoarsen/bunnylod/rab*
rm -rf doc/Doxygen/html
gzip tutorials/snappyHexMesh/motorBike/constant/triSurface/motorBike.stl
find . -name "wmkdep" | xargs rm
find . -name "dirToString" | xargs rm

git add -u
git add * .timeStamp .gitignore
git commit -m "Source code for OpenFOAM 1.5, without third-party source code and generated data ('lnInclude', 'html' and '*.jar')."
```
