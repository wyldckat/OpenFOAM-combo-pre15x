Introduction
============

This OpenFOAM-combo-pre15x repository was created as a complement to the unofficial OpenFOAM-combo repository, which tracks the git history of 1.5.x and up: https://github.com/wyldckat/OpenFOAM-combo

The objective is to have a single repository to track the changes that occurred since OpenFOAM 1.1 up to 1.5, without the git history for 1.5.x. This is because a git repository was not officially created for these versions and because the original _tarballs_ have code from third-parties.

Therefore, this repository will attempt to not include any source code that is unrelated to OpenFOAM itself, such as MICO (MICO Is COrba) and LAM/MPI, in an effort to have a more easily browseable source code history of OpenFOAM itself.

*Disclaimer*: This offering is not approved or endorsed by OpenCFD Limited, the producer of the OpenFOAM software and owner of the OPENFOAM®  and OpenCFD® trade marks.


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

