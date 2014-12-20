/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    triSurfaceMeshes

\*----------------------------------------------------------------------------*/

#include "triSurfaceMeshes.H"
#include "Random.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::triSurfaceMeshes::calcTrees() const
{
    // Random number generator
    Random rndGen(65431);

    treesPtr_.reset(new PtrList<indexedOctree<treeDataTriSurface> >(size()));
    PtrList<indexedOctree<treeDataTriSurface> >& trees = treesPtr_();

    forAll(*this, i)
    {
        const triSurface& s = operator[](i);

        // bb of surface
        treeBoundBox bb(s.localPoints());

        trees.set
        (
            i,
            new indexedOctree<treeDataTriSurface>
            (
                treeDataTriSurface(s),
                bb.extend(rndGen, 1E-3),    // slightly randomize bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::triSurfaceMeshes::triSurfaceMeshes
(
    const IOobject& io,
    const fileNameList& names
)
:
    PtrList<triSurfaceMesh>(names.size()),
    names_(names),
    treesPtr_(NULL)
{
    forAll(names, i)
    {
        autoPtr<IOobject> surfaceIO = io.clone();
        surfaceIO().rename(names[i]);

        Pout<< "Loading surface " << surfaceIO().filePath() << endl;

        set
        (
            i,
            new triSurfaceMesh
            (
                surfaceIO(),
                surfaceIO().filePath()
            )
        );

        string oldPrefix(Pout.prefix());
        Pout.prefix() += "    ";
        operator[](i).writeStats(Pout);
        Pout.prefix() = oldPrefix;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileNameList Foam::triSurfaceMeshes::allNames(const IOobject& io)
{
    return readDir(io.path(), fileName::FILE);
}


// ************************************************************************* //
