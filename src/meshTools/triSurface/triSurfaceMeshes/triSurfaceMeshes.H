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

Description
    Container for triSurfaces read from files.

SourceFiles
    triSurfaceMeshes.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceMeshes_H
#define triSurfaceMeshes_H

#include "triSurfaceMesh.H"
#include "fileNameList.H"
#include "treeDataTriSurface.H"
#include "indexedOctree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class triSurfaceMeshes Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceMeshes
:
    public PtrList<triSurfaceMesh>
{
    // Private data

        //- Names of surfaces
        fileNameList names_;

        //- Trees
        mutable autoPtr<PtrList<indexedOctree<treeDataTriSurface> > > treesPtr_;


    // Private Member Functions

        void calcTrees() const;

        //- Disallow default bitwise copy construct
        triSurfaceMeshes(const triSurfaceMeshes&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceMeshes&);


public:

    // Constructors

        //- Construct from directory (io.instance()) and list of local filenames
        triSurfaceMeshes(const IOobject& io, const fileNameList&);


    // Member Functions

        //- Get all surfaces in directory
        static fileNameList allNames(const IOobject&);

        //- Get current list of surfaces
        const fileNameList& names() const
        {
            return names_;
        }

        //- Get trees for searching
        const PtrList<indexedOctree<treeDataTriSurface> >& trees() const
        {
            if (!treesPtr_.valid())
            {
                calcTrees();
            }
            return treesPtr_();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
