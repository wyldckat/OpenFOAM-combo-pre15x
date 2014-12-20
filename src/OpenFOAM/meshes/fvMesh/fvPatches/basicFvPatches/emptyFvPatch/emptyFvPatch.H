/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    emptyFvPatch

Description
    A patch which will not exist in the fvMesh. Typical example is a front and
    back plane of a 2-D geometry

SourceFiles
    emptyFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef emptyFvPatch_H
#define emptyFvPatch_H

#include "fvPatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class emptyFvPatch Declaration
\*---------------------------------------------------------------------------*/

class emptyFvPatch
:
    public fvPatch
{
    // Private data

        static const vectorField emptyNf_;


public:

    //- Runtime type information
    TypeName(emptyPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        emptyFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm)
        {}


    // Member Functions

        virtual label size() const
        {
            return 0;
        }

        //- Return face normals. Over-riding base class return to get zero size
        // 
        virtual const vectorField& nf() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
