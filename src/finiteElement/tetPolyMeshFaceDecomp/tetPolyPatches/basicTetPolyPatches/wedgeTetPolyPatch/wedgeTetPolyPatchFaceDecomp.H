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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    wedgeTetPolyPatchFaceDecomp

Description
    Wedge front and back plane patch.

SourceFiles
    wedgeTetPolyPatchFaceDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeTetPolyPatchFaceDecomp_H
#define wedgeTetPolyPatchFaceDecomp_H

#include "faceTetPolyPatchFaceDecomp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class wedgeTetPolyPatchFaceDecomp Declaration
\*---------------------------------------------------------------------------*/

class wedgeTetPolyPatchFaceDecomp
:
    public faceTetPolyPatchFaceDecomp
{

public:

    //- Runtime type information
    TypeName("wedge");


    // Constructors

        //- Construct from polyPatch
        wedgeTetPolyPatchFaceDecomp
        (
            const polyPatch& patch,
            const tetPolyBoundaryMeshFaceDecomp& bm
        )
        :
            faceTetPolyPatchFaceDecomp(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
