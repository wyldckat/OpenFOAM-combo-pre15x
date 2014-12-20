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
    wallTetPolyPatchFaceDecomp

Description
    Wall patch.

SourceFiles
    wallTetPolyPatchFaceDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef wallTetPolyPatchFaceDecomp_H
#define wallTetPolyPatchFaceDecomp_H

#include "faceTetPolyPatchFaceDecomp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class wallTetPolyPatchFaceDecomp Declaration
\*---------------------------------------------------------------------------*/

class wallTetPolyPatchFaceDecomp
:
    public faceTetPolyPatchFaceDecomp
{

public:

    //- Runtime type information
    TypeName("wall");


    // Constructors

        //- Construct from polyPatch
        wallTetPolyPatchFaceDecomp
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
