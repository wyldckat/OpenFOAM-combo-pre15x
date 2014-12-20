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
    polyMeshReconstructor

Description
    Mesh reconstructor.

SourceFiles
    polyMeshReconstructor.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshReconstructor_H
#define polyMeshReconstructor_H

#include "PtrList.H"
#include "polyMesh.H"
#include "IOobjectList.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class polyMeshReconstructor Declaration
\*---------------------------------------------------------------------------*/

class polyMeshReconstructor
{
    // Private data

        //- Reconstructed mesh reference
        polyMesh& mesh_;

        //- List of processor meshes
        const PtrList<polyMesh>& procMeshes_;

        //- List of processor point addressing lists
        const PtrList<labelIOList>& pointProcAddressing_;

        //- List of processor face addressing lists
        const PtrList<labelIOList>& faceProcAddressing_;

        //- List of processor cell addressing lists
        const PtrList<labelIOList>& cellProcAddressing_;

        //- List of processor boundary addressing lists
        const PtrList<labelIOList>& boundaryProcAddressing_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        polyMeshReconstructor(const polyMeshReconstructor&);

        //- Disallow default bitwise assignment
        void operator=(const polyMeshReconstructor&);


public:

    // Constructors

        //- Construct from components
        polyMeshReconstructor
        (
            polyMesh& mesh,
            const PtrList<polyMesh>& procMeshes,
            const PtrList<labelIOList>& pointProcAddressing,
            const PtrList<labelIOList>& faceProcAddressing,
            const PtrList<labelIOList>& cellProcAddressing,
            const PtrList<labelIOList>& boundaryProcAddressing
        );


    // Member Functions

        //- Reconstruct point position after motion in parallel
        void reconstructPoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
