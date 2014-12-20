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
    polyMeshReconstructor

Description
    Mesh reconstructor.

SourceFiles
    polyMeshReconstructor.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshReconstructor_H
#define polyMeshReconstructor_H

#include "ptrList.H"
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
        const ptrList<polyMesh>& procMeshes_;

        //- List of processor point addressing lists
        const ptrList<labelIOList>& pointProcAddressing_;

        //- List of processor face addressing lists
        const ptrList<labelIOList>& faceProcAddressing_;

        //- List of processor cell addressing lists
        const ptrList<labelIOList>& cellProcAddressing_;

        //- List of processor boundary addressing lists
        const ptrList<labelIOList>& boundaryProcAddressing_;


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
            const ptrList<polyMesh>& procMeshes,
            const ptrList<labelIOList>& pointProcAddressing,
            const ptrList<labelIOList>& faceProcAddressing,
            const ptrList<labelIOList>& cellProcAddressing,
            const ptrList<labelIOList>& boundaryProcAddressing
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
