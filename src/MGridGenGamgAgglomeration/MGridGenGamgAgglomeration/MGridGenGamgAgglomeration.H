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
    MGridGenGamgAgglomeration

Description
    Agglomerate using the MGridGen algorithm.

SourceFiles
    MGridGenGamgAgglomeration.C
    MGridGenGamgAgglomerate.C

\*---------------------------------------------------------------------------*/

#ifndef MGridGenGamgAgglomeration_H
#define MGridGenGamgAgglomeration_H

#include "gamgAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class MGridGenGamgAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class MGridGenGamgAgglomeration
:
    public gamgAgglomeration
{
    // Private data

        const fvMesh& fvMesh_;


    // Private Member Functions

        //- Construct the CSR format addressing
        void makeCompactCellFaceAddressingAndFaceWeights
        (
            const lduAddressing& fineAddressing,
            labelList& cellCells,
            labelList& cellCellOffsets,
            const vectorField& Si,
            List<scalar>& faceWeights
        );

        //- Calculate and return agglomeration
        tmp<labelField> agglomerate
        (
            label& nCoarseCells,
            const label minSize,
            const label maxSize,
            const lduAddressing& fineAddressing,
            const scalarField& V,
            const vectorField& Sf,
            const scalarField& Sb
        );


        //- Disallow default bitwise copy construct
        MGridGenGamgAgglomeration(const MGridGenGamgAgglomeration&);

        //- Disallow default bitwise assignment
        void operator=(const MGridGenGamgAgglomeration&);


public:

    //- Runtime type information
    TypeName("MGridGen");


    // Constructors

        //- Construct given mesh and controls
        MGridGenGamgAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
