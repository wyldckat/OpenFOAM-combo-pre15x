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
    algebraicPairGamgAgglomeration

Description
    Agglomerate using the pair algorithm.

SourceFiles
    algebraicPairGamgAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef algebraicPairGamgAgglomeration_H
#define algebraicPairGamgAgglomeration_H

#include "pairGamgAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class algebraicPairGamgAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class algebraicPairGamgAgglomeration
:
    public pairGamgAgglomeration
{

public:

    //- Runtime type information
    TypeName("algebraicPair");


    // Constructors

        //- Construct given mesh and controls
        algebraicPairGamgAgglomeration
        (
            const lduMatrix& matrix,
            const dictionary& controlDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
