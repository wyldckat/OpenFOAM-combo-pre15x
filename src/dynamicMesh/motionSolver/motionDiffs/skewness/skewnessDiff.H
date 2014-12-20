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
    skewnessDiff

Description
    Skewness-based motion diffusion.

SourceFiles
    skewnessDiff.C

\*---------------------------------------------------------------------------*/

#ifndef skewnessDiff_H
#define skewnessDiff_H

#include "motionDiff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class skewnessDiff Declaration
\*---------------------------------------------------------------------------*/

class skewnessDiff
:
    public motionDiff
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        skewnessDiff(const skewnessDiff&);

        //- Disallow default bitwise assignment
        void operator=(const skewnessDiff&);


public:

    //- Runtime type information
    TypeName("skewness");


    // Constructors

        //- Construct from components
        skewnessDiff
        (
            const tetPolyMesh& tetMesh,
            const dictionary&
        );


    // Destructor

        virtual ~skewnessDiff();


    // Member Functions

        //- Return diffusivity field
        virtual tmp<elementScalarField> gamma() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
