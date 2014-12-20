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
    CrossPowerLaw

Description

SourceFiles
    CrossPowerLaw.C

\*---------------------------------------------------------------------------*/

#ifndef CrossPowerLaw_H
#define CrossPowerLaw_H

#include "transportModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace transportModels
{

/*---------------------------------------------------------------------------*\
                           Class CrossPowerLaw Declaration
\*---------------------------------------------------------------------------*/

class CrossPowerLaw
:
    public transportModel
{
    // Private data

        dictionary CrossPowerLawCoeffs_;

        dimensionedScalar nu0_;
        dimensionedScalar nuInf_;
        dimensionedScalar m_;
        dimensionedScalar n_;

        volScalarField nu_;

    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("CrossPowerLaw");


    // Constructors

        //- Construct from components
        CrossPowerLaw
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& phaseName = ""
        );


    // Destructor

        ~CrossPowerLaw()
        {}


    // Member Functions

        //- Return the laminar viscosity
        const volScalarField& nu() const
        {
            return nu_;
        }

        //- Correct the laminar viscosity
        void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace transportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
