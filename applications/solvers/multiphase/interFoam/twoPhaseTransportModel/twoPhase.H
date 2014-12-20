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
    twoPhase

Description

SourceFiles
    twoPhase.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhase_H
#define twoPhase_H

#include "incompressible/transportModel/transportModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoPhase Declaration
\*---------------------------------------------------------------------------*/

class twoPhase
:
    public transportModel
    
{
    // Private data

        word phase1Name_;
        word phase2Name_;

        autoPtr<transportModel> phase1_;
        autoPtr<transportModel> phase2_;

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;

        const volScalarField& gamma_;

        volScalarField nu_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();


public:

    //- Runtime type information
    TypeName("twoPhase");


    // Constructors

        //- construct from components
        twoPhase
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& phaseName = ""
        );


    // Destructor

        ~twoPhase()
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
            calcNu();
        }

        //- Read base transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
