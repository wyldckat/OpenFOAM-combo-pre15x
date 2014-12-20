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
    gammaContactAngleFvPatchScalarField

Description

SourceFiles
    gammaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef gammaContactAngleFvPatchScalarField_H
#define gammaContactAngleFvPatchScalarField_H

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class gammaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class gammaContactAngleFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
    // Private data

        //- Equilibrium contact angle
        scalar theta0_;

        //- Dynamic contact angle velocity scale
        scalar uTheta_;

        //- Limiting advancing contact angle
        scalar thetaA_;

        //- Limiting receeding contact angle
        scalar thetaR_;


public:

    //- Runtime type information
    TypeName("gammaContactAngle");


    // Constructors

        //- Construct from patch and internal field
        gammaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const scalarField&
        );

        //- Construct from patch, internal field and dictionary
        gammaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const scalarField&,
            const dictionary&
        );

        //- Construct by mapping given gammaContactAngleFvPatchScalarField
        //  onto a new patch
        gammaContactAngleFvPatchScalarField
        (
            const gammaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const scalarField&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new gammaContactAngleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        gammaContactAngleFvPatchScalarField
        (
            const gammaContactAngleFvPatchScalarField&,
            const scalarField&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(const scalarField& iF) const
        {
            return tmp<fvPatchScalarField>
            (
                new gammaContactAngleFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Return the equilibrium contact angle theta0
        scalar theta0() const
        {
            return theta0_;
        }

        //- Return the dynamic contact angle velocity scale
        scalar uTheta() const
        {
            return uTheta_;
        }

        //- Return the limiting advancing contact angle
        scalar thetaA() const
        {
            return thetaA_;
        }

        //- Return the limiting receeding contact angle
        scalar thetaR() const
        {
            return thetaR_;
        }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
