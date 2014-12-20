/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006 Mark Olesen
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
    massFluxInletVelocityFvPatchVectorField

Description
    Describes an massflux normal vector boundary condition by its magnitude
    as an integral over its area.  The current density is used to correct the
    velocity.

    Example of the BC specification:

        inlet
        {
            type            massFluxInletVelocity;
            massFlux        0.2;     // mass flow [kg/s]
        }

NOTE!
 * may not work correctly for transonic inlets!
 * strange behaviour with potentialFoam since the U equation is not solved

SourceFiles
    massFluxInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef massFluxInletVelocityFvPatchVectorField_H
#define massFluxInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class massFluxInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class massFluxInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet mass flux
        scalar massFlux_;

        word phiName_;
        word rhoName_;


public:

    //- Runtime type information
    TypeName("massFluxInletVelocity");

    // Constructors

        //- Construct from patch and internal field
        massFluxInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const vectorField&
        );

        //- Construct from patch, internal field and dictionary
        massFluxInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const vectorField&,
            const dictionary&
        );

        //- Construct by mapping given
        //  massFluxInletVelocityFvPatchVectorField
        //  onto a new patch
        massFluxInletVelocityFvPatchVectorField
        (
            const massFluxInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const vectorField&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        massFluxInletVelocityFvPatchVectorField
        (
            const massFluxInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new massFluxInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        massFluxInletVelocityFvPatchVectorField
        (
            const massFluxInletVelocityFvPatchVectorField&,
            const vectorField&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone(const vectorField& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new massFluxInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the mass flux
            scalar massFlux() const
            {
                return massFlux_;
            }

            //- Return reference to the mass flux to allow adjustment
            scalar& massFlux()
            {
                return massFlux_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
