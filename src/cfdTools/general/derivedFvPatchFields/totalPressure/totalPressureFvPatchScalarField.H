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
    totalPressureFvPatchScalarField

Description

SourceFiles
    totalPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef totalPressureFvPatchScalarField_H
#define totalPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class totalPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class totalPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Total pressure
        scalarField p0_;


public:

    //- Runtime type information
    TypeName("totalPressure");


    // Constructors

        //- Construct from patch and internal field
        totalPressureFvPatchScalarField
        (
            const fvPatch&,
            const scalarField&
        );

        //- Construct from patch, internal field and dictionary
        totalPressureFvPatchScalarField
        (
            const fvPatch&,
            const scalarField&,
            const dictionary&
        );

        //- Construct by mapping given totalPressureFvPatchScalarField onto a new patch
        totalPressureFvPatchScalarField
        (
            const totalPressureFvPatchScalarField&,
            const fvPatch&,
            const scalarField&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new totalPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        totalPressureFvPatchScalarField
        (
            const totalPressureFvPatchScalarField&,
            const scalarField&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(const scalarField& iF) const
        {
            return tmp<fvPatchScalarField>
            (
                new totalPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the total pressure
            const scalarField& p0() const
            {
                return p0_;
            }

            //- Return reference to the total pressure to allow adjustment
            scalarField& p0()
            {
                return p0_;
            }


        // Evaluation functions

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
