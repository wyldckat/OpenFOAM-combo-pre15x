/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2004 OpenCFD Ltd.
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
    fixedTempRhoEFvPatchScalarField

Description

SourceFiles
    fixedTempRhoEFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedTempRhoEFvPatchScalarField_H
#define fixedTempRhoEFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fixedTempRhoEFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedTempRhoEFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        scalarField Twall_;

public:

    //- Runtime type information
    TypeName("fixedTempRhoE");


    // Constructors

        //- Construct from patch and internal field
        fixedTempRhoEFvPatchScalarField
        (
            const fvPatch&,
            const scalarField&
        );

        //- Construct from patch, internal field and dictionary
        fixedTempRhoEFvPatchScalarField
        (
            const fvPatch&,
            const scalarField&,
            const dictionary&
        );

        //- Construct by mapping given fixedTempRhoEFvPatchScalarField
        //  onto a new patch
        fixedTempRhoEFvPatchScalarField
        (
            const fixedTempRhoEFvPatchScalarField&,
            const fvPatch&,
            const scalarField&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedTempRhoEFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedTempRhoEFvPatchScalarField
        (
            const fixedTempRhoEFvPatchScalarField&,
            const scalarField&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(const scalarField& iF) const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedTempRhoEFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

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
