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
    adiabaticRhoEFvPatchScalarField

Description

SourceFiles
    adiabaticRhoEFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef adiabaticRhoEFvPatchScalarFields_H
#define adiabaticRhoEFvPatchScalarFields_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class adiabaticRhoEFvPatch Declaration
\*---------------------------------------------------------------------------*/

class adiabaticRhoEFvPatchScalarField
:
    public mixedFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("adiabaticRhoE");


    // Constructors

        //- Construct from patch and internal field
        adiabaticRhoEFvPatchScalarField
        (
            const fvPatch&,
            const scalarField&
        );

        //- Construct from patch, internal field and dictionary
        adiabaticRhoEFvPatchScalarField
        (
            const fvPatch&,
            const scalarField&,
            const dictionary&
        );

        //- Construct by mapping given adiabaticRhoEFvPatchScalarField
        //  onto a new patch
        adiabaticRhoEFvPatchScalarField
        (
            const adiabaticRhoEFvPatchScalarField&,
            const fvPatch&,
            const scalarField&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new adiabaticRhoEFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        adiabaticRhoEFvPatchScalarField
        (
            const adiabaticRhoEFvPatchScalarField&,
            const scalarField&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(const scalarField& iF) const
        {
            return tmp<fvPatchScalarField>
            (
                new adiabaticRhoEFvPatchScalarField(*this, iF)
            );
        }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );


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
