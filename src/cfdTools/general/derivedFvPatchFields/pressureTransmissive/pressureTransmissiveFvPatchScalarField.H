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
    pressureTransmissiveFvPatchScalarField

Description

SourceFiles
    pressureTransmissiveFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureTransmissiveFvPatchScalarFields_H
#define pressureTransmissiveFvPatchScalarFields_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class pressureTransmissiveFvPatch Declaration
\*---------------------------------------------------------------------------*/

class pressureTransmissiveFvPatchScalarField
:
    public mixedFvPatchScalarField
{

    // Private data

        //- Pressure of the free stream
        scalar pInf_;

        //- Relaxation length-scale
        scalar lInf_;

        //- Time index
        label curTimeIndex_;

        //- Old time patch field
        Field<scalar> pb0_;

        //- Old-old time patch field
        Field<scalar> pb00_;


public:

    //- Runtime type information
    TypeName("pressureTransmissive");


    // Constructors

        //- Construct from patch and internal field
        pressureTransmissiveFvPatchScalarField
        (
            const fvPatch&,
            const scalarField&
        );

        //- Construct from patch, internal field and dictionary
        pressureTransmissiveFvPatchScalarField
        (
            const fvPatch&,
            const scalarField&,
            const dictionary&
        );

        //- Construct by mapping given pressureTransmissiveFvPatchScalarField
        //  onto a new patch
        pressureTransmissiveFvPatchScalarField
        (
            const pressureTransmissiveFvPatchScalarField&,
            const fvPatch&,
            const scalarField&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new pressureTransmissiveFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureTransmissiveFvPatchScalarField
        (
            const pressureTransmissiveFvPatchScalarField&,
            const scalarField&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(const scalarField& iF) const
        {
            return tmp<fvPatchScalarField>
            (
                new pressureTransmissiveFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the pressure at infinity
            scalar pInf() const
            {
                return pInf_;
            }

            //- Return reference to the pressure at infinity to allow adjustment
            scalar& pInf()
            {
                return pInf_;
            }

            //- Return the relaxation length-scale
            scalar lInf() const
            {
                return lInf_;
            }

            //- Return reference to the relaxation length-scale
            //  to allow adjustment
            scalar& lInf()
            {
                return lInf_;
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
