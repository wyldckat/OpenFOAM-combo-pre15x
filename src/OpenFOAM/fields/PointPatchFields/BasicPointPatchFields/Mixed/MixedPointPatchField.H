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
    MixedPointPatchField. The boundary condition is a mix of a FixedValue and
    a ZeroGradient boundary condition.  I am still not sure how to do the
    FixedValue-fixedGradient combination. 

Description

SourceFiles
    MixedPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef MixedPointPatchField_H
#define MixedPointPatchField_H

#include "ValueStoredPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class MixedPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<template<class> class PatchField, class PointPatch, class Type>
class MixedPointPatchField
:
    public ValueStoredPointPatchField<PatchField, PointPatch, Type>
{
    // Private data

        //- Value field
        Field<Type> refValue_;

        //- Fraction (0-1) of value used for boundary condition
        scalarField valueFraction_;

    // Private member functions

        void checkFieldSize() const;


public:

    //- Runtime type information
    TypeName("mixed");


    // Constructors

        //- Construct from patch and internal field
        MixedPointPatchField
        (
            const PointPatch&,
            const Field<Type>&
        );

        //- Construct from patch and internal field patch field
        // and value fraction
        MixedPointPatchField
        (
            const PointPatch&,
            const Field<Type>&,
            const Field<Type>&,
            const Field<Type>&,
            const scalarField&
        );

        //- Construct from patch, internal field and dictionary
        MixedPointPatchField
        (
            const PointPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        MixedPointPatchField
        (
            const MixedPointPatchField<PatchField, PointPatch, Type>&,
            const PointPatch&,
            const Field<Type>&,
            const PointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new MixedPointPatchField<PatchField, PointPatch, Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        MixedPointPatchField
        (
            const MixedPointPatchField<PatchField, PointPatch, Type>&,
            const Field<Type>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<Type> > clone(const Field<Type>& iF) const
        {
            return autoPtr<PatchField<Type> >
            (
                new MixedPointPatchField<PatchField, PointPatch, Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Return defining fields

            virtual Field<Type>& refValue()
            {
                return refValue_;
            }

            virtual const Field<Type>& refValue() const
            {
                return refValue_;
            }

            virtual scalarField& valueFraction()
            {
                return valueFraction_;
            }

            virtual const scalarField& valueFraction() const
            {
                return valueFraction_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const PointPatchFieldMapper&
            );

            //- Reverse map the given PointPatchField onto
            //  this PointPatchField
            virtual void rmap
            (
                const PointPatchField<PatchField, PointPatch, Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the patch field
            virtual void updateBoundaryField();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "MixedPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
