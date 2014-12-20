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
    ValueStoredPointPatchField

Description

SourceFiles
    ValueStoredPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef ValueStoredPointPatchField_H
#define ValueStoredPointPatchField_H

#include "PointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class ValueStoredPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<template<class> class PatchField, class PointPatch, class Type>
class ValueStoredPointPatchField
:
    public PatchField<Type>,
    public Field<Type>
{
    // Private member functions

        void checkFieldSize() const;


public:

    //- Runtime type information
    TypeName("valueStored");


    // Constructors

        //- Construct from patch and internal field
        ValueStoredPointPatchField
        (
            const PointPatch&,
            const Field<Type>&
        );

        //- Construct from patch and internal field and patch field
        ValueStoredPointPatchField
        (
            const PointPatch&,
            const Field<Type>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        ValueStoredPointPatchField
        (
            const PointPatch&,
            const Field<Type>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        ValueStoredPointPatchField
        (
            const ValueStoredPointPatchField<PatchField, PointPatch, Type>&,
            const PointPatch&,
            const Field<Type>&,
            const PointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<PatchField<Type> > clone() const
        {
            return autoPtr<PatchField<Type> >
            (
                new ValueStoredPointPatchField<PatchField, PointPatch, Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        ValueStoredPointPatchField
        (
            const ValueStoredPointPatchField<PatchField, PointPatch, Type>&,
            const Field<Type>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<PatchField<Type> > clone(const Field<Type>& iF) const
        {
            return autoPtr<PatchField<Type> >
            (
                new ValueStoredPointPatchField<PatchField, PointPatch, Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return size
            label size() const
            {
                return Field<Type>::size();
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

            //- Does it store field data
            virtual bool storesFieldData() const
            {
                return true;
            }

            //- Insert boundary value into the internal field
            virtual void initEvaluate(const bool=false)
            {}

            virtual void evaluate();

            //- Update boundary value
            virtual void updateBoundaryField()
            {}

            //- Set boundary condition to matrix
            virtual void setBoundaryCondition
            (
                Map<constraint<Type> > &
            ) const
            {}

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        //- Assignment operators
        virtual void operator=
        (
            const ValueStoredPointPatchField<PatchField, PointPatch, Type>&
        );

        virtual void operator=(const Field<Type>&);
        virtual void operator=(const Type&);

        // Force an assignment irrespective of form of patch

        virtual void operator==
        (
            const ValueStoredPointPatchField<PatchField, PointPatch, Type>&
        );

        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);

};


// This function is added to overide the hack in PointPatchField.H which enables
// simple backward compatability with verions using referenceLevel
// in GeometicField
template<template<class> class PatchField, class PointPatch, class Type>
tmp<Field<Type> > operator+
(
    const ValueStoredPointPatchField<PatchField, PointPatch, Type>& vsppf,
    const Type& t
)
{
    return static_cast<const Field<Type>&>(vsppf) + t;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "ValueStoredPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
