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

Description

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "ValueStoredPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
void ValueStoredPointPatchField<PatchField, PointPatch, Type>::
checkFieldSize() const
{
    if (size() != this->patchMesh().size())
    {
        FatalErrorIn
        (
            "void ValueStoredPointPatchField<PatchField, PointPatch, Type>::"
            "checkField() const"
        )   << "field does not correspond to patch. " << endl
            << "Field size: " << size() << " patch size: "
            << this->patchMesh().size()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
ValueStoredPointPatchField<PatchField, PointPatch, Type>::
ValueStoredPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    PatchField<Type>(p, iF),
    Field<Type>(p.size())
{}


template<template<class> class PatchField, class PointPatch, class Type>
ValueStoredPointPatchField<PatchField, PointPatch, Type>::
ValueStoredPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const Field<Type>& f
)
:
    PatchField<Type>(p, iF),
    Field<Type>(f)
{
    checkFieldSize();
}


template<template<class> class PatchField, class PointPatch, class Type>
ValueStoredPointPatchField<PatchField, PointPatch, Type>::
ValueStoredPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    PatchField<Type>(p, iF),
    Field<Type>("value", dict, p.size())
{}


template<template<class> class PatchField, class PointPatch, class Type>
ValueStoredPointPatchField<PatchField, PointPatch, Type>::
ValueStoredPointPatchField
(
    const ValueStoredPointPatchField<PatchField, PointPatch, Type>& ptf,
    const PointPatch& p,
    const Field<Type>& iF,
    const PointPatchFieldMapper& mapper
)
:
    PatchField<Type>(p, iF),
    Field<Type>((const Field<Type>&)ptf, (const FieldMapper&)mapper)
{}


template<template<class> class PatchField, class PointPatch, class Type>
ValueStoredPointPatchField<PatchField, PointPatch, Type>::
ValueStoredPointPatchField
(
    const ValueStoredPointPatchField<PatchField, PointPatch, Type>& ptf,
    const Field<Type>& iF
)
:
    PatchField<Type>(ptf, iF),
    Field<Type>(ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
template<template<class> class PatchField, class PointPatch, class Type>
void ValueStoredPointPatchField<PatchField, PointPatch, Type>::autoMap
(
    const PointPatchFieldMapper& m
)
{
    Field<Type>::autoMap((const FieldMapper&)(m));
}


// Grab the values using rmap
template<template<class> class PatchField, class PointPatch, class Type>
void ValueStoredPointPatchField<PatchField, PointPatch, Type>::rmap
(
    const PatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap
    (
        refCast<const ValueStoredPointPatchField<PatchField, PointPatch, Type> >
        (
            ptf
        ),
        addr
    );
}


// Insert boundary value into the internal field
template<template<class> class PatchField, class PointPatch, class Type>
void ValueStoredPointPatchField<PatchField, PointPatch, Type>::evaluate()
{
    // Evaluate boundary condition
    updateBoundaryField();

    if (this->isPointField())
    {
        // Get value
        const Field<Type>& values = *this;

        // Get internal field to insert values into
        Field<Type>& iF = ((Field<Type>&)(this->internalField()));

        setInInternalField(iF, values);
    }
}


// Write
template<template<class> class PatchField, class PointPatch, class Type>
void ValueStoredPointPatchField<PatchField, PointPatch, Type>::
write(Ostream& os) const
{
    PatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
void ValueStoredPointPatchField<PatchField, PointPatch, Type>::operator=
(
    const ValueStoredPointPatchField<PatchField, PointPatch, Type>& ptf
)
{
    (Field<Type>&)(*this) = (Field<Type>&)ptf;
}


template<template<class> class PatchField, class PointPatch, class Type>
void ValueStoredPointPatchField<PatchField, PointPatch, Type>::operator=
(
    const Field<Type>& tf
)
{
    (Field<Type>&)(*this) = tf;
}


template<template<class> class PatchField, class PointPatch, class Type>
void ValueStoredPointPatchField<PatchField, PointPatch, Type>::operator=
(
    const Type& t
)
{
    (Field<Type>&)(*this) = t;
}


// Force an assignment
template<template<class> class PatchField, class PointPatch, class Type>
void ValueStoredPointPatchField<PatchField, PointPatch, Type>::operator==
(
    const ValueStoredPointPatchField<PatchField, PointPatch, Type>& ptf
)
{
    (Field<Type>&)(*this) = (Field<Type>&)ptf;
}


template<template<class> class PatchField, class PointPatch, class Type>
void ValueStoredPointPatchField<PatchField, PointPatch, Type>::operator==
(
    const Field<Type>& tf
)
{
    (Field<Type>&)(*this) = tf;
}


template<template<class> class PatchField, class PointPatch, class Type>
void ValueStoredPointPatchField<PatchField, PointPatch, Type>::operator==
(
    const Type& t
)
{
    (Field<Type>&)(*this) = t;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
