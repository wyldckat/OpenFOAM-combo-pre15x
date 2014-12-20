/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "ValuePointPatchField.H"
#include "PointPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::
checkFieldSize() const
{
    if (size() != this->patch().size())
    {
        FatalErrorIn
        (
            "void ValuePointPatchField<PatchField, PointPatch, Type>::"
            "checkField() const"
        )   << "field does not correspond to patch. " << endl
            << "Field size: " << size() << " patch size: "
            << this->patch().size()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
ValuePointPatchField<PatchField, PointPatch, Type>::
ValuePointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    PatchField<Type>(p, iF),
    Field<Type>(p.size())
{}


template<template<class> class PatchField, class PointPatch, class Type>
ValuePointPatchField<PatchField, PointPatch, Type>::
ValuePointPatchField
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
ValuePointPatchField<PatchField, PointPatch, Type>::
ValuePointPatchField
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
ValuePointPatchField<PatchField, PointPatch, Type>::
ValuePointPatchField
(
    const ValuePointPatchField<PatchField, PointPatch, Type>& ptf,
    const PointPatch& p,
    const Field<Type>& iF,
    const PointPatchFieldMapper& mapper
)
:
    PatchField<Type>(p, iF),
    Field<Type>(ptf, mapper)
{}


template<template<class> class PatchField, class PointPatch, class Type>
ValuePointPatchField<PatchField, PointPatch, Type>::
ValuePointPatchField
(
    const ValuePointPatchField<PatchField, PointPatch, Type>& ptf,
    const Field<Type>& iF
)
:
    PatchField<Type>(ptf, iF),
    Field<Type>(ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::autoMap
(
    const PointPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
}


// Grab the values using rmap
template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::rmap
(
    const PointPatchField<PatchField, PointPatch, Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap
    (
        refCast<const ValuePointPatchField<PatchField, PointPatch, Type> >
        (
            ptf
        ),
        addr
    );
}


// Insert boundary value into the internal field
template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::evaluate()
{
    // Get internal field to insert values into
    Field<Type>& iF = const_cast<Field<Type>&>(this->internalField());

    setInInternalField(iF, *this);
}


// Write
template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::
write(Ostream& os) const
{
    PatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::operator=
(
    const ValuePointPatchField<PatchField, PointPatch, Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::operator=
(
    const PointPatchField<PatchField, PointPatch, Type>& ptf
)
{
    Field<Type>::operator=(ptf.patchInternalField());
}


template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::operator=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// Force an assignment
template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::operator==
(
    const ValuePointPatchField<PatchField, PointPatch, Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::operator==
(
    const PointPatchField<PatchField, PointPatch, Type>& ptf
)
{
    Field<Type>::operator=(ptf.patchInternalField());
}


template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<template<class> class PatchField, class PointPatch, class Type>
void ValuePointPatchField<PatchField, PointPatch, Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
