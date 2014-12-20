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

#include "MixedPointPatchField.H"
#include "PointPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
void MixedPointPatchField<PatchField, PointPatch, Type>::checkFieldSize() const
{
    if
    (
        this->size() != this->patch().size()
     || refValue_.size() != this->patch().size()
     || valueFraction_.size() != this->patch().size()
    )
    {
        FatalErrorIn
        (
            "void MixedPointPatchField<PatchField, PointPatch, Type>::"
            "checkField() const"
        )   << "field does not correspond to patch. " << endl
            << "Field size: " << this->size() << " value size: "
            << refValue_.size()
            << " valueFraction size: " << valueFraction_.size()
            << " patch size: " << this->patch().size()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
MixedPointPatchField<PatchField, PointPatch, Type>::MixedPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    ValuePointPatchField<PatchField, PointPatch, Type>(p, iF),
    refValue_(p.size()),
    valueFraction_(p.size())
{}


template<template<class> class PatchField, class PointPatch, class Type>
MixedPointPatchField<PatchField, PointPatch, Type>::MixedPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const Field<Type>& f,
    const Field<Type>& v,
    const scalarField& vf
)
:
    ValuePointPatchField<PatchField, PointPatch, Type>(p, iF, f),
    refValue_(v),
    valueFraction_(vf)
{
    checkFieldSize();
}


template<template<class> class PatchField, class PointPatch, class Type>
MixedPointPatchField<PatchField, PointPatch, Type>::MixedPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    ValuePointPatchField<PatchField, PointPatch, Type>(p, iF),
    refValue_("refValue", dict, p.size()),
    valueFraction_("valueFraction", dict, p.size())
{}


template<template<class> class PatchField, class PointPatch, class Type>
MixedPointPatchField<PatchField, PointPatch, Type>::MixedPointPatchField
(
    const MixedPointPatchField<PatchField, PointPatch, Type>& ptf,
    const PointPatch& p,
    const Field<Type>& iF,
    const PointPatchFieldMapper& mapper
)
:
    ValuePointPatchField<PatchField, PointPatch, Type>
    (
        ptf,
        p,
        iF,
        mapper
    ),
    refValue_(ptf.refValue_, mapper),
    valueFraction_(ptf.valueFraction_, mapper)

{}


template<template<class> class PatchField, class PointPatch, class Type>
MixedPointPatchField<PatchField, PointPatch, Type>::MixedPointPatchField
(
    const MixedPointPatchField<PatchField, PointPatch, Type>& ptf,
    const Field<Type>& iF
)
:
    ValuePointPatchField<PatchField, PointPatch, Type>(ptf, iF),
    refValue_(ptf.refValue_),
    valueFraction_(ptf.valueFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map and resize from self given a mapper
template<template<class> class PatchField, class PointPatch, class Type>
void MixedPointPatchField<PatchField, PointPatch, Type>::autoMap
(
    const PointPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
    refValue_.autoMap(m);
    valueFraction_.autoMap(m);
}


// Grab the values using rmap
template<template<class> class PatchField, class PointPatch, class Type>
void MixedPointPatchField<PatchField, PointPatch, Type>::rmap
(
    const PointPatchField<PatchField, PointPatch, Type>& ptf,
    const labelList& addr
)
{
    const MixedPointPatchField<PatchField, PointPatch, Type>& mptf =
        refCast<const MixedPointPatchField<PatchField, PointPatch, Type> >(ptf);

    Field<Type>::rmap(mptf, addr);
    refValue_.rmap(mptf.refValue_, addr);
    valueFraction_.rmap(mptf.valueFraction_, addr);
}


// Evaluate patch field
template<template<class> class PatchField, class PointPatch, class Type>
void MixedPointPatchField<PatchField, PointPatch, Type>::evaluate()
{
    Field<Type>::operator=
    (
        valueFraction_*refValue_
      + (1.0 - valueFraction_)*this->patchInternalField()
    );

    // Get internal field to insert values into
    Field<Type>& iF = const_cast<Field<Type>&>(this->internalField());

    setInInternalField(iF, *this);
}


// Write
template<template<class> class PatchField, class PointPatch, class Type>
void
MixedPointPatchField<PatchField, PointPatch, Type>::write(Ostream& os) const
{
    PatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    valueFraction_.writeEntry("valueFraction", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
