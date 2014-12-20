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

#include "FixedValuePointPatchField.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
FixedValuePointPatchField<PatchField, PointPatch, Type>::
FixedValuePointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    ValuePointPatchField<PatchField, PointPatch, Type>(p, iF)
{}


template<template<class> class PatchField, class PointPatch, class Type>
FixedValuePointPatchField<PatchField, PointPatch, Type>::
FixedValuePointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const Field<Type>& f
)
:
    ValuePointPatchField<PatchField, PointPatch, Type>(p, iF, f)
{}


template<template<class> class PatchField, class PointPatch, class Type>
FixedValuePointPatchField<PatchField, PointPatch, Type>::
FixedValuePointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    ValuePointPatchField<PatchField, PointPatch, Type>(p, iF, dict)
{}


template<template<class> class PatchField, class PointPatch, class Type>
FixedValuePointPatchField<PatchField, PointPatch, Type>::
FixedValuePointPatchField
(
    const FixedValuePointPatchField<PatchField, PointPatch, Type>& ptf,
    const PointPatch& p,
    const Field<Type>& iF,
    const PointPatchFieldMapper& mapper
)
:
    ValuePointPatchField<PatchField, PointPatch, Type>(ptf, p, iF, mapper)
{}


template<template<class> class PatchField, class PointPatch, class Type>
FixedValuePointPatchField<PatchField, PointPatch, Type>::
FixedValuePointPatchField
(
    const FixedValuePointPatchField<PatchField, PointPatch, Type>& ptf,
    const Field<Type>& iF
)
:
    ValuePointPatchField<PatchField, PointPatch, Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
