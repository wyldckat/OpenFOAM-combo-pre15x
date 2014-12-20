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

\*---------------------------------------------------------------------------*/

#include "FixedValueTetPointPatchField.H"
#include "boolList.H"
#include "Map.H"
#include "constraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
FixedValueTetPointPatchField<PatchField, PointPatch, Type>::
FixedValueTetPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    ValueStoredTetPointPatchField<PatchField, PointPatch, Type>(p, iF)
{}


template<template<class> class PatchField, class PointPatch, class Type>
FixedValueTetPointPatchField<PatchField, PointPatch, Type>::
FixedValueTetPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const Field<Type>& f
)
:
    ValueStoredTetPointPatchField<PatchField, PointPatch, Type>(p, iF, f)
{}


template<template<class> class PatchField, class PointPatch, class Type>
FixedValueTetPointPatchField<PatchField, PointPatch, Type>::
FixedValueTetPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    ValueStoredTetPointPatchField<PatchField, PointPatch, Type>(p, iF, dict)
{}


template<template<class> class PatchField, class PointPatch, class Type>
FixedValueTetPointPatchField<PatchField, PointPatch, Type>::
FixedValueTetPointPatchField
(
    const FixedValueTetPointPatchField<PatchField, PointPatch, Type>& ptf,
    const PointPatch& p,
    const Field<Type>& iF,
    const TetPointPatchFieldMapper& mapper
)
:
    ValueStoredTetPointPatchField<PatchField, PointPatch, Type>(ptf, p, iF, mapper)
{}


template<template<class> class PatchField, class PointPatch, class Type>
FixedValueTetPointPatchField<PatchField, PointPatch, Type>::
FixedValueTetPointPatchField
(
    const FixedValueTetPointPatchField<PatchField, PointPatch, Type>& ptf,
    const Field<Type>& iF
)
:
    ValueStoredTetPointPatchField<PatchField, PointPatch, Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set boundary condition to matrix
template<template<class> class PatchField, class PointPatch, class Type>
void FixedValueTetPointPatchField<PatchField, PointPatch, Type>::
setBoundaryCondition
(
    Map<constraint<Type> > & fix
) const
{
    const Field<Type>& values = *this;

    // get addressing
    const labelList& meshPoints = this->patch().meshPoints();

    forAll (meshPoints, pointI)
    {
        const label curPoint = meshPoints[pointI];

        // create a constraint
        constraint<Type> bc
        (
            curPoint,
            values[pointI],
            pTraits<Type>::one
        );

        // If not set add it, otherwise combine with
        // already existing value
        if (!fix.found(curPoint))
        {
            fix.insert(curPoint, bc);
        }
        else
        {
            fix[curPoint].combine(bc);
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Force an assignment
template<template<class> class PatchField, class PointPatch, class Type>
void FixedValueTetPointPatchField<PatchField, PointPatch, Type>::operator==
(
    const ValueStoredTetPointPatchField<PatchField, PointPatch, Type>& ptf
)
{
    Field<Type>::operator=(ptf);

    // insert the result into the internal field
    initEvaluate();
}


template<template<class> class PatchField, class PointPatch, class Type>
void FixedValueTetPointPatchField<PatchField, PointPatch, Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);

    // insert the result into the internal field
    initEvaluate();
}


template<template<class> class PatchField, class PointPatch, class Type>
void FixedValueTetPointPatchField<PatchField, PointPatch, Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);

    // insert the result into the internal field
    initEvaluate();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
