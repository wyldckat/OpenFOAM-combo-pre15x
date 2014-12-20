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

#include "SlipTetPointPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
SlipTetPointPatchField<PatchField, PointPatch, Type>::SlipTetPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    BasicSymmetryTetPointPatchField<PatchField, PointPatch, Type>(p, iF)
{
    this->checkPointField();
}


template<template<class> class PatchField, class PointPatch, class Type>
SlipTetPointPatchField<PatchField, PointPatch, Type>::SlipTetPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const Field<Type>& f
)
:
    BasicSymmetryTetPointPatchField<PatchField, PointPatch, Type>(p, iF, f)
{
    this->checkPointField();
}


template<template<class> class PatchField, class PointPatch, class Type>
SlipTetPointPatchField<PatchField, PointPatch, Type>::SlipTetPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const dictionary&
)
:
    BasicSymmetryTetPointPatchField<PatchField, PointPatch, Type>(p, iF)
{
    this->checkPointField();
}


template<template<class> class PatchField, class PointPatch, class Type>
SlipTetPointPatchField<PatchField, PointPatch, Type>::SlipTetPointPatchField
(
    const SlipTetPointPatchField<PatchField, PointPatch, Type>&,
    const PointPatch& p,
    const Field<Type>& iF,
    const TetPointPatchFieldMapper&
)
:
    BasicSymmetryTetPointPatchField<PatchField, PointPatch, Type>(p, iF)
{
    this->checkPointField();
}


template<template<class> class PatchField, class PointPatch, class Type>
SlipTetPointPatchField<PatchField, PointPatch, Type>::SlipTetPointPatchField
(
    const SlipTetPointPatchField<PatchField, PointPatch, Type>& ptf,
    const Field<Type>& iF
)
:
    BasicSymmetryTetPointPatchField<PatchField, PointPatch, Type>(ptf, iF)
{
    this->checkPointField();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
