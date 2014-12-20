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

#include "CoupledPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class CoupledPointPatch,
    class Type
>
CoupledPointPatchField
<
    PatchField,
    PointPatch,
    CoupledPointPatch,
    Type
>::CoupledPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    PatchField<Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class CoupledPointPatch,
    class Type
>
CoupledPointPatchField
<
    PatchField,
    PointPatch,
    CoupledPointPatch,
    Type
>::CoupledPointPatchField
(
    const CoupledPointPatchField
        <PatchField, PointPatch, CoupledPointPatch, Type>&,
    const PointPatch& p,
    const Field<Type>& iF,
    const PointPatchFieldMapper&
)
:
    PatchField<Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class CoupledPointPatch,
    class Type
>
CoupledPointPatchField
<
    PatchField,
    PointPatch,
    CoupledPointPatch,
    Type
>::CoupledPointPatchField
(
    const CoupledPointPatchField
        <PatchField, PointPatch, CoupledPointPatch, Type>& ptf,
    const Field<Type>& iF
)
:
    PatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
