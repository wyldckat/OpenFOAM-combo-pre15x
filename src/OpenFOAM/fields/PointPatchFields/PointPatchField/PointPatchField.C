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

#include "PointPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
void PointPatchField<PatchField, PointPatch, Type>::checkInternalField() const
{
    const label iFs = internalField_.size();

    const label nPoints = patch_.boundaryMesh().mesh().nPoints();
    const label nCells = patch_.boundaryMesh().mesh().nCells();

    if (iFs != nPoints && iFs != nCells)
    {
        FatalErrorIn
        (
            "void PointPatchField<PatchField, PointPatch, Type>::"
            "checkInternalField() const"
        )   << "internal field is not point or element field." << nl
            << "Size of internal field: " << iFs << ".  Number of points: "
            << nPoints << ".  Number of cells: " << nCells << "."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
PointPatchField<PatchField, PointPatch, Type>::PointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    patch_(p),
    internalField_(iF)
{
    checkInternalField();
}


template<template<class> class PatchField, class PointPatch, class Type>
PointPatchField<PatchField, PointPatch, Type>::PointPatchField
(
    const PointPatchField<PatchField, PointPatch, Type>& ptf
)
:
    lduCoupledInterface(),
    patch_(ptf.patch_),
    internalField_(ptf.internalField_)
{}


template<template<class> class PatchField, class PointPatch, class Type>
PointPatchField<PatchField, PointPatch, Type>::PointPatchField
(
    const PointPatchField<PatchField, PointPatch, Type>& ptf,
    const Field<Type>& iF
)
:
    patch_(ptf.patch_),
    internalField_(iF)
{
    checkInternalField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return objectRegistry
template<template<class> class PatchField, class PointPatch, class Type>
const objectRegistry& PointPatchField<PatchField, PointPatch, Type>::db() const
{
    return patch_.boundaryMesh().mesh();
}


// Write
template<template<class> class PatchField, class PointPatch, class Type>
void PointPatchField<PatchField, PointPatch, Type>::write(Ostream& os) const
{
    os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
}


// Return field created from appropriate internal field values
template<template<class> class PatchField, class PointPatch, class Type>
tmp<Field<Type> > PointPatchField<PatchField, PointPatch, Type>::
patchInternalField() const
{
    return patchInternalField(internalField());
}


// Return field created from appropriate internal field values
// given the internal field
template<template<class> class PatchField, class PointPatch, class Type>
template<class Type1>
tmp<Field<Type1> > PointPatchField<PatchField, PointPatch, Type>::
patchInternalField
(
    const Field<Type1>& iF
) const
{
    // Check size
    if (iF.size() != internalField().size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type1> > PointPatchField<PatchField, PointPatch, "
            "Type>::"
            "patchInternalField(const Field<Type1>& iF) const"
        )   << "given internal field does not correspond to the mesh. "
            << "Field size: " << iF.size()
            << " mesh size: " << internalField().size()
            << abort(FatalError);
    }

    // get addressing
    const labelList& meshPoints = patch().meshPoints();

    tmp<Field<Type1> > tvalues(new Field<Type1>(meshPoints.size()));
    Field<Type1>& values = tvalues();

    forAll (meshPoints, pointI)
    {
        values[pointI] = iF[meshPoints[pointI]];
    }

    return tvalues;
}


// Does this patchField correspond to a pointTypeField
template<template<class> class PatchField, class PointPatch, class Type>
bool PointPatchField<PatchField, PointPatch, Type>::isPointField() const
{
    return
        internalField().size()
     == patch().boundaryMesh().mesh().nPoints();
}

// Does this patchField correspond to a pointTypeField
template<template<class> class PatchField, class PointPatch, class Type>
void PointPatchField<PatchField, PointPatch, Type>::checkPointField() const
{
    if (!isPointField())
    {
        FatalErrorIn
        (
            "void PointPatchField<PatchField, PointPatch, Type>::"
            "::checkPointField() const"
        )   << "This " << type() << " patchField"
            << " is not part of a pointTypeField which may cause "
            << "undefined behaviour from the evaluate and other functions"
            << abort(FatalError);
    }
}


template<template<class> class PatchField, class PointPatch, class Type>
template<class Type1>
void PointPatchField<PatchField, PointPatch, Type>::addToInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF
) const
{
    // Check size
    if (iF.size() != internalField().size())
    {
        FatalErrorIn
        (
            "void PointPatchField<PatchField, PointPatch, Type>::"
            "addToInternalField("
            "Field<Type1>& iF, const Field<Type1>& iF) const"
        )   << "given internal field does not correspond to the mesh. "
            << "Field size: " << iF.size()
            << " mesh size: " << internalField().size()
            << abort(FatalError);
    }

    if (pF.size() != size())
    {
        FatalErrorIn
        (
            "void PointPatchField<PatchField, PointPatch, Type>::"
            "addToInternalField("
            "Field<Type1>& iF, const Field<Type1>& iF) const"
        )   << "given patch field does not correspond to the mesh. "
            << "Field size: " << pF.size()
            << " mesh size: " << size()
            << abort(FatalError);
    }

    // Get the addressing
    const labelList& mp = patch().meshPoints();

    forAll (mp, pointI)
    {
        iF[mp[pointI]] += pF[pointI];
    }
}


template<template<class> class PatchField, class PointPatch, class Type>
template<class Type1>
void PointPatchField<PatchField, PointPatch, Type>::setInInternalField
(
    Field<Type1>& iF,
    const Field<Type1>& pF
) const
{
    // Check size
    if (iF.size() != internalField().size())
    {
        FatalErrorIn
        (
            "void PointPatchField<PatchField, PointPatch, Type>::"
            "setInInternalField("
            "Field<Type1>& iF, const Field<Type1>& iF) const"
        )   << "given internal field does not correspond to the mesh. "
            << "Field size: " << iF.size()
            << " mesh size: " << internalField().size()
            << abort(FatalError);
    }

    if (pF.size() != size())
    {
        FatalErrorIn
        (
            "void PointPatchField<PatchField, PointPatch, Type>::"
            "setInInternalField("
            "Field<Type1>& iF, const Field<Type1>& iF) const"
        )   << "given patch field does not correspond to the mesh. "
            << "Field size: " << pF.size()
            << " mesh size: " << size()
            << abort(FatalError);
    }

    // Get the addressing
    const labelList& mp = patch().meshPoints();

    forAll (mp, pointI)
    {
        iF[mp[pointI]] = pF[pointI];
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<template<class> class PatchField, class PointPatch, class Type>
Ostream& operator<<
(
    Ostream& os,
    const PointPatchField<PatchField, PointPatch, Type>& ptf
)
{
    ptf.write(os);

    os.check
    (
            "Ostream& operator<<"
            "(Ostream&, const PointPatchField<PatchField, PointPatch, Type>&"
    );

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#    include "newPointPatchField.C"

// ************************************************************************* //
