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

#include "WedgeTetPointPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class WedgePointPatch,
    class Type
>
WedgeTetPointPatchField<PatchField, PointPatch, WedgePointPatch, Type>::
WedgeTetPointPatchField
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
    class WedgePointPatch,
    class Type
>
WedgeTetPointPatchField<PatchField, PointPatch, WedgePointPatch, Type>::
WedgeTetPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const Field<Type>& f
)
:
    PatchField<Type>(p, iF, f)
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class WedgePointPatch,
    class Type
>
WedgeTetPointPatchField<PatchField, PointPatch, WedgePointPatch, Type>::
WedgeTetPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    PatchField<Type>(p, iF)
{
    if (!isType<WedgePointPatch>(p))
    {
        FatalIOErrorIn
        (
            "WedgeTetPointPatchField"
            "<PatchField, PointPatch, WedgePointPatch, Type>::"
            "WedgeTetPointPatchField\n"
            "(\n"
            "    const PointPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not wedge type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class WedgePointPatch,
    class Type
>
WedgeTetPointPatchField<PatchField, PointPatch, WedgePointPatch, Type>::
WedgeTetPointPatchField
(
    const WedgeTetPointPatchField<PatchField, PointPatch, WedgePointPatch, Type>&,
    const PointPatch& p,
    const Field<Type>& iF,
    const TetPointPatchFieldMapper&
)
:
    PatchField<Type>(p, iF)
{
    if (!isType<WedgePointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "WedgeTetPointPatchField"
            "<PatchField, PointPatch, WedgePointPatch, Type>::"
            "WedgeTetPointPatchField\n"
            "(\n"
            "    const WedgeTetPointPatchField"
                 "<PatchField, PointPatch, WedgePointPatch, Type>& ptf,\n"
            "    const PointPatch& p,\n"
            "    const Field<Type>& iF,\n"
            "    const TetPointPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class WedgePointPatch,
    class Type
>
WedgeTetPointPatchField<PatchField, PointPatch, WedgePointPatch, Type>::
WedgeTetPointPatchField
(
    const WedgeTetPointPatchField<PatchField, PointPatch, WedgePointPatch, Type>&
        ptf,
    const Field<Type>& iF
)
:
    PatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Evaluate patch field
template
<
    template<class> class PatchField,
    class PointPatch,
    class WedgePointPatch,
    class Type
>
void WedgeTetPointPatchField<PatchField, PointPatch, WedgePointPatch, Type>::
evaluate()
{
    // In order to ensure that the wedge patch is always flat, take the
    // normal vector from the first point
    const vector& nHat = this->patch().pointNormals()[0];

    tmp<Field<Type> > tvalues =
        transform(I - nHat*nHat, this->patchInternalField());
    const Field<Type>& values = tvalues();

    // Get internal field to insert values into
    Field<Type>& iF = const_cast<Field<Type>&>(this->internalField());

    // Get addressing
    const labelList& meshPoints = this->patch().meshPoints();

    forAll (meshPoints, pointI)
    {
        iF[meshPoints[pointI]] = values[pointI];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
