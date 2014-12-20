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

#include "CyclicPointPatchField.H"
#include "Swap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class CyclicPointPatch,
    class Type
>
CyclicPointPatchField<PatchField, PointPatch, CyclicPointPatch, Type>::
CyclicPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        PointPatch,
        typename CyclicPointPatch::CoupledPointPatch,
        Type
    >(p, iF),
    cyclicPatch_(refCast<const CyclicPointPatch>(p))
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class CyclicPointPatch,
    class Type
>
CyclicPointPatchField<PatchField, PointPatch, CyclicPointPatch, Type>::
CyclicPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const Field<Type>& f
)
:
    CoupledPointPatchField
    <
        PatchField,
        PointPatch,
        typename CyclicPointPatch::CoupledPointPatch,
        Type
    >(p, iF, f),
    cyclicPatch_(refCast<const CyclicPointPatch>(p))
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class CyclicPointPatch,
    class Type
>
CyclicPointPatchField<PatchField, PointPatch, CyclicPointPatch, Type>::
CyclicPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    CoupledPointPatchField
    <
        PatchField,
        PointPatch,
        typename CyclicPointPatch::CoupledPointPatch,
        Type
    >(p, iF),
    cyclicPatch_(refCast<const CyclicPointPatch>(p))
{
    if (!isType<CyclicPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "CyclicPointPatchField"
            "<PatchField, PointPatch, CyclicPointPatch, Type>::"
            "CyclicPointPatchField\n"
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
    class CyclicPointPatch,
    class Type
>
CyclicPointPatchField<PatchField, PointPatch, CyclicPointPatch, Type>::
CyclicPointPatchField
(
    const CyclicPointPatchField
        <PatchField, PointPatch, CyclicPointPatch, Type>&,
    const PointPatch& p,
    const Field<Type>& iF,
    const PointPatchFieldMapper&
)
:
    CoupledPointPatchField
    <
        PatchField,
        PointPatch,
        typename CyclicPointPatch::CoupledPointPatch,
        Type
    >(p, iF),
    cyclicPatch_(refCast<const CyclicPointPatch>(p))
{
    if (!isType<CyclicPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "CyclicPointPatchField"
            "<PatchField, PointPatch, CyclicPointPatch, Type>::"
            "CyclicPointPatchField\n"
            "(\n"
            "    const CyclicPointPatchField"
                 "<PatchField, PointPatch, CyclicPointPatch, Type>& ptf,\n"
            "    const PointPatch& p,\n"
            "    const Field<Type>& iF,\n"
            "    const PointPatchFieldMapper& mapper\n"
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
    class CyclicPointPatch,
    class Type
>
CyclicPointPatchField<PatchField, PointPatch, CyclicPointPatch, Type>::
CyclicPointPatchField
(
    const CyclicPointPatchField
        <PatchField, PointPatch, CyclicPointPatch, Type>& ptf,
    const Field<Type>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        PointPatch,
        typename CyclicPointPatch::CoupledPointPatch,
        Type
    >(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class CyclicPointPatch,
    class Type
>
void CyclicPointPatchField<PatchField, PointPatch, CyclicPointPatch, Type>::
swapAdd(Field<Type>& pField) const
{
    Field<Type> pf(this->patchInternalField(pField));

    const edgeList& pairs = cyclicPatch_.transformPairs();

    if (doTransform())
    {
        forAll(pairs, pairi)
        {
            Type tmp = pf[pairs[pairi][0]];
            pf[pairs[pairi][0]] = transform(forwardT()[0], pf[pairs[pairi][1]]);
            pf[pairs[pairi][1]] = transform(reverseT()[0], tmp);
        }
    }
    else
    {
        forAll(pairs, pairi)
        {
            Swap(pf[pairs[pairi][0]], pf[pairs[pairi][1]]);
        }
    }

    addToInternalField(pField, pf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
