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

#include "GlobalPointPatchField.H"
#include "lduMatrix.H"
#include "PstreamCombineReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class GlobalPointPatch,
    class Type
>
GlobalPointPatchField
<
    PatchField,
    PointPatch,
    GlobalPointPatch,
    Type
>
::GlobalPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        PointPatch,
        typename GlobalPointPatch::CoupledPointPatch,
        Type
    >(p, iF),
    globalPointPatch_(refCast<const GlobalPointPatch>(p))
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class GlobalPointPatch,
    class Type
>
GlobalPointPatchField
<
    PatchField,
    PointPatch,
    GlobalPointPatch,
    Type
>
::GlobalPointPatchField
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
        typename GlobalPointPatch::CoupledPointPatch,
        Type
    >(p, iF),
    globalPointPatch_(refCast<const GlobalPointPatch>(p))
{
    if (!isType<GlobalPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "GlobalPointPatchField<PatchField, PointPatch, "
            "GlobalPointPatch, Type>::"
            "GlobalPointPatchField\n"
            "(\n" " const PointPatch& p,\n"
            " const Field<Type>& field,\n"
            " const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index()
            << " not processorPoint type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class GlobalPointPatch,
    class Type
>
GlobalPointPatchField
<
    PatchField,
    PointPatch,
    GlobalPointPatch,
    Type
>
::GlobalPointPatchField
(
    const GlobalPointPatchField
    <
        PatchField,
        PointPatch,
        GlobalPointPatch,
        Type
    >& ptf,
    const PointPatch& p,
    const Field<Type>& iF,
    const PointPatchFieldMapper&
)
:
    CoupledPointPatchField
    <
        PatchField,
        PointPatch,
        typename GlobalPointPatch::CoupledPointPatch,
        Type
    >(p, iF),
    globalPointPatch_(refCast<const GlobalPointPatch>(ptf.patch()))
{
    if (!isType<GlobalPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "GlobalPointPatchField<PatchField, PointPatch, "
            "GlobalPointPatch, Type>::"
            "GlobalPointPatchField\n"
            "(\n"
            " const GlobalPointPatchField<PatchField, "
            "PointPatch, GlobalPointPatch, Type>& ptf,\n"
            " const PointPatch& p,\n"
            " const Field<Type>& iF,\n"
            " const PointPatchFieldMapper& mapper\n"
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
    class GlobalPointPatch,
    class Type
>
GlobalPointPatchField
<
    PatchField,
    PointPatch,
    GlobalPointPatch,
    Type
>
::GlobalPointPatchField
(
    const GlobalPointPatchField
    <
        PatchField,
        PointPatch,
        GlobalPointPatch,
        Type
    >& ptf,
    const Field<Type>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        PointPatch,
        typename GlobalPointPatch::CoupledPointPatch,
        Type
    >(ptf, iF),
    globalPointPatch_(refCast<const GlobalPointPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class GlobalPointPatch,
    class Type
>
GlobalPointPatchField
<
    PatchField,
    PointPatch,
    GlobalPointPatch,
    Type
>
::~GlobalPointPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class GlobalPointPatch,
    class Type
>
void GlobalPointPatchField
<
    PatchField,
    PointPatch,
    GlobalPointPatch,
    Type
>::swapAdd(Field<Type>& pField) const
{
    // Create the global list and insert local values
    if (globalPointPatch_.globalPointSize() > 0)
    {
        Field<Type> lpf = patchInternalField(pField);
        const labelList& addr = globalPointPatch_.sharedPointAddr();

        Field<Type> gpf
        (
            globalPointPatch_.globalPointSize(),
            pTraits<Type>::zero
        );

        forAll(addr, i)
        {
            gpf[addr[i]] += lpf[i];
        }

        combineReduce(gpf, plusEqOp<Field<Type> >());

        // Extract local data
        forAll (addr, i)
        {
            lpf[i] = gpf[addr[i]];
        }

        setInInternalField(pField, lpf);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
