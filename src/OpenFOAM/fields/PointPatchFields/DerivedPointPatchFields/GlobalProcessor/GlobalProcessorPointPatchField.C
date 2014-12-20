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

#include "GlobalProcessorPointPatchField.H"
#include "lduMatrix.H"
#include "PstreamCombineReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Reduce the field and extract the local values
template
<
    template<class> class PatchField,
    class PointPatch,
    class GlobalProcessorPointPatch,
    class Type
>
template<class Type2>
tmp<Field<Type2> > GlobalProcessorPointPatchField
<
    PatchField,
    PointPatch,
    GlobalProcessorPointPatch,
    Type
>
::reduceExtractPoint
(
    const tmp<Field<Type2> >& tpField
) const
{
    // Create the global list and insert local values
    if (procPatch_.globalPointSize() > 0)
    {
        Field<Type2> gpf(procPatch_.globalPointSize(), pTraits<Type2>::zero);

        const labelList& addr = procPatch_.sharedPointAddr();
        const Field<Type2>& pField = tpField();

        forAll (addr, i)
        {
            gpf[addr[i]] = pField[i];
        }

        combineReduce(gpf, plusEqOp<Field<Type2> >());

        // Extract local data
        tmp<Field<Type2> > tlpf(new Field<Type2>(addr.size()));
        Field<Type2>& lpf = tlpf();

        forAll (addr, i)
        {
            lpf[i] = gpf[addr[i]];
        }

        return tlpf;
    }
    else
    {
        return tpField;
    }
}


// Reduce the field and extract the local values
template
<
    template<class> class PatchField,
    class PointPatch,
    class GlobalProcessorPointPatch,
    class Type
>
template<class Type2>
tmp<Field<Type2> > GlobalProcessorPointPatchField
<
    PatchField,
    PointPatch,
    GlobalProcessorPointPatch,
    Type
>
::reduceExtractEdge
(
    const tmp<Field<Type2> >& teField
) const
{
    if (procPatch_.globalEdgeSize() > 0)
    {
        // Create the global list and insert local values
        Field<Type2> gef(procPatch_.globalEdgeSize(), pTraits<Type2>::zero);

        const labelList& addr = procPatch_.sharedEdgeAddr();
        const Field<Type2>& eField = teField();

        forAll (addr, i)
        {
            gef[addr[i]] = eField[i];
        }

        combineReduce(gef, plusEqOp<Field<Type2> >());

        // Extract local data
        tmp<Field<Type2> > tlef(new Field<Type2>(addr.size()));
        Field<Type2>& lef = tlef();

        forAll (addr, i)
        {
            lef[i] = gef[addr[i]];
        }

        return tlef;
    }
    else
    {
        return teField;
    }
}


// Add the diagonal/source to the internal field.
template
<
    template<class> class PatchField,
    class PointPatch,
    class GlobalProcessorPointPatch,
    class Type
>
template<class Type2>
void GlobalProcessorPointPatchField
<
    PatchField,
    PointPatch,
    GlobalProcessorPointPatch,
    Type
>
::addFieldTempl
(
    Field<Type2>& pField
) const
{
    // Set the values from the global sum
    tmp<Field<Type2> > trpf =
        reduceExtractPoint<Type2>(patchInternalField(pField));
    Field<Type2>& rpf = trpf();

    // Get addressing
    const labelList& addr = procPatch_.meshPoints();

    forAll (addr, i)
    {
        pField[addr[i]] = rpf[i];
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class GlobalProcessorPointPatch,
    class Type
>
GlobalProcessorPointPatchField
<
    PatchField,
    PointPatch,
    GlobalProcessorPointPatch,
    Type
>
::GlobalProcessorPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    CoupledPointPatchField<PatchField, PointPatch, Type>(p, iF),
    procPatch_(refCast<const GlobalProcessorPointPatch>(p))
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class GlobalProcessorPointPatch,
    class Type
>
GlobalProcessorPointPatchField
<
    PatchField,
    PointPatch,
    GlobalProcessorPointPatch,
    Type
>
::GlobalProcessorPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    CoupledPointPatchField<PatchField, PointPatch, Type>(p, iF),
    procPatch_(refCast<const GlobalProcessorPointPatch>(p))
{
    if (!isType<GlobalProcessorPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "GlobalProcessorPointPatchField<PatchField, PointPatch, "
            "GlobalProcessorPointPatch, Type>::"
            "GlobalProcessorPointPatchField\n"
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
    class GlobalProcessorPointPatch,
    class Type
>
GlobalProcessorPointPatchField
<
    PatchField,
    PointPatch,
    GlobalProcessorPointPatch,
    Type
>
::GlobalProcessorPointPatchField
(
    const GlobalProcessorPointPatchField
    <
        PatchField,
        PointPatch,
        GlobalProcessorPointPatch,
        Type
    >& ptf,
    const PointPatch& p,
    const Field<Type>& iF,
    const PointPatchFieldMapper&
)
:
    CoupledPointPatchField<PatchField, PointPatch, Type>(p, iF),
    procPatch_(refCast<const GlobalProcessorPointPatch>(ptf.patch()))
{
    if (!isType<GlobalProcessorPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "GlobalProcessorPointPatchField<PatchField, PointPatch, "
            "GlobalProcessorPointPatch, Type>::"
            "GlobalProcessorPointPatchField\n"
            "(\n"
            " const GlobalProcessorPointPatchField<PatchField, "
            "PointPatch, GlobalProcessorPointPatch, Type>& ptf,\n"
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
    class GlobalProcessorPointPatch,
    class Type
>
GlobalProcessorPointPatchField
<
    PatchField,
    PointPatch,
    GlobalProcessorPointPatch,
    Type
>
::GlobalProcessorPointPatchField
(
    const GlobalProcessorPointPatchField
    <
        PatchField,
        PointPatch,
        GlobalProcessorPointPatch,
        Type
    >& ptf,
    const Field<Type>& iF
)
:
    CoupledPointPatchField<PatchField, PointPatch, Type>(ptf, iF),
    procPatch_(refCast<const GlobalProcessorPointPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class GlobalProcessorPointPatch,
    class Type
>
GlobalProcessorPointPatchField
<
    PatchField,
    PointPatch,
    GlobalProcessorPointPatch,
    Type
>
::~GlobalProcessorPointPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class GlobalProcessorPointPatch,
    class Type
>
void GlobalProcessorPointPatchField
<
    PatchField,
    PointPatch,
    GlobalProcessorPointPatch,
    Type
>
::addField(Field<Type>& d) const
{
    addFieldTempl(d);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
