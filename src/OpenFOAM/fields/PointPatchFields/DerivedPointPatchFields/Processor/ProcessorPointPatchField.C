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

#include "ProcessorPointPatchField.H"
#include "transformField.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
ProcessorPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        PointPatch,
        typename ProcessorPointPatch::CoupledPointPatch,
        Type
    >(p, iF),
    procPatch_(refCast<const ProcessorPointPatch>(p))
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
ProcessorPointPatchField
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
        typename ProcessorPointPatch::CoupledPointPatch,
        Type
    >(p, iF),
    procPatch_(refCast<const ProcessorPointPatch>(p))
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
ProcessorPointPatchField
(
    const ProcessorPointPatchField
        <PatchField, PointPatch, ProcessorPointPatch, Type>& ptf,
    const PointPatch& p,
    const Field<Type>& iF,
    const PointPatchFieldMapper&
)
:
    CoupledPointPatchField
    <
        PatchField,
        PointPatch,
        typename ProcessorPointPatch::CoupledPointPatch,
        Type
    >(p, iF),
    procPatch_(refCast<const ProcessorPointPatch>(ptf.patch()))
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
ProcessorPointPatchField
(
    const ProcessorPointPatchField
        <PatchField, PointPatch, ProcessorPointPatch, Type>& ptf,
    const Field<Type>& iF
)
:
    CoupledPointPatchField
    <
        PatchField,
        PointPatch,
        typename ProcessorPointPatch::CoupledPointPatch,
        Type
    >(ptf, iF),
    procPatch_(refCast<const ProcessorPointPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
~ProcessorPointPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
initSwapAdd(Field<Type>& pField) const
{
    Field<Type> pf(this->patchInternalField(pField));

    OPstream::write
    (
        procPatch_.neighbProcNo(),
        reinterpret_cast<const char*>(pf.begin()),
        pf.byteSize()
    );
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
swapAdd(Field<Type>& pField) const
{
    Field<Type> pnf(this->size());

    IPstream::read
    (
        procPatch_.neighbProcNo(),
        reinterpret_cast<char*>(pnf.begin()),
        pnf.byteSize()
    );

    if (doTransform())
    {
        const labelList& nonGlobalPatchPoints =
            procPatch_.nonGlobalPatchPoints();

        const processorPolyPatch& ppp = procPatch_.procPolyPatch();
        const labelListList& pointFaces = ppp.pointFaces();
        const tensorField& forwardT = ppp.forwardT();

        if (forwardT.size() == 1)
        {
            transform(pnf, forwardT[0], pnf);
        }
        else
        {
            forAll(nonGlobalPatchPoints, pfi)
            {
                pnf[pfi] = transform
                (
                    forwardT[pointFaces[nonGlobalPatchPoints[pfi]][0]],
                    pnf[pfi]
                );
            }
        }
    }

    addToInternalField(pField, pnf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
