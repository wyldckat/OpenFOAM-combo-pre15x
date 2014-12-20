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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "ProcessorPointPatchField.H"
#include "lduMatrix.H"
#include "Map.H"
#include "constraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Raw field sending and receiving
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
template<class Type2>
void ProcessorPointPatchField
<PatchField, PointPatch, ProcessorPointPatch, Type>::sendPointField
(
    const tmp<Field<Type2> >& tf
) const
{
    OPstream toNeighbProc
    (
        procPatch_.neighbProcNo(),
        this->size()*sizeof(Type2)
    );
    toNeighbProc << tf();
    tf.clear();
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
template<class Type2>
void ProcessorPointPatchField
<PatchField, PointPatch, ProcessorPointPatch, Type>::sendEdgeField
(
    const tmp<Field<Type2> >& tf
) const
{
    OPstream toNeighbProc
    (
        procPatch_.neighbProcNo(),
        procPatch_.localEdgeIndices().size()*sizeof(Type2)
    );
    toNeighbProc << tf();
    tf.clear();
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
template<class Type2>
tmp<Field<Type2> >
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
receivePointField() const
{
    IPstream fromNeighbProc
    (
        procPatch_.neighbProcNo(),
        this->size()*sizeof(Type2)
    );
    return tmp<Field<Type2> >(new Field<Type2>(fromNeighbProc));
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
template<class Type2>
tmp<Field<Type2> >
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
receiveEdgeField() const
{
    IPstream fromNeighbProc
    (
        procPatch_.neighbProcNo(),
        procPatch_.localEdgeIndices().size()*sizeof(Type2)
    );
    return tmp<Field<Type2> >(new Field<Type2>(fromNeighbProc));
}


// Initialise diagonal/source update.
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
template<class Type2>
void ProcessorPointPatchField
<PatchField, PointPatch, ProcessorPointPatch, Type>::
initAddFieldTempl
(
    const Field<Type2>& pField
) const
{
    sendPointField(patchInternalField(pField));
}


// Add the diagonal/source to the internal field.
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
template<class Type2>
void ProcessorPointPatchField
<PatchField, PointPatch, ProcessorPointPatch, Type>::
addFieldTempl
(
    Field<Type2>& pField
) const
{
    // Get the neighbour side values
    tmp<Field<Type2> > tpNeighbour = receivePointField<Type2>();
    addToInternalField(pField, tpNeighbour());
}


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
    CoupledPointPatchField<PatchField, PointPatch, Type>(p, iF),
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
    CoupledPointPatchField<PatchField, PointPatch, Type>(p, iF),
    procPatch_(refCast<const ProcessorPointPatch>(p))
{
    if (typeid(p) != typeid(ProcessorPointPatch))
    {
        FatalIOErrorIn
        (
            "ProcessorPointPatchField"
            "<PatchField, PointPatch, ProcessorPointPatch, Type>::"
            "ProcessorPointPatchField\n"
            "(\n"
            "    const PointPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patchMesh().index() << " not processor type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


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
    CoupledPointPatchField<PatchField, PointPatch, Type>(p, iF),
    procPatch_(refCast<const ProcessorPointPatch>(ptf.patchMesh()))
{
    if (typeid(this->patchMesh()) != typeid(ProcessorPointPatch))
    {
        FatalErrorIn
        (
            "ProcessorPointPatchField"
            "<PatchField, PointPatch, ProcessorPointPatch, Type>::"
            "ProcessorPointPatchField\n"
            "(\n"
            "    const ProcessorPointPatchField"
                 "<PatchField, PointPatch, ProcessorPointPatch, Type>& ptf,\n"
            "    const PointPatch& p,\n"
            "    const Field<Type>& iF,\n"
            "    const PointPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patchMesh().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patchMesh().type()
            << exit(FatalError);
    }
}


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
    CoupledPointPatchField<PatchField, PointPatch, Type>(ptf, iF),
    procPatch_(refCast<const ProcessorPointPatch>(ptf.patchMesh()))
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

// Initialise field transfer
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
initAddField() const
{
    initAddFieldTempl(this->internalField());
}


// Add field
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
addField(Field<Type>& f) const
{
    addFieldTempl(f);
}


// Set boundary condition to matrix
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
setBoundaryCondition
(
    Map<constraint<Type> > & fix
) const
{
    // get addressing
    const labelList& meshPoints = procPatch_.meshPoints();

    forAll (meshPoints, pointI)
    {
        const label curPoint = meshPoints[pointI];

        // create a constraint.  None of the components is fixed
        constraint<Type> bc
        (
            curPoint,
            pTraits<Type>::zero,
            pTraits<Type>::zero
        );

        // If pointer is not set, add it, otherwise combine with
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


// Initialise transfer of diagonal coefficients
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
initAddDiag(const scalarField& d) const
{
    initAddFieldTempl(d);
}


// Initialise transfer of source
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
initAddSource(const scalarField& s) const
{
    initAddFieldTempl(s);
}


// Add diagonal coefficients
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
addDiag(scalarField& d) const
{
    addFieldTempl(d);
}


// Add source
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
addSource(scalarField& s) const
{
    addFieldTempl(s);
}


// Initialise neighbour colouring transfer
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
initNbrColour
(
    const labelField& cField,
    const bool
) const
{
    sendPointField(this->patchInternalField(cField));
}


// Return neighbour colouring
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
tmp<labelField>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
nbrColour
(
    const labelField&
) const
{
    return receivePointField<label>();
}


// Initialise transfer of upper/lower coefficients
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
initAddUpperLower
(
    const scalarField& eField
) const
{
    // Gather the data from the given field and send it. Note that the
    // size of the field is equal to the number of edges on the field
    // and NOT the number of points.  

    // Get the addressing
    const labelList& me = procPatch_.localEdgeIndices();

    tmp<scalarField> tresult(new scalarField(me.size()));
    scalarField& result = tresult();

    forAll (me, edgeI)
    {
        result[edgeI] = eField[me[edgeI]];
    }

    // Send the result
    sendEdgeField(tresult);
}


// Add upper/lower coefficients
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
addUpperLower
(
    scalarField& eField
) const
{
    // Add the contribution for the local edge coefficients

    // Get the neighbour side values
    tmp<scalarField> teNeighbour = receiveEdgeField<scalar>();
    scalarField& eNeighbour = teNeighbour();

    // Get the addressing
    const labelList& me = procPatch_.localEdgeIndices();

    forAll (me, edgeI)
    {
        eField[me[edgeI]] += eNeighbour[edgeI];
    }
}


// Get the cut edge coefficients in Amul order
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
tmp<scalarField>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
cutBouCoeffs
(
    const lduMatrix& m
) const
{
    // Go through all the cut edges.  For all owners pick up the upper
    // and for all the neighbours pick up the lower.

    // Get the indices of cut edges
    const labelList& cutOwn = procPatch_.cutEdgeOwnerIndices();
    const labelList& cutNei = procPatch_.cutEdgeNeighbourIndices();
    const labelList& doubleCut = procPatch_.doubleCutEdgeIndices();

    // Get matrix coefficients
    const scalarField& Lower = m.lower();
    const scalarField& Upper = m.upper();

    tmp<scalarField> tcutCoeffs
    (
        new scalarField(cutOwn.size() + cutNei.size() + 2*doubleCut.size(), 0)
    );
    scalarField& cutCoeffs = tcutCoeffs();

    label coeffI = 0;
        
    // Owner side

    forAll (cutOwn, edgeI)
    {
        cutCoeffs[coeffI] = Upper[cutOwn[edgeI]];
        coeffI++;
    }

    // Neighbour side

    forAll (cutNei, edgeI)
    {
        cutCoeffs[coeffI] = Lower[cutNei[edgeI]];
        coeffI++;
    }

    // Doubly cut coeficients

    forAll (doubleCut, edgeI)
    {
        cutCoeffs[coeffI] = Upper[doubleCut[edgeI]];
        coeffI++;

        cutCoeffs[coeffI] = Lower[doubleCut[edgeI]];
        coeffI++;
    }

    return tcutCoeffs;
}


// Get the cut edge coefficients in Tmul order
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
tmp<scalarField>
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
cutIntCoeffs
(
    const lduMatrix& m
) const
{
    // Go through all the cut edges.  For all owners pick up the lower
    // and for all the neighbours pick up the upper.

    // Get the indices of cut edges
    const labelList& cutOwn = procPatch_.cutEdgeOwnerIndices();
    const labelList& cutNei = procPatch_.cutEdgeNeighbourIndices();
    const labelList& doubleCut = procPatch_.doubleCutEdgeIndices();

    // Get matrix coefficients
    const scalarField& Lower = m.lower();
    const scalarField& Upper = m.upper();

    tmp<scalarField> tcutCoeffs
    (
        new scalarField(cutOwn.size() + cutNei.size() + 2*doubleCut.size(), 0)
    );
    scalarField& cutCoeffs = tcutCoeffs();

    label coeffI = 0;

    // Owner side
    // ~~~~~~~~~~
    forAll (cutOwn, edgeI)
    {
        cutCoeffs[coeffI] = Lower[cutOwn[edgeI]];
        coeffI++;
    }

    // Neighbour side
    // ~~~~~~~~~~~~~~
    forAll (cutNei, edgeI)
    {
        cutCoeffs[coeffI] = Upper[cutNei[edgeI]];
        coeffI++;
    }

    // Doubly cut coeficients
    // ~~~~~~~~~~~~~~~~~~~~~~
    forAll (doubleCut, edgeI)
    {
        cutCoeffs[coeffI] = Lower[doubleCut[edgeI]];
        coeffI++;

        cutCoeffs[coeffI] = Upper[doubleCut[edgeI]];
        coeffI++;
    }

    return tcutCoeffs;
}


// Eliminate upper/lower coefficients
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
eliminateUpperLower
(
    scalarField& eField
) const
{
    // Kill the coefficient for cut edges.

    // Get the indices of cut edges
    const labelList& cutOwn = procPatch_.cutEdgeOwnerIndices();
    const labelList& cutNei = procPatch_.cutEdgeNeighbourIndices();
    const labelList& doubleCut = procPatch_.doubleCutEdgeIndices();

    // Owner side
    // ~~~~~~~~~~
    forAll (cutOwn, edgeI)
    {
        eField[cutOwn[edgeI]] = 0;
    }

    // Neighbour side
    // ~~~~~~~~~~~~~~
    forAll (cutNei, edgeI)
    {
        eField[cutNei[edgeI]] = 0;
    }

    // Doubly cut edges
    // ~~~~~~~~~~~~~~~~
    forAll (doubleCut, edgeI)
    {
        eField[doubleCut[edgeI]] = 0;
    }
}


// Initialise matrix update on coupled interfaces
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix& m,
    const scalarField& coeffs,
    const direction,
    const bool
) const
{
    tmp<scalarField> tlocalMult(new scalarField(this->size(), 0));
    scalarField& localMult = tlocalMult();

    const labelList& mp = procPatch_.meshPoints();

    // Get matrix addressing
    const unallocLabelList& L = m.lduAddr().lowerAddr();
    const unallocLabelList& U = m.lduAddr().upperAddr();

    // Note that the addressing is into the local points of the patch.
    // Mesh points is used only for size

    // Get the multiplication mask to exclude all unwanted local multiplies.
    // An example of this is an internal edge between two points which
    // belong to two different processor patches
    const scalarField& cutMask = procPatch_.ownNeiDoubleMask();

    // Coefficients are already ordered in the appropriate way. Just
    // use the counter.  
    label coeffI = 0;

    // Owner side
    // ~~~~~~~~~~
    {
        const labelList& cutOwn = procPatch_.cutEdgeOwnerIndices();
        const labelList& cutOwnStart = procPatch_.cutEdgeOwnerStart();

        forAll (mp, pointI)
        {
            label ownIndex = cutOwnStart[pointI];
            label endOwn = cutOwnStart[pointI + 1];

            for (; ownIndex < endOwn; ownIndex++)
            {
                localMult[pointI] +=
                    coeffs[coeffI]*psiInternal[U[cutOwn[ownIndex]]];

                // Multiply the internal side as well using the cut mask
                result[U[cutOwn[ownIndex]]] +=
                    cutMask[coeffI]*coeffs[coeffI]*psiInternal[mp[pointI]];

                coeffI++;
            }
        }
    }

    // Neighbour side
    // ~~~~~~~~~~~~~~
    {
        const labelList& cutNei = procPatch_.cutEdgeNeighbourIndices();
        const labelList& cutNeiStart = procPatch_.cutEdgeNeighbourStart();

        forAll (mp, pointI)
        {
            label neiIndex = cutNeiStart[pointI];
            label endNei = cutNeiStart[pointI + 1];

            for (; neiIndex < endNei; neiIndex++)
            {
                localMult[pointI] +=
                    coeffs[coeffI]*psiInternal[L[cutNei[neiIndex]]];

                // Multiply the internal side as well using the cut mask
                result[L[cutNei[neiIndex]]] +=
                    cutMask[coeffI]*coeffs[coeffI]*psiInternal[mp[pointI]];

                coeffI++;
            }
        }
    }

    // Doubly cut coefficients
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // There exists a possibility of having an internal edge for a
    // point on the processor patch which is in fact connected to
    // another point of the same patch.  This particular nastiness
    // introduces a deformation in the solution because the edge is
    // either multiplied twice or not at all.  For this purpose, the
    // offending edges need to be separated out and multiplied
    // appropriately.  This will only happen for cell tetrahedral
    // decomposition and is generally nasty.  
    // No need for cut mask here
    {
        const labelList& doubleCut = procPatch_.doubleCutEdgeIndices();

        const labelList& doubleCutOwner = procPatch_.doubleCutOwner();
        const labelList& doubleCutNeighbour = procPatch_.doubleCutNeighbour();

        forAll (doubleCut, edgeI)
        {
            // Owner side
            localMult[doubleCutOwner[edgeI]] +=
                coeffs[coeffI]*psiInternal[U[doubleCut[edgeI]]];
            coeffI++;

            // Neighbour side
            localMult[doubleCutNeighbour[edgeI]] +=
                coeffs[coeffI]*psiInternal[L[doubleCut[edgeI]]];
            coeffI++;
        }
    }

    // Add the local multiplication to this side as well

    forAll (mp, pointI)
    {
        result[mp[pointI]] += localMult[pointI];
    }

    // Send the localMult
    sendPointField(tlocalMult);
}


// Complete matrix update on coupled interfaces
template
<
    template<class> class PatchField,
    class PointPatch,
    class ProcessorPointPatch,
    class Type
>
void
ProcessorPointPatchField<PatchField, PointPatch, ProcessorPointPatch, Type>::
updateInterfaceMatrix
(
    const scalarField&,
    scalarField& result,
    const lduMatrix&,
    const scalarField&,
    const direction
) const
{
    // get the neighbour side multiplication
    tmp<scalarField> tneiMult = receivePointField<scalar>();
    this->addToInternalField(result, tneiMult());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
