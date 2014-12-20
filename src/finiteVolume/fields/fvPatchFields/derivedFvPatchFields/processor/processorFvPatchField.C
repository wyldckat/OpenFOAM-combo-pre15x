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

#include "processorFvPatchField.H"
#include "processorFvPatch.H"
#include "IPstream.H"
#include "OPstream.H"
#include "demandDrivenData.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFvPatch>(p))
{}


template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const Field<Type>& f
)
:
    coupledFvPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFvPatch>(p))
{}


// Construct by mapping given processorFvPatchField<Type>
template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const fvPatch& p,
    const Field<Type>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFvPatch>(p))
{
    if (!isType<processorFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "processorFvPatchField<Type>::processorFvPatchField\n"
            "(\n"
            "    const processorFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorFvPatch>(p))
{
    if (!isType<processorFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "processorFvPatchField<Type>::processorFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not processor type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf
)
:
    lduCoupledInterface(),
    processorLduCoupledInterface(),
    coupledFvPatchField<Type>(ptf),
    procPatch_(refCast<const processorFvPatch>(ptf.patch()))
{}


template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
processorFvPatchField<Type>::~processorFvPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Transfer internal cell data to neighbour processor
// and receive and return neighbour processor internal cell data
template<class Type>
tmp<Field<Type> > processorFvPatchField<Type>::patchNeighbourField() const
{
    return
        (*this - this->patch().weights()*this->patchInternalField())
       /(1.0 - this->patch().weights());
}



// Initialise neighbour processor internal cell data
template<class Type>
void processorFvPatchField<Type>::initEvaluate
(
    const bool bufferdTransfer
)
{
    compressedSend(this->patchInternalField(), bufferdTransfer);
}


// Transfer internal cell data to neighbour processor
// and receive and return neighbour processor internal cell data
template<class Type>
void processorFvPatchField<Type>::evaluate()
{
    tmp<Field<Type> > tpnf(compressedReceive<Type>(this->size()));

    if (doTransform())
    {
        tpnf = transform(procPatch_.forwardT(), tpnf);
    }

    Field<Type>::operator=
    (
        this->patch().weights()*this->patchInternalField()
      + (1.0 - this->patch().weights())*tpnf()
    );
}


// Initialise neighbour colouring transfer
template<class Type>
void processorFvPatchField<Type>::initNbrColour
(
    const labelField& cField,
    const bool bufferdTransfer
) const
{
    send(this->patchInternalField(cField), bufferdTransfer);
}


// Return neighbour colouring
template<class Type>
tmp<labelField> processorFvPatchField<Type>::nbrColour
(
    const labelField&
) const
{
    return receive<label>(this->size());
}


// Initialise neighbour processor internal cell data
template<class Type>
void processorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const bool bufferdTransfer
) const
{
    compressedSend(this->patchInternalField(psiInternal), bufferdTransfer);
}


// Transfer internal cell data to neighbour processor
// and receive and return neighbour processor internal cell data
template<class Type>
void processorFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField&,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt
) const
{
    scalarField pnf(compressedReceive<scalar>(this->size()));

    // Transform according to the transformation tensor
    transformProcCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result

    const labelList::subList FaceCells = this->patch().faceCells();

    forAll(FaceCells, elemI)
    {
        result[FaceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
