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

#include "processorFvPatchField.H"
#include "processorFvPatch.H"
#include "IPstream.H"
#include "OPstream.H"
#include "demandDrivenData.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Raw field sending and receiving
template<class Type>
template<class Type2>
void processorFvPatchField<Type>::send
(
    const tmp<Field<Type2> >& tf,
    const bool bufferdTransfer
) const
{
    OPstream toNeighbProc
    (
        neighbProcNo(),
        this->size()*sizeof(Type2),
        bufferdTransfer
    );
    toNeighbProc.write((const char*)tf().begin(), this->size()*sizeof(Type2));
}

template<class Type>
template<class Type2>
tmp<Field<Type2> > processorFvPatchField<Type>::receive() const
{
    tmp<Field<Type2> > tf(new Field<Type2>(this->size()));
    Field<Type2>& f = tf();

    IPstream fromNeighbProc(neighbProcNo(), this->size()*sizeof(Type2));
    fromNeighbProc.read((char*)f.begin(), this->size()*sizeof(Type2));

    return tf;
}


// Raw field sending and receiving
template<class Type>
template<class Type2>
void processorFvPatchField<Type>::compressedSend
(
    const tmp<Field<Type2> >& tf,
    const bool bufferdTransfer
) const
{
    if (Pstream::floatTransfer)
    {
        static const label nCmpts = sizeof(Type2)/sizeof(scalar);
        label nm1 = (this->size() - 1)*nCmpts;
        label nlast = sizeof(Type2)/sizeof(float);
        label nFloats = nm1 + nlast;
        label nBytes = nFloats*sizeof(float);

        const Field<Type2>& f = tf();
        const scalar *sArray = (const scalar*)f.begin();
        const scalar *slast = &sArray[nm1];
        float fArray[nFloats];

        for (register label i=0; i<nm1; i++)
        {
            fArray[i] = sArray[i] - slast[i%nCmpts];
        }

        (Type2&)fArray[nm1] = f[this->size() - 1];

        OPstream toNeighbProc(neighbProcNo(), nBytes, bufferdTransfer);
        toNeighbProc.write((const char*)fArray, nBytes);
        tf.clear();
    }
    else
    {
        this->send(tf, bufferdTransfer);
    }
}

template<class Type>
template<class Type2>
tmp<Field<Type2> > processorFvPatchField<Type>::compressedReceive() const
{
    if (Pstream::floatTransfer)
    {
        tmp<Field<Type2> > tf(new Field<Type2>(this->size()));
        Field<Type2>& f = tf();

        static const label nCmpts = sizeof(Type2)/sizeof(scalar);
        label nm1 = (this->size() - 1)*nCmpts;
        label nlast = sizeof(Type2)/sizeof(float);
        label nFloats = nm1 + nlast;
        label nBytes = nFloats*sizeof(float);

        float fArray[nFloats];

        IPstream fromNeighbProc(neighbProcNo(), nBytes);
        fromNeighbProc.read((char*)fArray, nBytes);

        f[this->size() - 1] = (const Type2&)fArray[nm1];
        scalar *sArray = (scalar*)f.begin();
        const scalar *slast = &sArray[nm1];

        for (register label i=0; i<nm1; i++)
        {
            sArray[i] = fArray[i] + slast[i%nCmpts];
        }

        return tf;
    }
    else
    {
        return this->receive<Type2>();
    }
}


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
    if (typeid(this->patchMesh()) != typeid(processorFvPatch))
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
            << this->patchMesh().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patchMesh().type()
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
    if (typeid(p) != typeid(processorFvPatch))
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
        )   << "patch " << this->patchMesh().index() << " not processor type. "
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
    procPatch_(refCast<const processorFvPatch>(ptf.patchMesh()))
{}


template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFvPatch>(ptf.patchMesh()))
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
        (*this - this->patchMesh().weights()*this->patchInternalField())
       /(1.0 - this->patchMesh().weights());
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
    tmp<Field<Type> > tpnf(compressedReceive<Type>());

    if (doTransform())
    {
        tpnf = transform(procPatch_.forwardT(), tpnf);
    }

    Field<Type>::operator=
    (
        this->patchMesh().weights()*this->patchInternalField()
      + (1.0 - this->patchMesh().weights())*tpnf()
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
    return receive<label>();
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
    scalarField pnf(compressedReceive<scalar>());

    // Transform according to the transformation tensor
    transformProcCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result

    const labelList::subList FaceCells = this->patchMesh().faceCells();

    forAll(FaceCells, elemI)
    {
        result[FaceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
