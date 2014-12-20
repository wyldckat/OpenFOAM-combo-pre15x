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

#include "cyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{}


template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const Field<Type>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{
    if (!isType<cyclicFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "cyclicFvPatchField<Type>::cyclicFvPatchField\n"
            "(\n"
            "    const cyclicFvPatchField<Type>& ptf,\n"
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
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{
    if (!isType<cyclicFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "cyclicFvPatchField<Type>::cyclicFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not cyclic type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }

    if (this->isVolField())
    {
        this->evaluate();
    }
}


template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    cyclicPatch_(refCast<const cyclicFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return neighbour coupled given internal cell data
template<class Type>
tmp<Field<Type> > cyclicFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->internalField();
    const labelList::subList FaceCells = cyclicPatch_.faceCells();

    tmp<Field<Type> > tpnf(new Field<Type>(this->size()));
    Field<Type>& pnf = tpnf();

    label sizeby2 = this->size()/2;

    if (doTransform())
    {
        for (label facei=0; facei<sizeby2; facei++)
        {
            pnf[facei] = transform
            (
                cyclicPatch_.forwardT()[0], iField[FaceCells[facei + sizeby2]]
            );

            pnf[facei + sizeby2] = transform
            (
                cyclicPatch_.reverseT()[0], iField[FaceCells[facei]]
            );
        }
    }
    else
    {
        for (label facei=0; facei<sizeby2; facei++)
        {
            pnf[facei] = iField[FaceCells[facei + sizeby2]];
            pnf[facei + sizeby2] = iField[FaceCells[facei]];
        }
    }

    return tpnf;
}


// Return neighbour colouring
template<class Type>
tmp<labelField> cyclicFvPatchField<Type>::nbrColour
(
    const labelField& cField
) const
{
    const labelList::subList FaceCells = this->patch().faceCells();

    tmp<labelField> tpnf(new labelField(this->size()));
    labelField& pnf = tpnf();

    label sizeby2 = this->size()/2;

    for (label facei=0; facei<sizeby2; facei++)
    {
        pnf[facei] = cField[FaceCells[facei + sizeby2]];
        pnf[facei + sizeby2] = cField[FaceCells[facei]];
    }

    return tpnf;
}


// Return matrix product for coupled boundary
template<class Type>
void cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt
) const
{
    scalarField pnf(this->size());

    label sizeby2 = this->size()/2;
    const labelList::subList FaceCells = cyclicPatch_.faceCells();

    for (label facei=0; facei<sizeby2; facei++)
    {
        pnf[facei] = psiInternal[FaceCells[facei + sizeby2]];
        pnf[facei + sizeby2] = psiInternal[FaceCells[facei]];
    }

    // Transform according to the transformation tensors
    transformCyclicCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result
    forAll(FaceCells, elemI)
    {
        result[FaceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
