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

#include "fixedNormalSlipFvPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF
)
:
    transformFvPatchField<Type>(p, iF),
    fixedValue_(p.size(), pTraits<Type>::zero)
{
    this->checkVolField();
}


template<class Type>
fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fixedNormalSlipFvPatchField<Type>& ptf,
    const fvPatch& p,
    const Field<Type>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<Type>(ptf, p, iF, mapper),
    fixedValue_(ptf.fixedValue_, mapper)
{
    this->checkVolField();
}


template<class Type>
fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<Type>(p, iF),
    fixedValue_("fixedValue", dict, p.size())
{
    this->checkVolField();
    evaluate();
}


template<class Type>
fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fixedNormalSlipFvPatchField<Type>& ptf
)
:
    transformFvPatchField<Type>(ptf),
    fixedValue_(ptf.fixedValue_)
{
    this->checkVolField();
}


template<class Type>
fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fixedNormalSlipFvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    transformFvPatchField<Type>(ptf, iF),
    fixedValue_(ptf.fixedValue_)
{
    this->checkVolField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
template<class Type>
void fixedNormalSlipFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
    fixedValue_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
template<class Type>
void fixedNormalSlipFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    transformFvPatchField<Type>::rmap(ptf, addr);

    const fixedNormalSlipFvPatchField<Type>& dmptf =
        refCast<const fixedNormalSlipFvPatchField<Type> >(ptf);

    fixedValue_.rmap(dmptf.fixedValue_, addr);
}


// Return gradient at boundary
template<class Type>
tmp<Field<Type> > fixedNormalSlipFvPatchField<Type>::snGrad() const
{
    vectorField nHat = this->patch().nf();
    Field<Type> pif = this->patchInternalField();

    return
    (
        (nHat*(nHat & fixedValue_) + transform(I - sqr(nHat), pif)) - pif
    )*this->patch().deltaCoeffs();
}


// Evaluate the field on the patch
template<class Type>
void fixedNormalSlipFvPatchField<Type>::evaluate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField nHat = this->patch().nf();

    Field<Type>::operator=
    (
        nHat*(nHat & fixedValue_)
      + transform(I - sqr(nHat), this->patchInternalField())
    );

    transformFvPatchField<Type>::evaluate();
}


// Return defining fields
template<class Type>
tmp<Field<Type> > fixedNormalSlipFvPatchField<Type>::snGradTransformDiag() const
{
    vectorField nHat = this->patch().nf();
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


// Write
template<class Type>
void fixedNormalSlipFvPatchField<Type>::write(Ostream& os) const
{
    transformFvPatchField<Type>::write(os);
    fixedValue_.writeEntry("fixedValue", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
