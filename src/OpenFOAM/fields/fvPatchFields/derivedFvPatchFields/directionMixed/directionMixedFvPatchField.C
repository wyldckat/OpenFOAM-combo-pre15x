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

Class
    directionMixedFvPatchField

Description

\*---------------------------------------------------------------------------*/

#include "directionMixedFvPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
directionMixedFvPatchField<Type>::directionMixedFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF
)
:
    transformFvPatchField<Type>(p, iF),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size())
{
    this->checkVolField();
}


template<class Type>
directionMixedFvPatchField<Type>::directionMixedFvPatchField
(
    const directionMixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const Field<Type>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<Type>(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    refGrad_(ptf.refGrad_, mapper),
    valueFraction_(ptf.valueFraction_, mapper)
{
    this->checkVolField();
}


template<class Type>
directionMixedFvPatchField<Type>::directionMixedFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<Type>(p, iF),
    refValue_("refValue", dict, p.size()),
    refGrad_("refGradient", dict, p.size()),
    valueFraction_("valueFraction", dict, p.size())
{
    this->checkVolField();
    evaluate();
}


template<class Type>
directionMixedFvPatchField<Type>::directionMixedFvPatchField
(
    const directionMixedFvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    transformFvPatchField<Type>(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{
    this->checkVolField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
template<class Type>
void directionMixedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (m.resizeOnly())
    {
        setSize(m.size());
        refValue_.setSize(m.size());
        refGrad_.setSize(m.size());
        valueFraction_.setSize(m.size());
    }
    else
    {
        Field<Type>::autoMap(m);
        refValue_.autoMap(m);
        refGrad_.autoMap(m);
        valueFraction_.autoMap(m);
    }
}


// Reverse-map the given fvPatchField onto this fvPatchField
template<class Type>
void directionMixedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    transformFvPatchField<Type>::rmap(ptf, addr);

    const directionMixedFvPatchField<Type>& dmptf =
        refCast<const directionMixedFvPatchField<Type> >(ptf);

    refValue_.rmap(dmptf.refValue_, addr);
    refGrad_.rmap(dmptf.refGrad_, addr);
    valueFraction_.rmap(dmptf.valueFraction_, addr);
}


// Return gradient at boundary
template<class Type>
tmp<Field<Type> > directionMixedFvPatchField<Type>::snGrad() const
{
    const vectorField& nHat = this->patch().nf();
    Field<Type> pif = this->patchInternalField();

    Field<Type> normalValue = valueFraction_*transform(nHat*nHat, refValue_);

    Field<Type> gradValue = pif + refGrad_/this->patch().deltaCoeffs();

    Field<Type> transformGradValue =
        transform(I - valueFraction_*nHat*nHat, gradValue);

    return
        (normalValue + transformGradValue - pif)*
        this->patch().deltaCoeffs();
}


// Evaluate the field on the patch
template<class Type>
void directionMixedFvPatchField<Type>::evaluate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const vectorField& nHat = this->patch().nf();

    Field<Type> normalValue =
        valueFraction_*transform(nHat*nHat, refValue_);

    Field<Type> gradValue =
        this->patchInternalField() + refGrad_/this->patch().deltaCoeffs();

    Field<Type> transformGradValue =
        transform(I - valueFraction_*nHat*nHat, gradValue);

    Field<Type>::operator=(normalValue + transformGradValue);

    transformFvPatchField<Type>::evaluate();
}


// Return defining fields
template<class Type>
tmp<Field<Type> > directionMixedFvPatchField<Type>::snGradTransformDiag() const
{
    const vectorField& nHat = this->patch().nf();
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return valueFraction_*pow<vector, pTraits<Type>::rank>(diag);
}


// Write
template<class Type>
void directionMixedFvPatchField<Type>::write(Ostream& os) const
{
    transformFvPatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    valueFraction_.writeEntry("valueFraction", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
