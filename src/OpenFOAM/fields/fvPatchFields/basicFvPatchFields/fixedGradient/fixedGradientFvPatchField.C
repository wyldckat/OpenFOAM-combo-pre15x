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

Description

\*---------------------------------------------------------------------------*/

#include "fixedGradientFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF
)
:
    fvPatchField<Type>(p, iF),
    gradient_(p.size(), pTraits<Type>::zero)
{
    this->checkVolField();
}


template<class Type>
fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
    const fixedGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const Field<Type>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    gradient_(ptf.gradient_, mapper)
{
    this->checkVolField();
}


template<class Type>
fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF),
    gradient_("gradient", dict, p.size())
{
    this->checkVolField();
    evaluate();
}


template<class Type>
fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
    const fixedGradientFvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    gradient_(ptf.gradient_)
{
    this->checkVolField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
template<class Type>
void fixedGradientFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (m.resizeOnly())
    {
        setSize(m.size());
        gradient_.setSize(m.size());
    }
    else
    {
        Field<Type>::autoMap(m);
        gradient_.autoMap(m);
    }
}


// Reverse-map the given fvPatchField onto this fvPatchField
template<class Type>
void fixedGradientFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const fixedGradientFvPatchField<Type>& fgptf =
        refCast<const fixedGradientFvPatchField<Type> >(ptf);

    gradient_.rmap(fgptf.gradient_, addr);
}


// Evaluate the field on the patch
template<class Type>
void fixedGradientFvPatchField<Type>::evaluate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        this->patchInternalField() + gradient_/this->patch().deltaCoeffs()
    );

    fvPatchField<Type>::evaluate();
}


//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > fixedGradientFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >(new Field<Type>(this->size(), pTraits<Type>::one));
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > fixedGradientFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return gradient()/this->patch().deltaCoeffs();
}

//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > fixedGradientFvPatchField<Type>::
gradientInternalCoeffs() const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > fixedGradientFvPatchField<Type>::
gradientBoundaryCoeffs() const
{
    return gradient();
}


// Write
template<class Type>
void fixedGradientFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    gradient_.writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
