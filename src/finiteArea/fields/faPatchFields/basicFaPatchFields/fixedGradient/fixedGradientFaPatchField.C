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

#include "fixedGradientFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
fixedGradientFaPatchField<Type>::fixedGradientFaPatchField
(
    const faPatch& p,
    const Field<Type>& iF
)
:
    faPatchField<Type>(p, iF),
    gradient_(p.size())
{}


template<class Type>
fixedGradientFaPatchField<Type>::fixedGradientFaPatchField
(
    const fixedGradientFaPatchField<Type>& ptf,
    const faPatch& p,
    const Field<Type>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper),
    gradient_(ptf.gradient_, mapper)
{}


template<class Type>
fixedGradientFaPatchField<Type>::fixedGradientFaPatchField
(
    const faPatch& p,
    const Field<Type>& iF,
    Istream& is
)
:
    faPatchField<Type>(p, iF, is),
    gradient_(is)
{
    evaluate();
}


template<class Type>
fixedGradientFaPatchField<Type>::fixedGradientFaPatchField
(
    const faPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF),
    gradient_("gradient", dict, p.size())
{
    evaluate();
}


template<class Type>
fixedGradientFaPatchField<Type>::fixedGradientFaPatchField
(
    const fixedGradientFaPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    faPatchField<Type>(ptf, iF),
    gradient_(ptf.gradient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Evaluate the field on the patch
template<class Type>
void fixedGradientFaPatchField<Type>::evaluate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        this->patchInternalField() + gradient_/this->patch().deltaCoeffs()
    );

    faPatchField<Type>::evaluate();
}


//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > fixedGradientFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >
        (
            new Field<Type>(this->size(), pTraits<Type>::one)
        );
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > fixedGradientFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return gradient()/this->patch().deltaCoeffs();
}

//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > fixedGradientFaPatchField<Type>::gradientInternalCoeffs() const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > fixedGradientFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return gradient();
}


// Write
template<class Type>
void fixedGradientFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    gradient_.writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
