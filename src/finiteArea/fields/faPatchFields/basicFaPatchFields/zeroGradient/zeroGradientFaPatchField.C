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

#include "zeroGradientFaPatchField.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
zeroGradientFaPatchField<Type>::zeroGradientFaPatchField
(
    const faPatch& p,
    const Field<Type>& iF
)
:
    faPatchField<Type>(p, iF)
{
    this->checkAreaField();
}


template<class Type>
zeroGradientFaPatchField<Type>::zeroGradientFaPatchField
(
    const zeroGradientFaPatchField<Type>& ptf,
    const faPatch& p,
    const Field<Type>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper)
{
    this->checkAreaField();
}


template<class Type>
zeroGradientFaPatchField<Type>::zeroGradientFaPatchField
(
    const faPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF)
{
    this->checkAreaField();
    faPatchField<Type>::operator=(this->patchInternalField());
}


template<class Type>
zeroGradientFaPatchField<Type>::zeroGradientFaPatchField
(
    const zeroGradientFaPatchField& zgpf,
    const Field<Type>& iF
)
:
    faPatchField<Type>(zgpf, iF)
{
    this->checkAreaField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Evaluate the field on the patch
template<class Type>
void zeroGradientFaPatchField<Type>::evaluate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    operator==(this->patchInternalField());
    faPatchField<Type>::evaluate();
}


//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > zeroGradientFaPatchField<Type>::valueInternalCoeffs
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
tmp<Field<Type> > zeroGradientFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}

//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > zeroGradientFaPatchField<Type>::gradientInternalCoeffs() const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > zeroGradientFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return tmp<Field<Type> >(new Field<Type>(this->size(), pTraits<Type>::zero));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
