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

#include "fixedValueFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
fixedValueFaPatchField<Type>::fixedValueFaPatchField
(
    const faPatch& p,
    const Field<Type>& iF
)
:
    faPatchField<Type>(p, iF)
{}


template<class Type>
fixedValueFaPatchField<Type>::fixedValueFaPatchField
(
    const fixedValueFaPatchField<Type>& ptf,
    const faPatch& p,
    const Field<Type>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
fixedValueFaPatchField<Type>::fixedValueFaPatchField
(
    const faPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, Field<Type>("value", dict, p.size()))
{}


template<class Type>
fixedValueFaPatchField<Type>::fixedValueFaPatchField
(
    const fixedValueFaPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    faPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > fixedValueFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > fixedValueFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return *this;
}

//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > fixedValueFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > fixedValueFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return this->patch().deltaCoeffs()*(*this);
}


// Write
template<class Type>
void fixedValueFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
