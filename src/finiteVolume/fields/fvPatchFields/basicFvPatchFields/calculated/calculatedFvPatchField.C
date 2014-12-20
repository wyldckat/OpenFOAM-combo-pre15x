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

#include "calculatedFvPatchField.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
const word& fvPatchField<Type>::calculatedType()
{
    return calculatedFvPatchField<Type>::typeName;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF
)
:
    fvPatchField<Type>(p, iF)
{}


template<class Type>
calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const calculatedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const Field<Type>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, Field<Type>("value", dict, p.size()))
{}


template<class Type>
calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const calculatedFvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    fvPatchField<Type>(ptf, iF)
{}


template<class Type>
template<class Type2>
tmp<fvPatchField<Type> > fvPatchField<Type>::NewCalculatedType
(
    const fvPatchField<Type2>& pf
)
{
    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTablePtr_->find(pf.patch().type());

    if (patchTypeCstrIter != patchConstructorTablePtr_->end())
    {
        return patchTypeCstrIter()
        (
            pf.patch(),
            Field<Type>::null()
        );
    }
    else
    {
        return tmp<fvPatchField<Type> >
        (
            new calculatedFvPatchField<Type>
            (
                pf.patch(),
                Field<Type>::null()
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > calculatedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorIn
    (
        "calculatedFvPatchField<Type>::"
        "valueInternalCoeffs(const tmp<scalarField>&) const"
    )   << "valueInternalCoeffs cannot be called for a "
           "calculatedFvPatchField.\n"
           "You are probably trying to solve for a field with a "
           "calculated boundary conditions."
        << exit(FatalError);

    return *this;
}


template<class Type>
tmp<Field<Type> > calculatedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorIn
    (
        "calculatedFvPatchField<Type>::"
        "valueBoundaryCoeffs(const tmp<scalarField>&) const"
    )   << "valueBoundaryCoeffs cannot be called for a "
           "calculatedFvPatchField.\n"
           "You are probably trying to solve for a field with a "
           "calculated boundary conditions."
        << exit(FatalError);

    return *this;
}

template<class Type>
tmp<Field<Type> > calculatedFvPatchField<Type>::gradientInternalCoeffs() const
{
    FatalErrorIn
    (
        "calculatedFvPatchField<Type>::"
        "gradientInternalCoeffs() const"
    )   << "gradientInternalCoeffs cannot be called for a "
           "calculatedFvPatchField.\n"
           "You are probably trying to solve for a field with a "
           "calculated boundary conditions."
        << exit(FatalError);

    return *this;
}

template<class Type>
tmp<Field<Type> > calculatedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    FatalErrorIn
    (
        "calculatedFvPatchField<Type>::"
        "gradientBoundaryCoeffs() const"
    )   << "gradientBoundaryCoeffs cannot be called for a "
           "calculatedFvPatchField.\n"
           "You are probably trying to solve for a field with a "
           "calculated boundary conditions."
        << exit(FatalError);

    return *this;
}


// Write
template<class Type>
void calculatedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
