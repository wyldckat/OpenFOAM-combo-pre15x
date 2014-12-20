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

#include "uniformFixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_(pTraits<Type>::zero)
{}


template<class Type>
uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
    const uniformFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const Field<Type>& iF,
    const fvPatchFieldMapper&
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_(ptf.uniformValue_)
{
    operator==(uniformValue_);
}


template<class Type>
uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_(pTraits<Type>(dict.lookup("uniformValue")))
{
    operator==(uniformValue_);
}


template<class Type>
uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
    const uniformFixedValueFvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    uniformValue_(ptf.uniformValue_)
{
    operator==(uniformValue_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
template<class Type>
void uniformFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    setSize(m.size());
    operator==(uniformValue_);
}


// Write
template<class Type>
void uniformFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("uniformValue")
        << uniformValue_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
