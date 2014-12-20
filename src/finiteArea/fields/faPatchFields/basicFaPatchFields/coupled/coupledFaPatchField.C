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

#include "coupledFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const Field<Type>& iF
)
:
    faPatchField<Type>(p, iF)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const Field<Type>& iF,
    const Field<Type>& f
)
:
    faPatchField<Type>(p, iF, f)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf,
    const faPatch& p,
    const Field<Type>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, dict)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    faPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// return gradient at boundary
template<class Type>
tmp<Field<Type> > coupledFaPatchField<Type>::snGrad() const
{
    return
        (patchNeighbourField() - this->patchInternalField())
       *this->patchMesh().deltaCoeffs();
}


// Return neighbour coupled internal cell data
template<class Type>
void coupledFaPatchField<Type>::initPatchNeighbourField() const
{
    this->initPatchNeighbourField(this->internalField());
}


// Return neighbour coupled internal cell data
template<class Type>
tmp<Field<Type> > coupledFaPatchField<Type>::patchNeighbourField() const
{
    return this->patchNeighbourField(this->internalField());
}


//- Initialise the evaluation of the patch field
template<class Type>
void coupledFaPatchField<Type>::initEvaluate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    initPatchNeighbourField();
}


// Evaluate the patch field
template<class Type>
void coupledFaPatchField<Type>::evaluate()
{
    Field<Type>::operator=
    (
        this->patchMesh().weights()*this->patchInternalField()
      + (1.0 - this->patchMesh().weights())*patchNeighbourField()
    );
}


//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > coupledFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*w;
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > coupledFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*(1.0 - w);
}

//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > coupledFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*this->patchMesh().deltaCoeffs();
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > coupledFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return -gradientInternalCoeffs();
}


// Write
template<class Type>
void coupledFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
