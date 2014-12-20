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

#include "slicedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
slicedFvPatchField<Type>::slicedFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const Field<Type>& completeField
)
:
    fvPatchField<Type>(p, iF, Field<Type>())
{
    // Set the fvPatchField to a slice of the given complete field
    UList<Type>::operator=(p.patchSlice(completeField));
}


template<class Type>
slicedFvPatchField<Type>::slicedFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF
)
:
    fvPatchField<Type>(p, iF)
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "slicedFvPatchField(const fvPatch&, const Field<Type>&)"
    );
}


template<class Type>
slicedFvPatchField<Type>::slicedFvPatchField
(
    const slicedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const Field<Type>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper)
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "slicedFvPatchField(const slicedFvPatchField<Type>&, "
        "const fvPatch&, const Field<Type>&, const fvPatchFieldMapper&)"
    );
}


template<class Type>
slicedFvPatchField<Type>::slicedFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, Field<Type>("value", dict, p.size()))
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "slicedFvPatchField(const Field<Type>&, const dictionary&)"
    );
}


template<class Type>
slicedFvPatchField<Type>::slicedFvPatchField
(
    const slicedFvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    fvPatchField<Type>(ptf.patch(), iF, Field<Type>())
{
    // Transfer the slice from the argument
    UList<Type>::operator=(ptf);
}

template<class Type>
tmp<fvPatchField<Type> > slicedFvPatchField<Type>::clone() const
{
    return tmp<fvPatchField<Type> >
    (
        new slicedFvPatchField<Type>(*this)
    );
}


template<class Type>
slicedFvPatchField<Type>::slicedFvPatchField
(
    const slicedFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf.patch(), ptf.internalField(), Field<Type>())
{
    // Transfer the slice from the argument
    UList<Type>::operator=(ptf);
}


template<class Type>
tmp<fvPatchField<Type> > slicedFvPatchField<Type>::clone
(
    const Field<Type>& iF
) const
{
    return tmp<fvPatchField<Type> >
    (
        new slicedFvPatchField<Type>(*this, iF)
    );
}


template<class Type>
slicedFvPatchField<Type>::~slicedFvPatchField<Type>()
{
    // Set the fvPatchField storage pointer to NULL before it's destruction
    // to protect the field it a slice of.
    UList<Type>::operator=(UList<Type>(NULL, 0));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return patch-normal gradient
template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::snGrad() const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "snGrad()"
    );

    return Field<Type>::null();
}

//- Update the coefficients associated with the patch field
//  Sets Updated to true
template<class Type>
void slicedFvPatchField<Type>::updateCoeffs()
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "updateCoeffs()"
    );
}

//- Return internal field next to patch as patch field
template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::patchInternalField() const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "patchInternalField()"
    );

    return Field<Type>::null();
}

//- Return neighbour coupled given internal cell data
template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::patchNeighbourField
(
    const Field<Type>& iField
) const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "patchNeighbourField(const Field<Type>& iField)"
    );

    return Field<Type>::null();
}


//- Return patchField of the values on the patch or on the
//  opposite patch
template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::patchNeighbourField() const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "patchNeighbourField()"
    );

    return Field<Type>::null();
}

//- Initialise the evaluation of the patch field
template<class Type>
void slicedFvPatchField<Type>::initEvaluate(const bool)
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "initEvaluate(const bool)"
    );
}

//- Evaluate the patch field, sets Updated to false
template<class Type>
void slicedFvPatchField<Type>::evaluate()
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "evaluate()"
    );
}


//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "valueInternalCoeffs(const tmp<scalarField>&)"
    );

    return Field<Type>::null();
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "valueBoundaryCoeffs(const tmp<scalarField>&)"
    );

    return Field<Type>::null();
}

//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::gradientInternalCoeffs() const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "gradientInternalCoeffs()"
    );

    return Field<Type>::null();
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "gradientBoundaryCoeffs()"
    );

    return Field<Type>::null();
}


// Write
template<class Type>
void slicedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
