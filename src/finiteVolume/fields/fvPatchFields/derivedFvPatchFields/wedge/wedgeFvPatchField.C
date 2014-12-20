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

#include "wedgeFvPatch.H"
#include "wedgeFvPatchField.H"
#include "transformField.H"
#include "symmTransform.H"
#include "diagTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
wedgeFvPatchField<Type>::wedgeFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF
)
:
    transformFvPatchField<Type>(p, iF)
{}


template<class Type>
wedgeFvPatchField<Type>::wedgeFvPatchField
(
    const wedgeFvPatchField<Type>& ptf,
    const fvPatch& p,
    const Field<Type>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<Type>(ptf, p, iF, mapper)
{
    if (!isType<wedgeFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "wedgeFvPatchField<Type>::wedgeFvPatchField\n"
            "(\n"
            "    const wedgeFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
wedgeFvPatchField<Type>::wedgeFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<Type>(p, iF, dict)
{
    if (!isType<wedgeFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "wedgeFvPatchField<Type>::wedgeFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not wedge type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
    if (this->isVolField())
    {
        evaluate();
    }
}


template<class Type>
wedgeFvPatchField<Type>::wedgeFvPatchField
(
    const wedgeFvPatchField<Type>& ptf
)
:
    transformFvPatchField<Type>(ptf)
{}


template<class Type>
wedgeFvPatchField<Type>::wedgeFvPatchField
(
    const wedgeFvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    transformFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return gradient at boundary
template<class Type>
tmp<Field<Type> > wedgeFvPatchField<Type>::snGrad() const
{
    Field<Type> pif = this->patchInternalField();
    return
    (
        transform(refCast<const wedgeFvPatch>(this->patch()).cellT(), pif) - pif
    )*(0.5*this->patch().deltaCoeffs());
}


// Evaluate the patch field
template<class Type>
void wedgeFvPatchField<Type>::evaluate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<Type>::operator==
    (
        transform
        (
            refCast<const wedgeFvPatch>(this->patch()).faceT(),
            this->patchInternalField()
        )
    );
}


// Return defining fields
template<class Type>
tmp<Field<Type> > wedgeFvPatchField<Type>::snGradTransformDiag() const
{
    diagTensor diagT = 
        0.5*diag(I - refCast<const wedgeFvPatch>(this->patch()).cellT());

    return tmp<Field<Type> >
    (
        new Field<Type>
        (
            this->size(),
            transformMask<Type>
            (
                pow
                (
                    reinterpret_cast<const vector&>(diagT),
                    pTraits<typename powProduct<vector, pTraits<Type>::rank>
                    ::type>::zero
                )
            )
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
