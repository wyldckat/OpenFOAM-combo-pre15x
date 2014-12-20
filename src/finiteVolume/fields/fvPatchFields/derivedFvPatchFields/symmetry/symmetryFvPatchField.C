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

#include "symmetryFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
symmetryFvPatchField<Type>::symmetryFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF
)
:
    basicSymmetryFvPatchField<Type>(p, iF)
{}


template<class Type>
symmetryFvPatchField<Type>::symmetryFvPatchField
(
    const symmetryFvPatchField<Type>& ptf,
    const fvPatch& p,
    const Field<Type>& iF,
    const fvPatchFieldMapper& mapper
)
:
    basicSymmetryFvPatchField<Type>(ptf, p, iF, mapper)
{
    if (!isType<symmetryFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "symmetryFvPatchField<Type>::symmetryFvPatchField\n"
            "(\n"
            "    const symmetryFvPatchField<Type>& ptf,\n"
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
symmetryFvPatchField<Type>::symmetryFvPatchField
(
    const fvPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    basicSymmetryFvPatchField<Type>(p, iF, dict)
{
    if (!isType<symmetryFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "symmetryFvPatchField<Type>::symmetryFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not symmetry type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
symmetryFvPatchField<Type>::symmetryFvPatchField
(
    const symmetryFvPatchField<Type>& ptf
)
:
    basicSymmetryFvPatchField<Type>(ptf)
{}


template<class Type>
symmetryFvPatchField<Type>::symmetryFvPatchField
(
    const symmetryFvPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    basicSymmetryFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
