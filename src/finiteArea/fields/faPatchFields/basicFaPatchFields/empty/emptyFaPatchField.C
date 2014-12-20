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

#include "emptyFaPatchField.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
emptyFaPatchField<Type>::emptyFaPatchField
(
    const faPatch& p,
    const Field<Type>& iF
)
:
    faPatchField<Type>(p, iF, Field<Type>(0))
{}


template<class Type>
emptyFaPatchField<Type>::emptyFaPatchField
(
    const emptyFaPatchField<Type>&,
    const faPatch& p,
    const Field<Type>& iF,
    const faPatchFieldMapper&
)
:
    faPatchField<Type>(p, iF, Field<Type>(0))
{
    if (typeid(this->patchMesh()) != typeid(emptyFaPatch))
    {
        FatalErrorIn
        (
            "emptyFaPatchField<Type>::emptyFaPatchField\n"
            "(\n"
            "    const emptyFaPatchField<Type>&,\n"
            "    const faPatch& p,\n"
            "    const Field<Type>& iF,\n"
            "    const faPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patchMesh().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patchMesh().type()
            << exit(FatalError);
    }
}


template<class Type>
emptyFaPatchField<Type>::emptyFaPatchField
(
    const faPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, Field<Type>(0))
{
    if (typeid(p) != typeid(emptyFaPatch))
    {
        FatalIOErrorIn
        (
            "emptyFaPatchField<Type>::emptyFaPatchField\n"
            "(\n"
            "    const faPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patchMesh().index() << " not empty type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
emptyFaPatchField<Type>::emptyFaPatchField
(
    const emptyFaPatchField<Type>& ptf,
    const Field<Type>& iF
)
:
    faPatchField<Type>(ptf.patchMesh(), iF, Field<Type>(0))
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
