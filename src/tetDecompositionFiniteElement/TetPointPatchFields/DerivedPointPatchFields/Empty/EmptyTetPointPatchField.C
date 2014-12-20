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

\*---------------------------------------------------------------------------*/

#include "EmptyTetPointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch,
    class EmptyPointPatch,
    class Type
>
EmptyTetPointPatchField<PatchField, PointPatch, EmptyPointPatch, Type>::
EmptyTetPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF
)
:
    PatchField<Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class EmptyPointPatch,
    class Type
>
EmptyTetPointPatchField<PatchField, PointPatch, EmptyPointPatch, Type>::
EmptyTetPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const Field<Type>& f
)
:
    PatchField<Type>(p, iF)
{}


template
<
    template<class> class PatchField,
    class PointPatch,
    class EmptyPointPatch,
    class Type
>
EmptyTetPointPatchField<PatchField, PointPatch, EmptyPointPatch, Type>::
EmptyTetPointPatchField
(
    const PointPatch& p,
    const Field<Type>& iF,
    const dictionary& dict
)
:
    PatchField<Type>(p, iF)
{
    if (!isType<EmptyPointPatch>(p))
    {
        FatalIOErrorIn
        (
            "EmptyTetPointPatchField"
            "<PatchField, PointPatch, EmptyPointPatch, Type>::"
            "EmptyTetPointPatchField\n"
            "(\n"
            "    const PointPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not empty type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class EmptyPointPatch,
    class Type
>
EmptyTetPointPatchField<PatchField, PointPatch, EmptyPointPatch, Type>::
EmptyTetPointPatchField
(
    const EmptyTetPointPatchField<PatchField, PointPatch, EmptyPointPatch, Type>&,
    const PointPatch& p,
    const Field<Type>& iF,
    const TetPointPatchFieldMapper&
)
:
    PatchField<Type>(p, iF)
{
    if (!isType<EmptyPointPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "EmptyTetPointPatchField"
            "<PatchField, PointPatch, EmptyPointPatch, Type>::"
            "EmptyTetPointPatchField\n"
            "(\n"
            "    const EmptyTetPointPatchField<PatchField, PointPatch, "
                 "EmptyPointPatch, Type>& ptf,\n"
            "    const PointPatch& p,\n"
            "    const Field<Type>& iF,\n"
            "    const TetPointPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template
<
    template<class> class PatchField,
    class PointPatch,
    class EmptyPointPatch,
    class Type
>
EmptyTetPointPatchField<PatchField, PointPatch, EmptyPointPatch, Type>::
EmptyTetPointPatchField
(
    const EmptyTetPointPatchField<PatchField, PointPatch, EmptyPointPatch, Type>&
        ptf,
    const Field<Type>& iF
)
:
    PatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
