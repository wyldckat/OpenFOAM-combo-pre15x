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

Class
    select

Description

\*---------------------------------------------------------------------------*/

#include "interpolation.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
autoPtr<interpolation<Type> > interpolation<Type>::New
(
    const dictionary& interpolationSchemes,
    const volPointInterpolation& pInterp,
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
    word interpolationType
    (
        interpolationSchemes.lookup(psi.name())
    );

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_
            ->find(interpolationType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "interpolation::New(const dictionary&, "
            "const GeometricField<Type, fvPatchField, volMesh>&)",
            interpolationSchemes
        )   << "Unknown interpolation type " << interpolationType
            << endl << endl
            << "Valid interpolation types : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<interpolation<Type> >(cstrIter()(pInterp, psi));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
