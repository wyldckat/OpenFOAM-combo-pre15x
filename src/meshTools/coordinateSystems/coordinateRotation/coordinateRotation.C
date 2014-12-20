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

#include "coordinateRotation.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateRotation, 0);
    defineRunTimeSelectionTable(coordinateRotation, dictionary);
}


const Foam::vector Foam::coordinateRotation::z = Foam::vector(0, 0, 1);
const Foam::vector Foam::coordinateRotation::x = Foam::vector(1, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateRotation::coordinateRotation()
:
    R_(tensor::zero)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateRotation> Foam::coordinateRotation::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        Pout<< "coordinateRotation::New(const dictionary&) : "
            << "constructing coordinateRotation"
            << endl;
    }

    word rotType(dict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(rotType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "coordinateRotation::New(const dictionary&)",
            dict
        )   << "Unknown coordinateRotation type " << rotType << nl << nl
            << "Valid coordinateRotation types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<coordinateRotation>(cstrIter()(dict));
}


// ************************************************************************* //
