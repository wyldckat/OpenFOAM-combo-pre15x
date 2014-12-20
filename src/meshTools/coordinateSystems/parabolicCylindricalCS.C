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

#include "parabolicCylindricalCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(parabolicCylindricalCS, 0);
    addToRunTimeSelectionTable
    (
        coordinateSystem,
        parabolicCylindricalCS,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parabolicCylindricalCS::parabolicCylindricalCS
(
    const word& name,
    const vector& origin,
    const vector& axis,
    const vector& direction
)
:
    coordinateSystem(name, origin, axis, direction)
{}


Foam::parabolicCylindricalCS::parabolicCylindricalCS
(
    const word& name,
    const vector& origin,
    const coordinateRotation& cr
)
:
    coordinateSystem(name, origin, cr)
{}


Foam::parabolicCylindricalCS::parabolicCylindricalCS
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::parabolicCylindricalCS::toGlobal(const vector& localV) const
{
    // Notation: u = localV.x() v = localV.y() z = localV.z();
    if (localV.y() < 0.0)
    {
        FatalErrorIn("parabolicCylindricalCS::toGlobal(const vector&) const")
            << "parabolic cylindrical coordinates v < 0"
            << abort(FatalError);
    }

    return coordinateSystem::toGlobal
    (
        vector
        (
            0.5*(sqr(localV.x()) - sqr(localV.y())),
            localV.x()*localV.y(),
            localV.z()
        )
    );
}


Foam::tmp<Foam::vectorField> Foam::parabolicCylindricalCS::toGlobal
(
    const vectorField& localV
) const
{
    if (min(localV.component(vector::Y)) < 0.0)
    {
        FatalErrorIn
        (
            "parabolicCylindricalCS::toGlobal(const vectorField&) const"
        )   << "parabolic cylindrical coordinates v < 0"
            << abort(FatalError);
    }

    vectorField lc(localV.size());
    lc.replace
    (
        vector::X,
        0.5*
        (
            sqr(localV.component(vector::X))
          - sqr(localV.component(vector::Y))
        )
    );

    lc.replace
    (
        vector::Y,
        localV.component(vector::X)*localV.component(vector::Y)
    );

    lc.replace
    (
        vector::Z,
        localV.component(vector::Z)
    );

    return coordinateSystem::toGlobal(lc);
}


Foam::vector Foam::parabolicCylindricalCS::toLocal(const vector& globalV) const
{
    notImplemented
    (
        "vector parabolicCylindricalCS::toLocal(const vector& globalV) const"
    );

    return vector::zero;
}


Foam::tmp<Foam::vectorField> Foam::parabolicCylindricalCS::toLocal
(
    const vectorField& globalV
) const
{
    notImplemented
    (
        "tmp<vectorField> parabolicCylindricalCS::toLocal(const "
        "vectorField& globalV) const"
    );

    return tmp<vectorField>(&vectorField::null());
}


// ************************************************************************* //
