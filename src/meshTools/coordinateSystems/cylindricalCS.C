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

#include "cylindricalCS.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cylindricalCS, 0);

    addToRunTimeSelectionTable(coordinateSystem, cylindricalCS, origAxisDir);
    addToRunTimeSelectionTable(coordinateSystem, cylindricalCS, origRotation);
    addToRunTimeSelectionTable(coordinateSystem, cylindricalCS, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylindricalCS::cylindricalCS
(
    const word& name,
    const vector& origin,
    const vector& axis,
    const vector& direction
)
:
    coordinateSystem(name, origin, axis, direction)
{}


Foam::cylindricalCS::cylindricalCS
(
    const word& name,
    const vector& origin,
    const coordinateRotation& cr
)
:
    coordinateSystem(name, origin, cr)
{}


Foam::cylindricalCS::cylindricalCS
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::cylindricalCS::toGlobal(const vector& localV) const
{
    scalar theta = 
	localV.y()*mathematicalConstant::pi/180.0;

    return coordinateSystem::toGlobal
    (
        vector(localV.x()*cos(theta), localV.x()*sin(theta), localV.z())
    );
}


Foam::tmp<Foam::vectorField> Foam::cylindricalCS::toGlobal
(
    const vectorField& localV
) const
{
    scalarField theta = 
	localV.component(vector::Y)*mathematicalConstant::pi/180.0;

    vectorField lc(localV.size());
    lc.replace(vector::X, localV.component(vector::X)*cos(theta));
    lc.replace(vector::Y, localV.component(vector::X)*sin(theta));
    lc.replace(vector::Z, localV.component(vector::Z));

    return coordinateSystem::toGlobal(lc);
}


Foam::vector Foam::cylindricalCS::toLocal(const vector& globalV) const
{
    const vector lc = coordinateSystem::toLocal(globalV);

    return
        vector
        (
            sqrt(sqr(lc.x()) + sqr(lc.y())),
            atan2(lc.y(),lc.x())*180.0/mathematicalConstant::pi,
            lc.z()
        );
}


Foam::tmp<Foam::vectorField> Foam::cylindricalCS::toLocal
(
    const vectorField& globalV
) const
{
    const vectorField lc = coordinateSystem::toLocal(globalV);

    tmp<vectorField> tresult(new vectorField(lc.size()));
    vectorField& result = tresult();

    result.replace
    (
        vector::X,
        sqrt(sqr(lc.component(vector::X)) + sqr(lc.component(vector::Y)))
    );

    result.replace
    (
        vector::Y,
        atan2(lc.component(vector::Y), lc.component(vector::X))*
        180.0/mathematicalConstant::pi
    );

    result.replace(vector::Z, lc.component(vector::Z));

    return tresult;
}


// ************************************************************************* //
