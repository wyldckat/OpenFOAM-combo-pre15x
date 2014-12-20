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
    Cylindrical coordinate system: axis is the local z-axis, direction the
    local x-axis

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "cylindricalCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cylindricalCS, 0);

addToRunTimeSelectionTable(coordinateSystem, cylindricalCS, origAxisDir);
addToRunTimeSelectionTable(coordinateSystem, cylindricalCS, origRotation);
addToRunTimeSelectionTable(coordinateSystem, cylindricalCS, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cylindricalCS::cylindricalCS
(
    const word& name,
    const vector& origin,
    const vector& axis,
    const vector& direction
)
:
    cartesianCS(name, origin, axis, direction)
{}


// Construct from origin and a coordinate rotation
cylindricalCS::cylindricalCS
(
    const word& name,
    const vector& origin,
    const coordinateRotation& cr
)
:
    cartesianCS(name, origin, cr)
{}


cylindricalCS::cylindricalCS
(
    const word& name,
    const dictionary& dict
)
:
    cartesianCS(name, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Convert from local coordinate system to the global Cartesian system
vector cylindricalCS::toGlobal(const vector& localV) const
{
    scalar theta = localV.y()*physicalConstant::pi/180.0;

    return cartesianCS::toGlobal
    (
        vector(localV.x()*cos(theta), localV.x()*sin(theta), localV.z())
    );
}


tmp<vectorField> cylindricalCS::toGlobal
(
    const vectorField& localV
) const
{
    scalarField theta = localV.component(vector::Y)*physicalConstant::pi/180.0;

    vectorField lc(localV.size());
    lc.replace(vector::X, localV.component(vector::X)*cos(theta));
    lc.replace(vector::Y, localV.component(vector::X)*sin(theta));
    lc.replace(vector::Z, localV.component(vector::Z));

    return cartesianCS::toGlobal(lc);
}


// Convert from global Cartesian coordinate system to the local system
vector cylindricalCS::toLocal(const vector& globalV) const
{
    const vector lc = cartesianCS::toLocal(globalV);

    return
        vector
        (
            sqrt(sqr(lc.x()) + sqr(lc.y())),
            atan2(lc.y(),lc.x())*180.0/physicalConstant::pi,
            lc.z()
        );
}


tmp<vectorField> cylindricalCS::toLocal
(
    const vectorField& globalV
) const
{
    const vectorField lc = cartesianCS::toLocal(globalV);

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
        180.0/physicalConstant::pi
    );

    result.replace(vector::Z, lc.component(vector::Z));

    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
