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

Description
    Spherical polar coordinate system.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "sphericalCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sphericalCS, 0);

addToRunTimeSelectionTable(coordinateSystem, sphericalCS, origAxisDir);
addToRunTimeSelectionTable(coordinateSystem, sphericalCS, origRotation);
addToRunTimeSelectionTable(coordinateSystem, sphericalCS, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
sphericalCS::sphericalCS
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
sphericalCS::sphericalCS
(
    const word& name,
    const vector& origin,
    const coordinateRotation& cr
)
:
    cartesianCS(name, origin, cr)
{}


sphericalCS::sphericalCS
(
    const word& name,
    const dictionary& dict
)
:
    cartesianCS(name, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Convert from local coordinate system to the global Cartesian system
vector sphericalCS::toGlobal(const vector& localV) const
{
    scalar r = localV.x();
    scalar theta = localV.y()*mathematicalConstant::pi/180.0;
    scalar phi = localV.z()*mathematicalConstant::pi/180.0;

    return cartesianCS::toGlobal
    (
        vector(r*cos(theta)*sin(phi), r*sin(theta)*sin(phi), r*cos(phi))
    );
}


tmp<vectorField> sphericalCS::toGlobal
(
    const vectorField& localV
) const
{
    const scalarField r = localV.component(vector::X);

    const scalarField theta =
        localV.component(vector::Y)*mathematicalConstant::pi/180.0;

    const scalarField phi =
        localV.component(vector::Z)*mathematicalConstant::pi/180.0;

    vectorField lc(localV.size());
    lc.replace(vector::X, r*cos(theta)*sin(phi));
    lc.replace(vector::Y, r*sin(theta)*sin(phi));
    lc.replace(vector::Z, r*cos(phi));

    return cartesianCS::toGlobal(lc);
}


// Convert from global Cartesian coordinate system to the local system
vector sphericalCS::toLocal(const vector& globalV) const
{
    const vector lc = cartesianCS::toLocal(globalV);
    const scalar r = mag(lc);

    return
        vector
        (
            r,
            atan2(lc.y(), lc.x())*180.0/mathematicalConstant::pi,
            acos(lc.z()/(r + SMALL))*180.0/mathematicalConstant::pi
        );
}


tmp<vectorField> sphericalCS::toLocal
(
    const vectorField& globalV
) const
{
    const vectorField lc = cartesianCS::toLocal(globalV);

    const scalarField r = mag(lc);

    tmp<vectorField> tresult(new vectorField(lc.size()));
    vectorField& result = tresult();

    result.replace
    (
        vector::X, r
        
    );

    result.replace
    (
        vector::Y,
        atan2
        (
            lc.component(vector::Y),
            lc.component(vector::X)
        )*180.0/mathematicalConstant::pi
    );

    result.replace
    (
        vector::Z,
        acos(lc.component(vector::Z)/(r + SMALL))*180.0/mathematicalConstant::pi
    );

    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
