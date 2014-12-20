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
    Source file implementing general coordinate transformations
    from arbitary coordinate systems to the global Cartesian system.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "cartesianCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cartesianCS, 0);

addToRunTimeSelectionTable(coordinateSystem, cartesianCS, origAxisDir);
addToRunTimeSelectionTable(coordinateSystem, cartesianCS, origRotation);
addToRunTimeSelectionTable(coordinateSystem, cartesianCS, dictionary);

const scalar cartesianCS::nonOrthogonalError = 1.0e-8;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void cartesianCS::calcTransformations()
{
    vector ip = direction()/mag(direction());
    vector kp = axis()/mag(axis());
    vector jp = kp ^ ip;

    if (mag(kp&ip)/(mag(kp)*mag(ip)) >= nonOrthogonalError)
    {
        FatalErrorIn("void cartesianCS::calcTransformations()")
            << "Coordinate system not orthogonal" << endl
            << "mag(kp & ip) = " << mag(kp & ip)
            << abort(FatalError);
    }

    R_ = tensor
    (
        ip.x(), ip.y(), ip.z(),
        jp.x(), jp.y(), jp.z(),
        kp.x(), kp.y(), kp.z()
    ).T();

    Rtrans_ = R_.T();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from origin and two axes
cartesianCS::cartesianCS
(
    const word& name,
    const vector& origin,
    const vector& axis,
    const vector& dir
)
:
    coordinateSystem(name, origin, axis, dir),
    R_(tensor::zero),
    Rtrans_(tensor::zero)
{
    calcTransformations();
}


// Construct from origin and a coordinate rotation
cartesianCS::cartesianCS
(
    const word& name,
    const vector& origin,
    const coordinateRotation& cr
)
:
    coordinateSystem(name, origin, cr),
    R_(cr.R()),
    Rtrans_(R_.T())
{}


cartesianCS::cartesianCS
(
    const word& name,
    const dictionary& dict
)
:
    coordinateSystem(name, dict),
    R_(tensor::zero),
    Rtrans_(tensor::zero)
{
    calcTransformations();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector cartesianCS::toGlobal(const vector& localV) const
{
    return (R_ & localV) + origin();
}


tmp<vectorField> cartesianCS::toGlobal(const vectorField& localV) const
{
    return (R_ & localV) + origin();
}


vector cartesianCS::toLocal(const vector& globalV) const
{
    return (Rtrans_ & (globalV - origin()));
}


tmp<vectorField> cartesianCS::toLocal(const vectorField& globalV) const
{
    return (Rtrans_ & (globalV - origin()));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
