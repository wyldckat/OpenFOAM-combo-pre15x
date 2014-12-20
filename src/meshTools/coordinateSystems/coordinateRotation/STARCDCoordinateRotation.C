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

#include "STARCDCoordinateRotation.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "Switch.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(STARCDCoordinateRotation, 0);

    addToRunTimeSelectionTable
    (
        coordinateRotation, 
        STARCDCoordinateRotation, 
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::STARCDCoordinateRotation::calcTransformations
(
    const scalar rotZ,
    const scalar rotX,
    const scalar rotY,
    const bool inDegrees
)
{
    scalar x = rotX;
    scalar y = rotY;
    scalar z = rotZ;

    if (inDegrees)
    {
        x *= mathematicalConstant::pi/180.0;
        y *= mathematicalConstant::pi/180.0;
        z *= mathematicalConstant::pi/180.0;
    }

    R_ = tensor
    (
        cos(y)*cos(z) - sin(x)*sin(y)*sin(z),
        -cos(x)*sin(z),
        sin(x)*cos(y)*sin(z) + sin(y)*cos(z),

        cos(y)*sin(z) + sin(x)*sin(y)*cos(z),
        cos(x)*cos(z),
        sin(y)*sin(z) - sin(x)*cos(y)*cos(z),

        -cos(x)*sin(y),
        sin(x),
        cos(x)*cos(y)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const vector& rotZrotXrotY,
    const bool inDegrees
)
{
    calcTransformations
    (
        rotZrotXrotY.component(0),
        rotZrotXrotY.component(1),
        rotZrotXrotY.component(2),
        inDegrees
    );
}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const scalar rotZ,
    const scalar rotX,
    const scalar rotY,
    const bool inDegrees
)
{
    calcTransformations( rotZ, rotX, rotY, inDegrees );
}


Foam::STARCDCoordinateRotation::STARCDCoordinateRotation
(
    const dictionary& dict
)
{
    vector rotation(dict.lookup("rotation"));

    bool inDegrees = true;
    if (dict.found("degrees"))
    {
        inDegrees = Switch(dict.lookup("degrees"));
    }

    calcTransformations
    (
        rotation.component(0),
        rotation.component(1),
        rotation.component(2),
        inDegrees
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
