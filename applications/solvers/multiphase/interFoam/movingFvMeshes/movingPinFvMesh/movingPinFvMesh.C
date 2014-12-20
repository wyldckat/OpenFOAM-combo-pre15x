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

#include "movingPinFvMesh.H"
#include "addToRunTimeSelectionTable.H"

#include "tetFem.H"
#include "SubField.H"
#include "fvc.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(movingPinFvMesh, 0);

addToRunTimeSelectionTable(movingFvMesh, movingPinFvMesh, IOobject);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from objectRegistry, and read/write options
movingPinFvMesh::movingPinFvMesh(const IOobject& io)
:
    movingFvMesh(io),
    ms_(*this)
{
    // calculate zero-time cell volumes for correct ddt in the first time-sep
    V();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

movingPinFvMesh::~movingPinFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void movingPinFvMesh::move()
{
    vector motionVelocity(0.1*Foam::cos(14.306*time().value()), 0, 0);

    Info << "New boundary velocity: " << motionVelocity.x() << endl;

    ms_.motionU().boundaryField()[0] == motionVelocity;
    ms_.motionU().boundaryField()[1] == motionVelocity;

    movePoints(ms_.newPoints());

    const Time& runTime = time();
    const fvMesh& mesh = *this;

#   include "volContinuity.H"
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
