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

Application
    engineSwirl

Description
    Generates a swirling flow for engine calulations

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "engineTime.H"
#include "physicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createEngineTime.H"
#   include "createMesh.H"
#   include "createEngineMovingMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Reading combustion properties\n" << endl;

    IOdictionary engineGeometry
    (
        IOobject
        (
            "engineGeometry",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedVector swirlAxis
    (
        engineGeometry.lookup("swirlAxis")
    );
    swirlAxis = swirlAxis/mag(swirlAxis);

    dimensionedVector swirlCenter
    (
        engineGeometry.lookup("swirlCenter")
    );

    dimensionedScalar swirlRPMRatio
    (
        engineGeometry.lookup("swirlRPMRatio")
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    dimensionedVector Omega =
        (physicalConstant::pi*swirlRPMRatio*runTime.rpm()/30)*swirlAxis;

    volVectorField r = 
        (mesh.C() - swirlCenter)
      - swirlAxis*(swirlAxis & (mesh.C() - swirlCenter));

    volScalarField Uz = runTime.pistonSpeed()
        *min((deckHeight - mesh.C().component(vector::Z))
        /(deckHeight - pistonPosition), 1.0);
    
    U = (Omega ^ r) + vector(0, 0, 1)*Uz;

    U.write();

    Info << "\n end\n";

    return(0);
}


// ************************************************************************* //
