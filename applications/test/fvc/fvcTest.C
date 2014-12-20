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
    test

Description
    Finite volume method test code.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    Info<< nl << "Reading field U" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Starting calculations" << endl;

    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ
        ),
        mesh,
        U.dimensions()*dimLength
    );

    volTensorField gradU = fvc::grad(U);

    volScalarField S = dimensionedScalar(0.5)*mag(gradU + gradU.T());

    dimensionedScalar Cs(0.1);

    //volScalarField lSgs = Cs*pow(mesh.V(), 1.0/3.0);

    //nut = sqr(lSgs)*S;

    //nut.write();

    Info<< "end" << endl;
}


// ************************************************************************* //
