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

Global
    createPhi

Description
    Creates and initialises the face-flux field phi.

\*---------------------------------------------------------------------------*/

#ifndef createPhi_H
#define createPhi_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    IOobject phiHeader
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );


    surfaceScalarField* phiPtr;

    if (phiHeader.headerOk())
    {
        Info<< "Reading face flux field phi\n" << endl;

        phiPtr = new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
    }
    else
    {
        Info<< "Calculating face flux field phi\n" << endl;

        wordList phiTypes
        (
            U.boundaryField().size(),
            calculatedPatchScalarField::typeName
        );

        for (label i=0; i<U.boundaryField().size(); i++)
        {
            if (!p.boundaryField()[i].fixedValue())
            {
                phiTypes[i] = fixedValuePatchScalarField::typeName;
            }
        }

        phiPtr = new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(rho*U) & mesh.Sf(),
            phiTypes
        );
    }

    surfaceScalarField& phi = *phiPtr;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
