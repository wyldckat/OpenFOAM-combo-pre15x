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

#include "liftDrag.H"

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::dragCoefficient
(
    const volVectorField& U,
    const volScalarField& p,
    const dimensionedScalar& mu,
    const word& patchName,
    const vector& Uinf,
    const scalar& Aref
)
{
    if (mag(Uinf) < VSMALL)
    {
        FatalErrorIn
        (
            "scalar dragCoefficient\n"
            "(\n"
            "    const volVectorField& U,\n"
            "    const volScalarField& p,\n"
            "    const dimensionedScalar& mu,\n"
            "    const word& patchName,\n"
            "    const vector& Uinf,\n"
            "    const scalar& Aref\n"
            ")\n"
        )   << "Uinf is zero."
            << exit(FatalError);
    }

    vector flowDirection = Uinf/mag(Uinf);

    const fvMesh& mesh = p.mesh();

    label patchLabel = -1;

    forAll (mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].name() == patchName)
        {
            patchLabel = patchi;
            break;
        }
    }

    vector pressureForce;
    vector viscousForce;

    if (patchLabel != -1)
    {
        pressureForce = sum
        (
            p.boundaryField()[patchLabel]
            *mesh.Sf().boundaryField()[patchLabel]
        );

        viscousForce = sum
        (
          - mu.value()*U.boundaryField()[patchLabel].snGrad()*
            mesh.magSf().boundaryField()[patchLabel]
        );
    }
    else
    {
        pressureForce = sum(vectorField(0));
        viscousForce = sum(vectorField(0));
    }

    scalar qRef = 0.5*magSqr(Uinf);
    scalar Fref = qRef*Aref;

    vector pressureCoeff = pressureForce/Fref;
    vector viscousCoeff = viscousForce/Fref;


    return ((pressureCoeff + viscousCoeff) & flowDirection);
}


Foam::scalar Foam::turbDragCoefficient
(
    const autoPtr<turbulenceModel>& turbulence,
    const volVectorField& U,
    const volScalarField& p,
    const dimensionedScalar& mu,
    const word& patchName,
    const vector& Uinf,
    const scalar& Aref
)
{
    if (mag(Uinf) < VSMALL)
    {
        FatalErrorIn
        (
            "scalar turbDragCoefficient\n"
            "(\n"
            "    const autoPtr<turbulenceModel>& turbulence,\n"
            "    const volVectorField& U,\n"
            "    const volScalarField& p,\n"
            "    const dimensionedScalar& mu,\n"
            "    const word& patchName,\n"
            "    const vector& Uinf,\n"
            "    const scalar& Aref\n"
            ")\n"
        )   << "Uinf is zero."
            << exit(FatalError);
    }

    vector flowDirection = Uinf/mag(Uinf);

    const fvMesh& mesh = p.mesh();

    label patchLabel = -1;

    forAll (mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].name() == patchName)
        {
            patchLabel = patchi;
            break;
        }
    }

    vector pressureForce;
    vector viscousForce;
    vector turbForce;

    if (patchLabel != -1)
    {
        pressureForce = sum
        (
            p.boundaryField()[patchLabel]
            *mesh.Sf().boundaryField()[patchLabel]
        );

        viscousForce = sum
        (
          - mu.value()*U.boundaryField()[patchLabel].snGrad()*
            mesh.magSf().boundaryField()[patchLabel]
        );

        turbForce = sum
        (
          - mesh.Sf().boundaryField()[patchLabel]
          & turbulence->R()().boundaryField()[patchLabel]
        );
    }
    else
    {
        pressureForce = sum(vectorField(0));
        viscousForce = sum(vectorField(0));
        turbForce = sum(vectorField(0));
    }

    scalar qRef = 0.5*magSqr(Uinf);
    scalar Fref = qRef*Aref;

    vector pressureCoeff = pressureForce/Fref;
    vector viscousCoeff = viscousForce/Fref;
    vector turbCoeff = turbForce/Fref;


    return ((pressureCoeff + viscousCoeff + turbCoeff) & flowDirection);
}


Foam::scalar Foam::lesDragCoefficient
(
    const autoPtr<LESmodel>& sgsModel,
    const volVectorField& U,
    const volScalarField& p,
    const dimensionedScalar& mu,
    const word& patchName,
    const vector& Uinf,
    const scalar& Aref
)
{
    if (mag(Uinf) < VSMALL)
    {
        FatalErrorIn
        (
            "scalar lesDragCoefficient\n"
            "(\n"
            "    const autoPtr<LESmodel>& sgsModel,\n"
            "    const volVectorField& U,\n"
            "    const volScalarField& p,\n"
            "    const dimensionedScalar& mu,\n"
            "    const word& patchName,\n"
            "    const vector& Uinf,\n"
            "    const scalar& Aref\n"
            ")\n"
        )   << "Uinf is zero."
            << exit(FatalError);
    }

    vector flowDirection = Uinf/mag(Uinf);

    const fvMesh& mesh = p.mesh();

    label patchLabel = -1;

    forAll (mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].name() == patchName)
        {
            patchLabel = patchi;
            break;
        }
    }

    vector pressureForce;
    vector viscousForce;
    vector turbForce;

    if (patchLabel != -1)
    {
        pressureForce = sum
        (
            p.boundaryField()[patchLabel]
            *mesh.Sf().boundaryField()[patchLabel]
        );

        viscousForce = sum
        (
          - mu.value()*U.boundaryField()[patchLabel].snGrad()*
            mesh.magSf().boundaryField()[patchLabel]
        );

        turbForce = sum
        (
          - mesh.Sf().boundaryField()[patchLabel]
          & sgsModel->B()().boundaryField()[patchLabel]
        );
    }
    else
    {
        pressureForce = sum(vectorField(0));
        viscousForce = sum(vectorField(0));
        turbForce = sum(vectorField(0));
    }

    scalar qRef = 0.5*magSqr(Uinf);
    scalar Fref = qRef*Aref;

    vector pressureCoeff = pressureForce/Fref;
    vector viscousCoeff = viscousForce/Fref;
    vector turbCoeff = turbForce/Fref;


    return ((pressureCoeff + viscousCoeff + turbCoeff) & flowDirection);
}


Foam::vector Foam::liftCoefficient
(
    const volVectorField& U,
    const volScalarField& p,
    const dimensionedScalar& mu,
    const word& patchName,
    const vector& Uinf,
    const scalar& Aref
)
{
    if (mag(Uinf) < VSMALL)
    {
        FatalErrorIn
        (
            "vector liftCoefficient\n"
            "(\n"
            "    const volVectorField& U,\n"
            "    const volScalarField& p,\n"
            "    const dimensionedScalar& mu,\n"
            "    const word& patchName,\n"
            "    const vector& Uinf,\n"
            "    const scalar& Aref\n"
            ")\n"
        )   << "Uinf is zero."
            << exit(FatalError);
    }

    vector flowDirection = Uinf/mag(Uinf);

    const fvMesh& mesh = p.mesh();

    label patchLabel = -1;

    forAll (mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].name() == patchName)
        {
            patchLabel = patchi;
            break;
        }
    }

    vector pressureForce;
    vector viscousForce;

    if (patchLabel != -1)
    {
        pressureForce = sum
        (
            p.boundaryField()[patchLabel]
            *mesh.Sf().boundaryField()[patchLabel]
        );

        viscousForce = sum
        (
            -mu.value()*U.boundaryField()[patchLabel].snGrad()
            *mesh.magSf().boundaryField()[patchLabel]
        );
    }
    else
    {
        pressureForce = sum(vectorField(0));
        viscousForce = sum(vectorField(0));
    }

    scalar qRef = 0.5*magSqr(Uinf);
    scalar Fref = qRef*Aref;

    vector pressureCoeff = pressureForce/Fref;
    vector viscousCoeff = viscousForce/Fref;

    return
        (pressureCoeff + viscousCoeff)
      - flowDirection*(flowDirection & (pressureCoeff + viscousCoeff));
}


Foam::vector Foam::turbLiftCoefficient
(
    const autoPtr<Foam::turbulenceModel>& turbulence,
    const volVectorField& U,
    const volScalarField& p,
    const dimensionedScalar& mu,
    const word& patchName,
    const vector& Uinf,
    const scalar& Aref
)
{
    if (mag(Uinf) < VSMALL)
    {
        FatalErrorIn
        (
            "vector turbLiftCoefficient\n"
            "(\n"
            "    const autoPtr<Foam::turbulenceModel>& turbulence,\n"
            "    const volVectorField& U,\n"
            "    const volScalarField& p,\n"
            "    const dimensionedScalar& mu,\n"
            "    const word& patchName,\n"
            "    const vector& Uinf,\n"
            "    const scalar& Aref\n"
            ")\n"
        )   << "Uinf is zero."
            << exit(FatalError);
    }

    vector flowDirection = Uinf/mag(Uinf);

    const fvMesh& mesh = p.mesh();

    label patchLabel = -1;

    forAll (mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].name() == patchName)
        {
            patchLabel = patchi;
            break;
        }
    }

    vector pressureForce;
    vector viscousForce;
    vector turbForce;

    if (patchLabel != -1)
    {
        pressureForce = sum
        (
            p.boundaryField()[patchLabel]
            *mesh.Sf().boundaryField()[patchLabel]
        );

        viscousForce = sum
        (
          - mu.value()*U.boundaryField()[patchLabel].snGrad()*
            mesh.magSf().boundaryField()[patchLabel]
        );

        turbForce = sum
        (
          - mesh.Sf().boundaryField()[patchLabel]
          & turbulence->R()().boundaryField()[patchLabel]
        );
    }
    else
    {
        pressureForce = sum(vectorField(0));
        viscousForce = sum(vectorField(0));
        turbForce = sum(vectorField(0));
    }

    scalar qRef = 0.5*magSqr(Uinf);
    scalar Fref = qRef*Aref;

    vector pressureCoeff = pressureForce/Fref;
    vector viscousCoeff = viscousForce/Fref;
    vector turbCoeff = turbForce/Fref;

    return
        (pressureCoeff + viscousCoeff + turbCoeff)
      - flowDirection
      * (flowDirection & (pressureCoeff + viscousCoeff + turbForce));
}


Foam::vector Foam::lesLiftCoefficient
(
    const autoPtr<LESmodel>& sgsModel,
    const volVectorField& U,
    const volScalarField& p,
    const dimensionedScalar& mu,
    const word& patchName,
    const vector& Uinf,
    const scalar& Aref
)
{
    if (mag(Uinf) < VSMALL)
    {
        FatalErrorIn
        (
            "vector lesLiftCoefficient\n"
            "(\n"
            "    const autoPtr<LESmodel>& sgsModel,\n"
            "    const volVectorField& U,\n"
            "    const volScalarField& p,\n"
            "    const dimensionedScalar& mu,\n"
            "    const word& patchName,\n"
            "    const vector& Uinf,\n"
            "    const scalar& Aref\n"
            ")\n"
        )   << "Uinf is zero."
            << exit(FatalError);
    }

    vector flowDirection = Uinf/mag(Uinf);

    const fvMesh& mesh = p.mesh();

    label patchLabel = -1;

    forAll (mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].name() == patchName)
        {
            patchLabel = patchi;
            break;
        }
    }

    vector pressureForce;
    vector viscousForce;
    vector turbForce;

    if (patchLabel != -1)
    {
        pressureForce = sum
        (
            p.boundaryField()[patchLabel]
            *mesh.Sf().boundaryField()[patchLabel]
        );

        viscousForce = sum
        (
          - mu.value()*U.boundaryField()[patchLabel].snGrad()*
            mesh.magSf().boundaryField()[patchLabel]
        );

        turbForce = sum
        (
          - mesh.Sf().boundaryField()[patchLabel]
          & sgsModel->B()().boundaryField()[patchLabel]
        );
    }
    else
    {
        pressureForce = sum(vectorField(0));
        viscousForce = sum(vectorField(0));
        turbForce = sum(vectorField(0));
    }

    scalar qRef = 0.5*magSqr(Uinf);
    scalar Fref = qRef*Aref;

    vector pressureCoeff = pressureForce/Fref;
    vector viscousCoeff = viscousForce/Fref;
    vector turbCoeff = turbForce/Fref;

    return
        (pressureCoeff + viscousCoeff + turbCoeff)
      - flowDirection
      * (flowDirection & (pressureCoeff + viscousCoeff + turbForce));
}


Foam::vector Foam::momentCoefficient
(
    const volVectorField& U,
    const volScalarField& p,
    const dimensionedScalar& mu,
    const word& patchName,
    const vector& Uinf,
    const scalar& Aref,
    const scalar& Lref
)
{
    if (mag(Uinf) < VSMALL)
    {
        FatalErrorIn
        (
            "vector momentCoefficient\n"
            "(\n"
            "    const volVectorField& U,\n"
            "    const volScalarField& p,\n"
            "    const dimensionedScalar& mu,\n"
            "    const word& patchName,\n"
            "    const vector& Uinf,\n"
            "    const scalar& Aref,\n"
            "    const scalar& Lref\n"
            ")\n"
        )   << "Uinf is zero."
            << exit(FatalError);
    }

    const fvMesh& mesh = p.mesh();

    label patchLabel = -1;

    forAll (mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].name() == patchName)
        {
            patchLabel = patchi;
            break;
        }
    }

    vector pressureForceMoment;
    vector viscousForceMoment;

    if (patchLabel != -1)
    {
        pressureForceMoment = sum
        (
            mesh.Cf().boundaryField()[patchLabel]^
            (p.boundaryField()[patchLabel]
            *mesh.Sf().boundaryField()[patchLabel])
        );

        viscousForceMoment = sum
        (
            mesh.Cf().boundaryField()[patchLabel]^
            (-mu.value()*U.boundaryField()[patchLabel].snGrad()
            *mesh.magSf().boundaryField()[patchLabel])
        );
    }
    else
    {
        pressureForceMoment = sum(vectorField(0));
        viscousForceMoment = sum(vectorField(0));
    }

    scalar qRef = 0.5*magSqr(Uinf);
    scalar Fref = qRef*Aref;
    scalar Mref = Fref*Lref;

    vector pressureCoeff = pressureForceMoment/Mref;
    vector viscousCoeff = viscousForceMoment/Mref;

    return (pressureCoeff + viscousCoeff);
}


// ************************************************************************* //
