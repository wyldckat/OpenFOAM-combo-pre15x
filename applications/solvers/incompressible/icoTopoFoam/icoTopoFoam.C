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

Application
    icoTopoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids
    with mesh motion and topological mesh changes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mixerFvMesh.H"

#include "movingWallVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createTopoMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readPISOControls.H"
#       include "readTimeControls.H"
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.moveAndMorph();

#       include "checkTotalVolume.H"

        A.correctBoundaryConditions();
        H.correctBoundaryConditions();

        surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh
            ),
            fvc::interpolate(H/A, "interpolate(HbyA)") & mesh.Sf()
        );

        {
            volScalarField p2 = p;

            fvScalarMatrix pEqn
            (
                fvm::laplacian(1.0/A, p2) == fvc::div(phi)
            );

            if (p.needReference())
            {
                pEqn.source() *= totalVolRatio;
            }

            fvScalarMatrix::reference pRef =
                pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();
            pEqn.unsetReference(pRef);

            phi -= pEqn.flux();
            phi -= fvc::meshPhi(U);
        }

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::Sp(fvc::div(phi), U)
          - fvm::laplacian(nu, U)
        );

        //solve(UEqn == -fvc::grad(p));

        A = UEqn.A();

        solve
        (
            UEqn
          ==
           -A
           *fvc::reconstruct
            (
                fvc::interpolate(1.0/A)*mesh.magSf()*fvc::snGrad(p)
            )
        );


        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            A = UEqn.A();
            H = UEqn.H();

            U = H/A;
            phi = fvc::interpolate(U, "interpolate(HbyA)") & mesh.Sf();
            surfaceScalarField phiU("phiU", phi);

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(1.0/A, p) == fvc::div(phi)
                );

                fvScalarMatrix::reference pRef =
                    pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();
                pEqn.unsetReference(pRef);
                
                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            //U -= fvc::grad(p)/A;
            U += fvc::reconstruct(phi - phiU);
            U.correctBoundaryConditions();

            phiField = phi;

            // Make the fluxes relative
            phi -= fvc::meshPhi(U);
        }

#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
