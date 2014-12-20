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
    refinePolyMesh

Description
    Refine a polyhedral mesh

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "vertexMarkup.H"
#include "movingMesh.H"

#include "movingWallVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMovingMesh.H"
#   include "initMovingMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readPISOControls.H"
#       include "readTimeControls.H"
        runTime++;

        Info<< "Time = " << runTime.value() << nl << endl;

#       include "moveMesh.H"
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
            phi -= mesh.phi();
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

            phi -= mesh.phi();
        }

#       include "CourantNo.H"
#       include "setDeltaT.H"

        // To be able to evaluate streamFunction

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
