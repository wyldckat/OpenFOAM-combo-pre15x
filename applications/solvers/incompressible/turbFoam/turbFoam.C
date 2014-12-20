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
    turbFoam

Description
    Transient solver for incompressible, turbulent flow.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\n Starting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "\nTime = " << runTime.timeName() << nl << endl;

#       include "readPISOControls.H"
#       include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            // Momentum predictor

            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi, U)
              + turbulence->divR(U)
            );

            solve(UEqn == -fvc::grad(p));

            // --- PISO loop

            for (int corr=0; corr<nCorr; corr++)
            {
                volScalarField A = UEqn.A();

                U = UEqn.H()/A;
                phi = 
                (
                    fvc::interpolate
                    (
                        U + UphiCoeff*fvc::ddt0(U)/A, "interpolate((H(U)|A(U)))"
                    ) & mesh.Sf()
                ) - UphiCoeff*fvc::interpolate(1.0/A)*fvc::ddt0(phi);

                // Non-orthogonal pressure corrector loop
                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    // Pressure corrector

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

#               include "continuityErrs.H"

                U -= fvc::grad(p)/A;
                U.correctBoundaryConditions();
            }
        }

        turbulence->correct();

        runTime.write();

        Info<< "\n    ExecutionTime = "
            << runTime.elapsedCpuTime() << " s\n" << endl;
    }

    Info<< "End" << endl;

    return(0);
}


// ************************************************************************* //
