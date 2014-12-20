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
    sonicFoamAutoMotion

Description
    Transient solver for trans-sonic/supersonic, laminar flow of a 
    compressible gas with mesh motion..

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readThermodynamicProperties.H"
#   include "readTransportProperties.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    motionSolver ms(mesh);

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "\n Time = " << runTime.timeName() << nl << endl;

#       include "readPISOControls.H"
#       include "CourantNo.H"

        mesh.movePoints(ms.newPoints());

#       include "rhoEqn.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(rho, U)
          + fvm::div(phi, U)
          - fvm::laplacian(mu, U)
        );

        solve(UEqn == -fvc::grad(p));

        solve
        (
            fvm::ddt(rho, e)
          + fvm::div(phi, e)
          - fvm::laplacian(mu, e)
          ==
          - p*fvc::div(phi/fvc::interpolate(rho) + mesh.phi())
          + mu*magSqr(symm(fvc::grad(U)))
        );

        T = e/Cv;
        psi = 1.0/(R*T);

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            U = UEqn.H()/UEqn.A();

            surfaceScalarField phid =
                fvc::interpolate(psi)*
                (
                    (fvc::interpolate(U) & mesh.Sf()) - mesh.phi()
                );

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::ddt(psi, p)
                  + fvm::div(phid, p, "div(phid,p)")
                  - fvm::laplacian(rho/UEqn.A(), p)
                );

                pEqn.solve();

                phi = pEqn.flux();
            }

#           include "continuityErrs.H"

            U -= fvc::grad(p)/UEqn.A();
            U.correctBoundaryConditions();

            rho = psi*p;
        }

        runTime.write();

        Info<< "\n    ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl;
    }

    Info<< "\nEnd" << endl;

    return(0);
}


// ************************************************************************* //
