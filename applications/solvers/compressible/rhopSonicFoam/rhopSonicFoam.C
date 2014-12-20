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
    rhopSonicFoam

Description
    Pressure-based compressible flow solver using density-weighted variables.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readThermodynamicProperties.H"
#   include "createFields.H"
#   include "readPISOControls.H"
#   include "initContinuityErrs.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    cpuTime executionTime;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "\n Time = " << runTime.value() << nl << endl;

        surfaceScalarField phiv
        (
            IOobject
            (
                "phiv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phi/fvc::interpolate(rho)
        );

        scalar CoNum = max
        (
            mesh.surfaceInterpolation::deltaCoeffs()
           *mag(phiv)/mesh.magSf()
        ).value()*runTime.deltaT().value();

        Info<< "\nMax Courant Number = " << CoNum << endl;

#       include "rhoEqn.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(rhoU)
          + fvm::div(phiv, rhoU)
        );

        solve(UEqn == -fvc::grad(p));

        solve
        (
            fvm::ddt(rhoE)
          + fvm::div(phiv, rhoE)
         == 
          - fvc::div(phiv, p)
        );

        T = (rhoE - 0.5*rho*magSqr(rhoU/rho))/Cv/rho;
        psi = 1.0/(R*T);

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            rhoU = UEqn.H()/UEqn.A();

            surfaceScalarField phid =
                (fvc::interpolate(rhoU)/fvc::interpolate(p) & mesh.Sf());

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::ddt(psi, p)
                  + fvm::div(phid, p, "div(phid,p)")
                  - fvm::laplacian(1.0/UEqn.A(), p)
                );

                pEqn.solve();

                phi = pEqn.flux();
            }

#           include "continuityErrs.H"

            rhoU -= fvc::grad(p)/UEqn.A();
            rhoU.correctBoundaryConditions();
        }

        runTime.write();

        Info<< "\n    ExecutionTime = "
            << executionTime.elapsedCpuTime()
            << " s\n" << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
