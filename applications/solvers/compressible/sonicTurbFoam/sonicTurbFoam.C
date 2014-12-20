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
    sonicFoam

Description
    Transient solver for trans-sonic/supersonic, turbulent flow of a 
    compressible gas.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicThermo.H"
#include "compressible/turbulenceModel/turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "\nTime = " << runTime.timeName() << nl << endl;

#       include "readPISOControls.H"
#       include "CourantNo.H"

#       include "rhoEqn.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(rho, U)
          + fvm::div(phi, U)
          + turbulence->divRhoR(U)
        );

        solve(UEqn == -fvc::grad(p));

#       include "hEqn.H"


        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            volScalarField A = UEqn.A();
            U = UEqn.H()/A;

            surfaceScalarField phid =
            (
                fvc::interpolate
                (
                    psi*(U + UphiCoeff*fvc::ddt0(rho, U)/A),
                    "interpolate((H(U)|A(U)))"
                ) & mesh.Sf()
            ) - UphiCoeff*fvc::interpolate(rho/A)*fvc::ddt0(phi)
               /fvc::interpolate(p);

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::ddt(psi, p)
                  + fvm::div(phid, p, "div(phid,p)")
                  - fvm::laplacian(rho/A, p)
                );

                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi = pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            U -= fvc::grad(p)/A;
            U.correctBoundaryConditions();
        }

        DpDt = 
            fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

        turbulence->correct();

        runTime.write();

        Info<< "\n    ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
