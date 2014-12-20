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
    Pressure-density-based compressible flow solver.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "weighted.H"
#include "gaussConvectionScheme.H"
#include "multivariateGaussConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readThermodynamicProperties.H"
#   include "createFields.H"

    rhoU.correctBoundaryConditions();

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    volScalarField magRhoU = mag(rhoU);
    volScalarField H("H", (rhoE + p)/rho);

    fields.add(rho);
    fields.add(magRhoU);
    fields.add(H);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.value() << nl << endl;

#       include "readPISOControls.H"

        scalar CoNum = max
        (
            mesh.surfaceInterpolation::deltaCoeffs()
           *mag(phiv)/mesh.magSf()
        ).value()*runTime.deltaT().value();

        Info<< "Max Courant Number = " << CoNum << endl;

        magRhoU = mag(rhoU);
        H = (rhoE + p)/rho;

        fv::multivariateGaussConvectionScheme<scalar> mvConvection
        (
            mesh,
            fields,
            phiv,
            mesh.divScheme("div(phiv,rhoUH)")
        );

        for (int outerCorr=0; outerCorr<nOuterCorr; outerCorr++)
        {
            solve
            (
                fvm::ddt(rho)
              + mvConvection.fvmDiv(phiv, rho)
            );

            surfaceScalarField rhoUWeights =
                mvConvection.interpolationScheme()()(magRhoU)()
               .weights(magRhoU);

            weighted<vector> rhoUScheme(rhoUWeights);

            fvVectorMatrix rhoUEqn
            (
                fvm::ddt(rhoU)
              + fv::gaussConvectionScheme<vector>(mesh, phiv, rhoUScheme)
                   .fvmDiv(phiv, rhoU)
            );

            solve(rhoUEqn == -fvc::grad(p));

            solve
            (
                fvm::ddt(rhoE)
              + mvConvection.fvmDiv(phiv, rhoE)
             ==
              - mvConvection.fvcDiv(phiv, p)
            );

            T == (rhoE - 0.5*rho*magSqr(rhoU/rho))/Cv/rho;
            psi = 1.0/(R*T);

            p = rho/psi;

            for (int corr=0; corr<nCorr; corr++)
            {
                volScalarField rrhoUA = 1.0/rhoUEqn.A();
                surfaceScalarField rrhoUAf = fvc::interpolate(rrhoUA);

                volVectorField HbyA = rrhoUA*rhoUEqn.H();

                phi = (fvc::interpolate(HbyA) & mesh.Sf())
                    + fvc::ddtPhiCorr(rrhoUA, rho, rhoU, phi);

                p.boundaryField().updateCoeffs();

                surfaceScalarField phiGradp =
                    rrhoUAf*mesh.magSf()*fvc::snGrad(p);

                phi -= phiGradp;

                phi.boundaryField() =
                    rhoU.boundaryField() & mesh.Sf().boundaryField();

                phiv = phi/fvc::interpolate(rho);

                fvScalarMatrix pEqn
                (
                    fvm::ddt(psi, p)
                  + mvConvection.fvcDiv(phiv, rho)
                  + fvc::div(phiGradp)
                  - fvm::laplacian(rrhoUAf, p)
                );

                pEqn.solve();

                phi += phiGradp + pEqn.flux();
                phiv = phi/fvc::interpolate(rho);
                rho = psi*p;

                rhoU = HbyA - rrhoUA*fvc::grad(p);
                rhoU.correctBoundaryConditions();
            }
        }

        U == rhoU/rho;

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
