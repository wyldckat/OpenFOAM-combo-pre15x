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
    dnsFoam

Description
    Direct numerical simulation solver for boxes of isotropic turbulence

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Kmesh.H"
#include "UOprocess.H"
#include "fft.H"
#include "calcEk.H"
#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMeshNoClear.H"
#   include "readTransportProperties.H"
#   include "createFields.H"
#   include "readTurbulenceProperties.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << endl;

#       include "readPISOControls.H"

        force.internalField() = ReImSum
        (
            fft::reverseTransform
            (
                K/(mag(K) + 1.0e-6) ^ forceGen.newField(), K.nn()
            )
        );

#       include "globalProperties.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(U) 
          + fvm::div(phi, U) 
          - fvm::laplacian(nu, U)
         ==
            force
        );

        solve(UEqn == -fvc::grad(p));


        // --- PISO loop

        for (int corr=1; corr<=1; corr++)
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

            fvScalarMatrix pEqn
            (
                fvm::laplacian(1.0/A, p) == fvc::div(phi)
            );

            pEqn.solve();

            phi -= pEqn.flux();

#           include "continuityErrs.H"

            U -= fvc::grad(p)/A;
            U.correctBoundaryConditions();
        }

        runTime.write();

        if (runTime.outputTime())
        {
            calcEk(U, K).write(runTime.timePath()/"Ek", runTime.graphFormat());
        }

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl << endl;
    }

    Info<< "end" << endl;

    return(0);
}


// ************************************************************************* //
