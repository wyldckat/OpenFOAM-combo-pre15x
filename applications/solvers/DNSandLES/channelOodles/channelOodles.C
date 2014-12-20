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
    oodles

Description
    Incompressible LES solver for flow in a channel.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/transportModel/transportModel.H"
#include "incompressible/LESmodel/LESmodel.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Random.H"
#include "Probe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readTransportProperties.H"
#   include "createFields.H"
#   include "createAverages.H"
#   include "createProbes.H"
#   include "initContinuityErrs.H"
#   include "createGradP.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for(runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readPISOControls.H"

#       include "CourantNo.H"

        sgsModel->correct();

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          + sgsModel->divB(U)
         ==
            gradP
        );

        if (momentumPredictor)
        {
            solve(UEqn == -fvc::grad(p));
        }


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

            U -= fvc::grad(p)/A;
            U.correctBoundaryConditions();
        }


        // Correct driving force for a constant mass flow rate

        dimensionedVector UbarStar = flowMask & U.weightedAverage(mesh.V());

        U += (Ubar - UbarStar);
        gradP += (Ubar - UbarStar)/(1.0/UEqn.A())().weightedAverage(mesh.V());

        Info<< "Uncorrected Ubar = " << (flowDirection & UbarStar.value())<< tab
            << "pressure gradient = " << (flowDirection & gradP.value())<< endl;


#       include "calculateAverages.H"

        runTime.write();

#       include "writeNaveragingSteps.H"

#       include "writeGradP.H"

#       include "writeProbes.H"

        Info<< "ExecutionTime = "
             << runTime.elapsedCpuTime()
             << " s\n" << endl;
    }

    Info<< "Final Execution Time = "
         << runTime.elapsedCpuTime()
         << " s\n" << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
