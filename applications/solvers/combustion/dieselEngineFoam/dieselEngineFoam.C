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
    dieselFoam

Description
    Diesel engine spray and combustion code.

\*---------------------------------------------------------------------------*/

#include "engineTime.H"
#include "fvCFD.H"
#include "tetFem.H"
#include "hCombustionThermo.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "spray.H"
#include "chemistryModel.H"
#include "chemistrySolver.H"

#include "multivariateScheme.H"
#include "fixedValueTetPolyPatchFields.H"
#include "slipTetPolyPatchFields.H"
#include "symmetryFvPatch.H"
#include "wedgeFvPatch.H"
#include "emptyFvPatch.H"
#include "Switch.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createEngineTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "readEnvironmentalProperties.H"
#   include "readCombustionProperties.H"
#   include "createSpray.H"
#   include "initContinuityErrs.H"
#   include "createEngineMovingMesh.H"
#   include "readEngineTimeControls.H"
#   include "setInitialDeltaT.H"
#   include "startSummary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readPISOControls.H"
#       include "readEngineTimeControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "\nCrank angle = " << runTime.theta() << " CA-deg"
            << endl;

#       include "movePiston.H"

        dieselSpray.evolve();

        Info << "Solving chemistry" << endl;

        chemistry.solve
        (
            runTime.value() - runTime.deltaT().value(),
            runTime.deltaT().value()
        );

        // turbulent time scale
        {
            volScalarField tk =
                Cmix*sqrt(turbulence->muEff()/rho/turbulence->epsilon());
            volScalarField tc = chemistry.tc();

            //Chalmers PaSR model
            kappa = (runTime.deltaT() + tc)/(runTime.deltaT() + tc + tk);
        }

#       include "rhoEqn.H"
#       include "UEqn.H"

        for (label ocorr=1; ocorr <= nOuterCorr; ocorr++)
        {
#           include "YEqn.H"
#           include "hEqn.H"
            
            // --- PISO loop
            for (int corr=1; corr<=nCorr; corr++)
            {
#               include "pEqn.H"
            }
        }

        turbulence->correct();

#       include "logSummary.H"
#       include "spraySummary.H"

        runTime.write();

        Info<< "\nExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl;
    }

    Info << "\n end\n";

    return(0);
}


// ************************************************************************* //
