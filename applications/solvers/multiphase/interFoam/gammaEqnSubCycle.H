if (nGammaSubCycles > 1)
{
    surfaceScalarField rhoPhiSum = 0.0*rhoPhi;
    dimensionedScalar totalDeltaT = runTime.deltaT();

    for
    (
        subCycle<volScalarField> gammaSubCycle(gamma, nGammaSubCycles);
        !(++gammaSubCycle).end();
    )
    {
#       include "gammaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
#       include "gammaEqn.H"
}

interface.correct();

rho == gamma*rho1 + (1.0 - gamma)*rho2;
