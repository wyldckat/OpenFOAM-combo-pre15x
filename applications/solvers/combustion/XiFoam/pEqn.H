rho = thermo->rho();

volScalarField A = UEqn.A();
U = UEqn.H()/A;

if (transSonic)
{
    surfaceScalarField phid = 
    (
        fvc::interpolate
        (
            psi*(U + UphiCoeff*fvc::ddt0(rho, U)/A),
            "interpolate((H(U)|A(U)))"
        ) & mesh.Sf()
    ) - UphiCoeff*fvc::interpolate(rho/A)*fvc::ddt0(phi)/fvc::interpolate(p);

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
            phi == pEqn.flux();
        }
    }
}
else
{
    phi =
    (
        fvc::interpolate
        (
            rho*(U + UphiCoeff*fvc::ddt0(rho, U)/A), "interpolate((H(U)|A(U)))"
        ) & mesh.Sf()
    ) - UphiCoeff*fvc::interpolate(rho/A)*fvc::ddt0(phi);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho/A, p)
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "continuityErrs.H"

U -= fvc::grad(p)/A;
U.correctBoundaryConditions();

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
