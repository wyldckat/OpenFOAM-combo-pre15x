volScalarField A = UEqn.A();

surfaceScalarField rUAf
(
    "(rho|A(U))",
    fvc::interpolate(rho)*fvc::interpolate(1/A)
);

U = UEqn.H()/A;
phi = fvc::interpolate(rho)*
(
    (
        fvc::interpolate
        (
            U + UphiCoeff*fvc::ddt0(rho, U)/A, "interpolate((H(U)|A(U)))"
        ) & mesh.Sf()
    ) - UphiCoeff*fvc::interpolate(1.0/A)*fvc::ddt0(phi)
);

surfaceScalarField phiU("phiU", phi);

phi -= ghf*fvc::snGrad(rho)*rUAf*mesh.magSf();


for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rUAf, p) == fvc::ddt(rho) + fvc::div(phi)
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

#include "rhoEqn.H"
#include "continuityErrs.H"

U += fvc::reconstruct((phi - phiU)/fvc::interpolate(rho));

U.correctBoundaryConditions();
