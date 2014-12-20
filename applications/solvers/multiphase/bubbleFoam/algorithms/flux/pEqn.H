{
    surfaceScalarField alphaf = fvc::interpolate(alpha);
    surfaceScalarField betaf = 1.0 - alphaf;

    Ua = UaEqn.H()/UaEqn.A();
    Ub = UbEqn.H()/UbEqn.A();

    surfaceScalarField phiHa = (fvc::interpolate(Ua) & mesh.Sf());
    surfaceScalarField phiHb = (fvc::interpolate(Ub) & mesh.Sf());

    phia = phiHa
      + fvc::interpolate(beta/rhoa*dragCoef/UaEqn.A())*phib
      + (g*fvc::interpolate(1.0/UaEqn.A()) & mesh.Sf());

    phib = phiHb
      + fvc::interpolate(alpha/rhob*dragCoef/UbEqn.A())*phia
      + (g*fvc::interpolate(1.0/UbEqn.A()) & mesh.Sf());

    phi = alphaf*phia + betaf*phib;

    surfaceScalarField Dpa = fvc::interpolate((1/rhoa)/UaEqn.A());
    surfaceScalarField Dpb = fvc::interpolate((1/rhob)/UbEqn.A());
    surfaceScalarField Dp("(rho|A(U))", alphaf*Dpa + betaf*Dpb);

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phi)
        );

        fvScalarMatrix::reference pRef =
            pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();
        pEqn.unsetReference(pRef);

        if (nonOrth == nNonOrthCorr)
        {
            surfaceScalarField SfGradp = pEqn.flux()/Dp;

            phia -= Dpa*SfGradp;
            phib -= Dpb*SfGradp;
            phi = alphaf*phia + betaf*phib;

            p.relax();

            SfGradp -= pEqn.flux()/Dp;

            Ua += fvc::reconstruct(phia - phiHa + Dpa*SfGradp);
            //Ua = fvc::reconstruct(phia);
            Ua.correctBoundaryConditions();

            Ub += fvc::reconstruct(phib - phiHb + Dpb*SfGradp);
            //Ub = fvc::reconstruct(phib);
            Ub.correctBoundaryConditions();

            U = alpha*Ua + beta*Ub;
        }
    }
}

#include "continuityErrs.H"
