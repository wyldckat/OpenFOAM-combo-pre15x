{
    surfaceScalarField alphaf = fvc::interpolate(alpha);
    surfaceScalarField betaf = 1.0 - alphaf;

    Ua = UaEqn.H()/UaEqn.A();
    Ub = UbEqn.H()/UbEqn.A();

    phia = (fvc::interpolate(Ua) & mesh.Sf());
    phib = (fvc::interpolate(Ub) & mesh.Sf());
    phi = alphaf*phia + betaf*phib;

    surfaceScalarField Dpa = fvc::interpolate((1.0/rhoa)/UaEqn.A());
    surfaceScalarField Dpb = fvc::interpolate((1.0/rhob)/UbEqn.A());
    surfaceScalarField Dp("(rho|A(U))", alphaf*Dpa + betaf*Dpb);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
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

            volVectorField gradp = fvc::grad(p);

            Ua -= gradp/(rhoa*UaEqn.A());
            Ua.correctBoundaryConditions();

            Ub -= gradp/(rhob*UbEqn.A());
            Ub.correctBoundaryConditions();

            U = alpha*Ua + beta*Ub;
        }
    }
}

#include "continuityErrs.H"
