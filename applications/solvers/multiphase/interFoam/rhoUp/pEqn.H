{
    volScalarField A = UEqn.A();

    surfaceScalarField rUAf("(rho|A(U))", fvc::interpolate(1.0/A));

    {
        U = UEqn.H()/A;

        // U-phi coupling algorithm flux predictor.
        // Note: this version will probabbly only work using Euler time scheme.
        phi = 
        (
            fvc::interpolate
            (
                U + UphiCoeff*fvc::ddt0(rho, U)/A, "interpolate((H(U)|A(U)))"
            ) & mesh.Sf()
        ) - UphiCoeff*fvc::interpolate(rho.oldTime()/A)
           *(fvc::ddt0(phi) + fvc::ddt0(mesh.phi()));
    }

    surfaceScalarField phiU("phiU", phi);

    phi +=
    (
        fvc::interpolate(interface.sigmaK())*fvc::snGrad(gamma)
      - ghf*fvc::snGrad(rho)
    )*rUAf*mesh.magSf();

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pd) == fvc::div(phi)
        );

        fvScalarMatrix::reference pRef =
            pdEqn.setReference(pRefCell, pRefValue);
        pdEqn.solve();
        pdEqn.unsetReference(pRef);

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pdEqn.flux();
        }
    }

    U += fvc::reconstruct(phi - phiU);
    U.correctBoundaryConditions();

    // Make the fluxes relative
    phi -= mesh.phi();
}
