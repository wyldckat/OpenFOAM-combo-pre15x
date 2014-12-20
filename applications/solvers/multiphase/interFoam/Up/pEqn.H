{
    volScalarField A = UEqn.A();

    surfaceScalarField rUAf
    (
        "1|A(U)",
        fvc::interpolate(1.0/A)/fvc::interpolate(rho)
    );

    {
        U = UEqn.H()/A;

        // Create absolute predicted fluxes
        phi = 
        (
            fvc::interpolate
            (
                U + UphiCoeff*fvc::ddt0(U)/A, "interpolate((H(U)|A(U)))"
            ) & mesh.Sf()
        ) - UphiCoeff*fvc::interpolate(1.0/A)
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
