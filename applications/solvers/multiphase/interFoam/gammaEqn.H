{
    surfaceScalarField phir
    (
        IOobject
        (
            "phir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        interface.cGamma()*mag(phi/mesh.magSf())*interface.nHatf()
    );

    for (int gCorr=0; gCorr<nGammaCorr; gCorr++)
    {
        word gammarScheme("div(phirb,gamma)");
        surfaceScalarField phirb = -fvc::flux(-phir, 1.0 - gamma, gammarScheme);

        fvScalarMatrix gammaEqn
        (
            fvm::ddt(gamma)
          + fvm::div(phi, gamma)
          + fvm::div(phirb, gamma, gammarScheme)
        );

        gammaEqn.solve();

        rhoPhi = gammaEqn.flux()*(rho1 - rho2) + phi*rho2;
    }

    Info<< "Liquid phase volume fraction = "
        << gamma.weightedAverage(mesh.V()).value()
        << "  Min(gamma) = " << min(gamma).value()
        << "  Max(gamma) = " << max(gamma).value()
        << endl;
}
