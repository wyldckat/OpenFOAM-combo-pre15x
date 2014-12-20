{
    word gammarScheme("div(phirb,gamma)");

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(interface.cGamma()*phic, max(phic));
    surfaceScalarField phir = phic*interface.nHatf();

    for (int gCorr=0; gCorr<nGammaCorr; gCorr++)
    {
        fvScalarMatrix gammaEqn
        (
            fvm::ddt(gamma)
          + fvm::div(phi, gamma)
          + fvm::div
            (
                -fvc::flux(-phir, scalar(1) - gamma, gammarScheme),
                gamma,
                gammarScheme
            )
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
