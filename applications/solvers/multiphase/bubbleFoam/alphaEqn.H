{
    word scheme("div(phi,alpha)");

    surfaceScalarField phir = phia - phib;

    Info<< "Max Ur Courant Number = "
        << (
               max
               (
                   mesh.surfaceInterpolation::deltaCoeffs()*mag(phir)
                  /mesh.magSf()
                )*runTime.deltaT()
            ).value()
        << endl;

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        fvScalarMatrix alphaEqn
        (
             fvm::ddt(alpha)
           + fvm::div(phi, alpha, scheme)
           + fvm::div(-fvc::flux(-phir, beta, scheme), alpha, scheme)
        );
        alphaEqn.relax();
        alphaEqn.solve();

        /*
        fvScalarMatrix betaEqn
        (
              fvm::ddt(beta)
            + fvm::div(phi, beta, scheme)
            + fvm::div(-fvc::flux(phir, 1.0 - beta, scheme), beta, scheme)
        );
        betaEqn.relax();
        betaEqn.solve();

        alpha = 0.5*(1.0 + sqr(1.0 - beta) - sqr(1.0 - alpha));
        */

        beta = 1.0 - alpha;
    }

    Info<< "Dispersed phase volume fraction = "
        << alpha.weightedAverage(mesh.V()).value()
        << "  Min(alpha) = " << min(alpha).value()
        << "  Max(alpha) = " << max(alpha).value()
        << endl;
}
