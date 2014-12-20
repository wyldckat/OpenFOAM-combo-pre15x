    surfaceScalarField gammaf = fvc::interpolate(gamma);
    surfaceScalarField rhof = gammaf*rho1 + (1.0 - gammaf)*rho2;

    surfaceScalarField muf
    (
        "muEff",
        rhof*
        (
            gammaf*fvc::interpolate(phase1->nu())
          + (1.0 - gammaf)*fvc::interpolate(phase2->nu())
        )
    );

    surfaceScalarField UBlendingFactor
    (
        IOobject
        (
            "UBlendingFactor",
            runTime.timeName(),
            mesh
        ),
        max
        (
            1.0 - mag(1.0 - Foam::pow(4.0*(gammaf - 0.5)*(0.5 - gammaf), 4.0)),
            0.0
        )
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muf, U)
      - (fvc::grad(U) & fvc::grad(muf))
    );
