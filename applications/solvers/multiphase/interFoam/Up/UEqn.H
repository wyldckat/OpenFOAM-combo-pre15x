    surfaceScalarField gammaf = fvc::interpolate(gamma);

    surfaceScalarField nuf
    (
        "nuEff",
        gammaf*fvc::interpolate(phase1->nu())
      + (1.0 - gammaf)*fvc::interpolate(phase2->nu())
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
        fvm::ddt(U)
      + fvm::div(phi, U)
      - fvm::laplacian(nuf, U)
      - (fvc::grad(U) & fvc::grad(nuf))
    );
