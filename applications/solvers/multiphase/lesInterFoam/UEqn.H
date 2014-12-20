    surfaceScalarField gammaf = fvc::interpolate(gamma);

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
      + turbulence->divB(U)
    );
