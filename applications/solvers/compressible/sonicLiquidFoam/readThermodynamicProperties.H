    Info<< "Reading thermodynamicProperties\n" << endl;

    IOdictionary thermodynamicProperties
    (
        IOobject
        (
            "thermodynamicProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rho0
    (
        thermodynamicProperties.lookup("rho0")
    );

    dimensionedScalar p0
    (
        thermodynamicProperties.lookup("p0")
    );

    dimensionedScalar psi
    (
        thermodynamicProperties.lookup("psi")
    );
