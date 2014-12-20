    Info<< "\nReading transportProperties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );


    //  Read centerline velocity for channel simulations
    dimensionedVector Ubar
    (
        transportProperties.lookup("Ubar")
    );

    vector flowDirection = (Ubar/mag(Ubar)).value();
    tensor flowMask = sqr(flowDirection);
