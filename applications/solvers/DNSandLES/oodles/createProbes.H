    Info<< "Reading probeLocations\n" << endl;

    IOdictionary probeLocations
    (
        IOobject
        (
            "probeLocations",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Info<< "Constructing probes\n" << endl;

    Probe<volVectorField> UProbe(probeLocations, U);
    Probe<volScalarField> pProbe(probeLocations, p);
