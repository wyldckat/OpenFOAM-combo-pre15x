    tetPolyMesh tetMesh(mesh);

    Info<< "Creating field motionDiff\n" << endl;
    elementScalarField motionDiff
    (
        IOobject
        (
            "motionDiff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tetMesh,
        dimensionedScalar("1.0", dimless, 1.0)
    );
