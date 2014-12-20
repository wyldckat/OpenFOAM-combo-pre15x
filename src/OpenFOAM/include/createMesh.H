    Info<< "Create mesh for time = "
        << runTime.value() << nl << endl;

    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    mesh.constructAndClear();
