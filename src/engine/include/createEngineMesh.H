    Info<< "Create mesh\n" << endl;
    autoPtr<engineMesh> meshPtr
    (
        engineMesh::New
        (
            IOobject
            (
                engineMesh::defaultRegion,
                runTime.timeName(),
                runTime
            )
        )
    );

    engineMesh& mesh = meshPtr();
