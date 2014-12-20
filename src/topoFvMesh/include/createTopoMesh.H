    Info<< "Create mesh\n" << endl;
    autoPtr<topoFvMesh> meshPtr
    (
        topoFvMesh::New
        (
            IOobject
            (
                topoFvMesh::defaultRegion,
                runTime.timeName(),
                runTime
            )
        )
    );

    topoFvMesh& mesh = meshPtr();

    mesh.constructAndClear();
