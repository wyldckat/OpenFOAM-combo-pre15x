    Info<< "Create mesh\n" << endl;
    autoPtr<movingFvMesh> meshPtr
    (
        movingFvMesh::New
        (
            IOobject
            (
                movingFvMesh::defaultRegion,
                runTime.timeName(),
                runTime
            )
        )
    );

    movingFvMesh& mesh = meshPtr();

    mesh.constructAndClear();
