    Info<< "Create mesh\n" << endl;
    autoPtr<dynamicFvMesh> meshPtr
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                dynamicFvMesh::defaultRegion,
                runTime.timeName(),
                runTime
            )
        )
    );

    dynamicFvMesh& mesh = meshPtr();
