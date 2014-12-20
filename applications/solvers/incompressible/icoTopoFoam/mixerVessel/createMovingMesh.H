    IOdictionary mixerProperties
    (
        IOobject
        (
            "mixerProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    mixerMesh mesh
    (
        IOobject
        (
            mixerMesh::defaultRegion,
            runTime.timeName(),
            runTime
        ),
        mixerProperties
    );
