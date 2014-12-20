{
    IOobject ioPoints
    (
        "points", 
        runTime.timeName(),
        mesh.name(),
        mesh
    );

    if (ioPoints.headerOk())
    {
        // Reading new points
        pointIOField newPoints
        (
            IOobject
            (
                "points",
                mesh.time().timeName(),
                mesh.name(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        
        mesh.movePoints(newPoints);
    }
}
