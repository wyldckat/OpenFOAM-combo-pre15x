    curLeft =
        average
        (
            mesh.faceZones()
            [
                mesh.faceZones().findZoneID("leftExtrusionFaces")
            ]().localPoints()
        ).x() - edgeTol;

    curRight =
        average
        (
            mesh.faceZones()
            [
                mesh.faceZones().findZoneID("rightExtrusionFaces")
            ]().localPoints()
        ).x() + edgeTol;

    Info << "Updating vertex markup" << endl;

    motionMarkup =
    (
        vertexMarkup
        (
            mesh.points(),
            curLeft,
            curRight
        )
    );
