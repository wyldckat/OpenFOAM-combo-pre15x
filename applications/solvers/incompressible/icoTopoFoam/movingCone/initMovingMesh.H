    // Reading mesh modifiers dictionary
    Info << "Reading mesh motion control" << endl;
    IOdictionary meshMotionControl
    (
        IOobject
        (
            "meshMotionControl",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    vector motionU(meshMotionControl.lookup("motionU"));

    Info<< "Motion velocity: " << motionU << endl;

    // Vertex markup is used to facilitate the motion setup.
    // For all vertices left of the moving body, the markup is set to -1
    // For the vertices above the body, the markup is set to zero and for
    // the vertices to the right it is set to 1
    scalar leftEdge = -0.0075;
    scalar edgeTol = 1e-9;

    scalar curLeft = 0;
    scalar curRight = 0;

    scalarField motionMarkup(mesh.nPoints(), 0);

#   include "resetLeftRight.H"
