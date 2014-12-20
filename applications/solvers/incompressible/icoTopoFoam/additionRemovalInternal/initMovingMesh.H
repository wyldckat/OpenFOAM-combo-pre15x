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
    scalar motionRefX(readScalar(meshMotionControl.lookup("motionRefX")));
