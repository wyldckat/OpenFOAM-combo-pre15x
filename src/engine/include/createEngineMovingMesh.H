    label pistonIndex = -1;
    bool foundPiston = false;

    label linerIndex = -1;
    bool foundLiner = false;

    label cylinderHeadIndex = -1;
    bool foundCylinderHead = false;

    forAll(mesh.boundary(), i)
    {
        if (mesh.boundary()[i].name() == "piston")
        {
            pistonIndex = i;
            foundPiston = true;
        }
        else if (mesh.boundary()[i].name() == "liner")
        {
            linerIndex = i;
            foundLiner = true;
        }
        else if (mesh.boundary()[i].name() == "cylinderHead")
        {
            cylinderHeadIndex = i;
            foundCylinderHead = true;
        }
    }

    reduce(foundPiston, orOp<bool>());
    reduce(foundLiner, orOp<bool>());
    reduce(foundCylinderHead, orOp<bool>());

    if (!foundPiston)
    {
        FatalErrorIn("main")
            << args.executable()
            << " : cannot find piston patch"
            << exit(FatalError);
    }

    if (!foundLiner)
    { 
        FatalErrorIn("main")
            << args.executable()
            << " : cannot find liner patch"
            << exit(FatalError);
    }

    if (!foundCylinderHead)
    { 
        FatalErrorIn("main")
            << args.executable()
            << " : cannot find cylinderHead patch"
            << exit(FatalError);
    }

    dimensionedScalar deckHeight("deckHeight", dimLength, GREAT);
    dimensionedScalar pistonPosition("deckHeight", dimLength, GREAT);

    {
        if (linerIndex != -1)
        {
            //pistonPosition = 
            //    min(mesh.boundary()[linerIndex].patch().localPoints()).z();
            pistonPosition.value() =
                max(mesh.boundary()[pistonIndex].patch().localPoints()).z();
        }
        reduce(pistonPosition.value(), minOp<scalar>());

        if (cylinderHeadIndex != -1)
        {
            deckHeight.value() = min
            (
                mesh.boundary()[cylinderHeadIndex].patch().localPoints()
            ).z();
        }
        reduce(deckHeight.value(), minOp<scalar>());

        Info<< "deckHeight: " << deckHeight.value() << nl
            << "piston position: " << pistonPosition.value() << endl;
    }

    Switch layeredMesh(false);

    if (runTime.engineDict().found("layeredMesh"))
    {
        runTime.engineDict().lookup("layeredMesh") >> layeredMesh;
    }

    dimensionedScalar pistonLayers("pistonLayers", dimLength, 0.0);
    if (runTime.engineDict().found("pistonLayers"))
    {
        runTime.engineDict().lookup("pistonLayers") >> pistonLayers;
    }
