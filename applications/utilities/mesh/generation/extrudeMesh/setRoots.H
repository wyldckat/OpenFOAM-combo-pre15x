    argList::validArgs.clear();

    argList::validArgs.append("source root");
    argList::validArgs.append("source case");

    argList::validArgs.append("source patch");

    argList::validArgs.append("target root");
    argList::validArgs.append("target case");

    argList::validArgs.append("nLayers");
    argList::validArgs.append("overal thickness");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    fileName rootDirSource(args.args()[1]);
    fileName caseDirSource(args.args()[2]);

    word patchName(args.args()[3]);

    fileName rootDirTarget(args.args()[4]);
    fileName caseDirTarget(args.args()[5]);

    label nLayers(readLabel(IStringStream(args.args()[6])()));
    scalar thickness(readScalar(IStringStream(args.args()[7])()));

    Info<< "Source: " << rootDirSource << " " << caseDirSource << nl
        << "Target: " << rootDirTarget << " " << caseDirTarget << endl;

    Info<< "Extruding layers:" << nl
        << "    template source patch " << patchName << nl
        << "    number of layers " << nLayers << nl
        << "    overall thickness " << thickness << nl
        << endl;
