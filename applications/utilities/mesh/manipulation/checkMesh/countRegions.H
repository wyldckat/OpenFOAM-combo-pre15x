    regionSplit rs(mesh);

    if (rs.nRegions() == 1)
    {
        Info<< "Number of regions: " << rs.nRegions() << " (OK)."<< endl;
        
    }
    else
    {
        Info<< "Number of regions: " << rs.nRegions() << endl;

        Info<< "The mesh has multiple regions which are not connected by "
            << "any face." << endl
            << "Dumping region information to "
            << runTime.timeName() << "/cellToRegion" << nl
            << "This can be used by e.g. the manual"
            << " decomposition option in decomposePar."
            << endl;

        labelIOList ctr
        (
            IOobject
            (
                "cellToRegion",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rs.cellToRegion()
        );
        ctr.write();
    }

