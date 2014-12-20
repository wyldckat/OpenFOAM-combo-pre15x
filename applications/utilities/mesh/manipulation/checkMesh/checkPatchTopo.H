    if (!Pstream::parRun())
    {
        Info<< "\nChecking patch topology for multiply connected surfaces ..."
            << endl;

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // Non-manifold points
        pointSet points
        (
            mesh,
            "nonManifoldPoints",
            mesh.nPoints()/100
        );

        Info << "\nPatch\tFaces\tPoints\tSurface" << endl;
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];


            primitivePatch::surfaceTopo pTyp = pp.surfaceType();

            if (pTyp == primitivePatch::MANIFOLD)
            {
                if (pp.checkPointManifold(true, &points))
                {
                    Info<< pp.name() << '\t' << pp.size()
                        << '\t' << pp.nPoints()
                        << "\tmultiply connected (shared point)" << endl;
                }
                else
                {
                    Info<< pp.name() << '\t' << pp.size()
                        << '\t' << pp.nPoints()
                        << "\tok (closed singly connected surface)" << endl;
                }
            }
            else
            {
                pp.checkTopology(false, &points);

                if (pTyp == primitivePatch::OPEN)
                {
                    Info<< pp.name() << '\t' << pp.size()
                        << '\t' << pp.nPoints()
                        << "\tok (not multiply connected)" << endl;
                }
                else
                {
                    Info<< pp.name() << '\t' << pp.size()
                        << '\t' << pp.nPoints()
                        << "\tmultiply connected surface (shared edge)" << endl;
                }
            }
        }

        if (points.size() > 0)
        {
            Info<< "Written " << points.size() << " conflicting points to set "
                << points.name() << endl;

            points.write();
        }
        Info<< endl;
    }
