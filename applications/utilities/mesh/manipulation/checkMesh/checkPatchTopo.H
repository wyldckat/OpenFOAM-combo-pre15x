    if (!Pstream::parRun())
    {
        Pout<< "\nChecking patch topology for multiply connected surfaces ..."
            << endl;

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // Non-manifold points
        pointSet points
        (
            mesh,
            "nonManifoldPoints",
            mesh.nPoints()/100
        );

        Pout << "\nPatch\tFaces\tPoints\tSurface" << endl;
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];


            primitivePatch::surfaceTopo pTyp = pp.surfaceType();

            if (pp.size() == 0)
            {
                Pout<< pp.name() << '\t' << pp.size()
                    << '\t' << pp.nPoints()
                    << "\tok (empty)" << endl;
            }
            else if (pTyp == primitivePatch::MANIFOLD)
            {
                if (pp.checkPointManifold(true, &points))
                {
                    Pout<< pp.name() << '\t' << pp.size()
                        << '\t' << pp.nPoints()
                        << "\tmultiply connected (shared point)" << endl;
                }
                else
                {
                    Pout<< pp.name() << '\t' << pp.size()
                        << '\t' << pp.nPoints()
                        << "\tok (closed singly connected surface)" << endl;
                }

                // Add points on non-manifold edges to make set complete
                pp.checkTopology(false, &points);
            }
            else
            {
                pp.checkTopology(false, &points);

                if (pTyp == primitivePatch::OPEN)
                {
                    Pout<< pp.name() << '\t' << pp.size()
                        << '\t' << pp.nPoints()
                        << "\tok (not multiply connected)" << endl;
                }
                else
                {
                    Pout<< pp.name() << '\t' << pp.size()
                        << '\t' << pp.nPoints()
                        << "\tmultiply connected surface (shared edge)" << endl;
                }
            }
        }

        if (points.size() > 0)
        {
            Pout<< "Written " << points.size() << " conflicting points to set "
                << points.name() << endl;

            points.write();
        }
        Pout<< endl;
    }
