    
    if (mesh.checkClosedBoundary(true)) noFailedChecks++;

    {
        cellSet cells(mesh, "nonClosedCells", mesh.nCells()/100 + 1);
        if (mesh.checkClosedCells(true, &cells))
        {
            noFailedChecks++;

            if (cells.size() > 0)
            {
                Info<< "Writing " << cells.size() << " non closed cells to set "
                    << cells.name() << endl << endl;
                cells.write();
            }
        }
    }

    {
        faceSet faces(mesh, "zeroAreaFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceAreas(true, &faces))
        {
            noFailedChecks++;

            if (faces.size() > 0)
            {
                Info<< "Writing " << faces.size() << " zero area faces to set "
                    << faces.name() << endl << endl;
                faces.write();
            }
        }
    }

    {
        cellSet cells(mesh, "zeroVolumeCells", mesh.nCells()/100 + 1);
        if (mesh.checkCellVolumes(true, &cells))
        {
            noFailedChecks++;

            if (cells.size() > 0)
            {
                Info<< "Writing " << cells.size() << " zero volume cells to"
                    << " set " << cells.name() << endl << endl;
                cells.write();
            }
        }
    }

    {
        faceSet faces(mesh, "nonOrthoFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceDotProduct(true, &faces))
        {
            noFailedChecks++;
        }

        if (faces.size() > 0)
        {
            Info<< "Writing " << faces.size() << " non-orthogonal faces to"
                << " set " << faces.name() << endl << endl;
            faces.write();
        }
    }


    {
        faceSet faces(mesh, "wrongOrientedFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFacePyramids(true, -SMALL, &faces))
        {
            noFailedChecks++;

            if (faces.size() > 0)
            {
                Info<< "Writing " << faces.size() << " faces with incorrect"
                    << " orientation to set " << faces.name() << endl << endl;
                faces.write();
            }
        }
    }

    {
        faceSet faces(mesh, "skewFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceSkewness(true, &faces))
        {
            noFailedChecks++;

            if (faces.size() > 0)
            {
                Info<< "Writing " << faces.size() << " skew faces to set "
                    << faces.name() << endl << endl;
                faces.write();
            }
        }
    }

    if (!args.options().found("noCoords"))
    {
        pointSet points(mesh, "nearPoints", mesh.nPoints()/100 + 1);
        if (checkCoords(mesh, true, SMALL, &points))
        {
            //noFailedChecks++;

            if (points.size() > 0)
            {
                Info<< "Writing " << points.size() << " near points to set "
                    << points.name() << endl << endl;
                points.write();
            }
        }
    }

    {
        pointSet points(mesh, "shortEdges", mesh.nEdges()/1000 + 1);
        if (checkEdges(mesh, true, SMALL, &points))
        {
            //noFailedChecks++;

            if (points.size() > 0)
            {
                Info<< "Writing " << points.size()
                    << " points on short edges to set "
                    << points.name() << endl << endl;
                points.write();
            }
        }
    }

    {
        faceSet faces(mesh, "concaveFaces", mesh.nFaces()/100 + 1);
        if
        (
            mesh.checkFaceAngles
            (
                true,
                10,
                &faces
            )
        )
        {
            //noFailedChecks++;

            if (faces.size() > 0)
            {
                Info<< "Writing " << faces.size()
                    << " faces with concave angles to set "
                    << faces.name() << endl << endl;
                faces.write();
            }
        }
    }
