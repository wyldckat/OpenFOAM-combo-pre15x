
    {
        pointSet points(mesh, "unusedPoints", mesh.nPoints()/100);
        if (mesh.checkPoints(true, &points))
        {
            noFailedChecks++;

            Info<< "Writing " << points.size() << " unused points to set "
                << points.name() << endl << endl;
            points.write();
        }
    }

    {
        faceSet faces(mesh, "upperTriangularFace", mesh.nFaces()/100);
        if (mesh.checkUpperTriangular(true, &faces))
        {
            noFailedChecks++;

            Info<< "Writing " << faces.size() << " unordered faces to set "
                << faces.name() << endl << endl;
            faces.write();
        }
    }

    {
        cellSet cells(mesh, "zipUpCells", mesh.nCells()/100);
        if (mesh.checkCellsZipUp(true, &cells))
        {
            noFailedChecks++;

            Info<< "Writing " << cells.size() << " cells with over used edges"
                << " to set " << cells.name() << endl << endl;
            cells.write();
        }
    }

    {
        faceSet faces(mesh, "outOfRangeFaces", mesh.nFaces()/100);
        if (mesh.checkFaceVertices(true, &faces))
        {
            noFailedChecks++;

            Info<< "Writing " << faces.size() << " faces with out-of-range"
                << " vertices to set " << faces.name() << endl << endl;
            faces.write();
        }
    }

    {
        faceSet faces(mesh, "edgeFaces", mesh.nFaces()/100);
        if (mesh.checkFaceFaces(true, &faces))
        {
            noFailedChecks++;

            Info<< "Writing " << faces.size() << " faces with incorrect"
                << " edges to set " << faces.name() << endl << endl;
            faces.write();
        }
    }

    {
        cellSet cells(mesh, "floatingCells", mesh.nCells()/100);
        if (mesh.checkFloatingCells(true, &cells))
        {
            noFailedChecks++;

            Info<< "Writing " << cells.size() << " cells connected to too"
                << " few other cells to set "
                << cells.name() << endl << endl;
            cells.write();
        }
    }
