        mesh.updateTopology();

        // Calculate the new point positions depending on whether the
        // topological change has happened or not
        pointField newPoints;

        if (mesh.morphing())
        {
            Info << "Topology change. Calculating motion points" << endl;
            const pointField& preMotionPoints =
                mesh.morphMap().preMotionPoints();

            // Create new points by moving old points but using the
            // pre-motion points at the motion selector for the moving
            // region
            newPoints =
                mesh.basicMesh().points()
              + pos(preMotionPoints.component(vector::X) - motionRefX)*
                motionU*Foam::cos(runTime.value()*M_PI/20.0)*
                runTime.deltaT().value();

            Info << "number of cells: " << mesh.basicMesh().nCells() << endl;
            Info << "number of faces: " << mesh.basicMesh().nFaces() << endl;
//             Info << "owner: " << mesh.owner() << endl;
//             Info << "neighbour: " << mesh.neighbour() << endl;
        }
        else
        {
            Info << "No topology change" << endl;
            // Set the mesh motion

            newPoints =
                mesh.basicMesh().points()
              + pos
                (
                    mesh.basicMesh().points().component(vector::X)
                  - motionRefX
                )*motionU*Foam::cos(runTime.value()*M_PI/20.0)*
                runTime.deltaT().value();
        }

        // The mesh now contains the cells with zero volume

        Info << "Executing mesh motion" << endl;
        mesh.movePoints(newPoints);

        mesh.V();
        Info<< "minV " << min(mesh.V()) << endl;
        Info<< "minV0 " << min(mesh.V0()) << endl;

//         Info << "newPoints: " << newPoints << endl;
//         Info << "mesh.Sf(): " << mesh.Sf() << endl;
//         Info << "mesh.V(): " << mesh.V() << endl;
        //  The mesh now has got non-zero volume cells

//         Info << "mesh.V() after motion: " << mesh.V() << endl;
