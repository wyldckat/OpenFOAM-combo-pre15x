        mesh.updateTopology();

        // Calculate the new point positions depending on whether the
        // topological change has happened or not
        pointField newPoints;

        vector curVel = motionU*Foam::sin(runTime.value()*M_PI/10.0);
// #       include "hrvsDebug.H"
        if (mesh.morphing())
        {
            Info << "Topology change. Calculating motion points" << endl;

            motionMarkup =
                vertexMarkup
                (
                    mesh.morphMap().preMotionPoints(),
                    curLeft,
                    curRight
                );

            // Create new points by moving old points but using the
            // pre-motion points at the motion selector for the moving
            // region
            newPoints =
                mesh.points()
              + (
                    pos(0.5 - mag(motionMarkup)) // cells above the body
                    /*
                  + pos(motionMarkup - 0.5)*      // cells in front of the body
                    (
                        mesh.points().component(vector::X)/curRight
                    )
                    */
//                + pos(-motionMarkup - 0.5)*      // cells behind the body
//                  (
//                      (
//                          mesh.points().component(vector::X)
//                        - leftEdge
//                      )/
//                      (curLeft - leftEdge)
//                  )
                )*curVel*runTime.deltaT().value();
        }
        else
        {
            Info << "No topology change" << endl;
            // Set the mesh motion
            newPoints =
                mesh.points()
              + (
                    pos(0.5 - mag(motionMarkup)) // cells above the body
                    /*
                  + pos(motionMarkup - 0.5)*  // cells in front of the body
                    (
                        mesh.points().component(vector::X)/curRight
                    )
                    */
                    /*
                  + pos(-motionMarkup - 0.5)*  // cells behind the body
                    (
                        (
                            mesh.points().component(vector::X)
                          - leftEdge
                        )/
                        (curLeft - leftEdge)
                    )
                    */
                )*curVel*runTime.deltaT().value();
        }

        curLeft += curVel.x()*runTime.deltaT().value();
        curRight += curVel.x()*runTime.deltaT().value();
        Info << "curLeft: " << curLeft << " curRight: " << curRight << endl;

        // The mesh now contains the cells with zero volume

        Info << "Executing mesh motion" << endl;

        mesh.movePoints(newPoints);

        Info<< "minV " << min(mesh.V()) << endl;
        Info<< "minV0 " << min(mesh.V0()) << endl;

#       include "volContinuity.H"

        Info << "Reset wall velocity for moving wall." << endl;
        refCast<movingWallVelocityFvPatchVectorField>(U.boundaryField()[0]) ==
            curVel;
