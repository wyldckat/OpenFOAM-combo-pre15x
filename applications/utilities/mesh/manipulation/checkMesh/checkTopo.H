    // Check if the boundary definition is unique
    mesh.boundaryMesh().checkDefinition(true);
    Info<< endl;

    Info<< "Number of points:           "
        << mesh.points().size() << nl
        << "          edges:            "
        << mesh.edges().size() << nl
        << "          faces:            "
        << mesh.faces().size() << nl
        << "          internal faces:   "
        << mesh.faceNeighbour().size() << nl
        << "          cells:            "
        << mesh.cells().size() << nl
        << "          boundary patches: "
        << mesh.boundaryMesh().size() << nl

        << "          point zones:      "
        << mesh.pointZones().size() << nl

        << "          face zones:       "
        << mesh.faceZones().size() << nl

        << "          cell zones:       "
        << mesh.cellZones().size() << nl;

    Info << "\nChecking topology and geometry ..." << endl;

#include "checkBasicTopo.H"
#include "checkPatchTopo.H"

    mesh.cells().size();
    mesh.faces().size();
    mesh.edges().size();
    mesh.points().size();
    mesh.faceOwner().size();
    mesh.faceNeighbour().size();
    mesh.cellCells().size();
    mesh.edgeCells().size();
    mesh.pointCells().size();
    mesh.edgeFaces().size();
    mesh.pointFaces().size();
    mesh.cellEdges().size();
    mesh.faceEdges().size();
    mesh.pointEdges().size();

    Info << "Topology check done.\n" << endl;
