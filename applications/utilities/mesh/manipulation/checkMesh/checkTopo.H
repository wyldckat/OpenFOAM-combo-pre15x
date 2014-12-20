    // Check if the boundary definition is unique
    mesh.boundaryMesh().checkDefinition(true);
    Pout<< endl;
    mesh.boundaryMesh().checkParallelSync(true);
    Pout<< endl;

    Pout<< "Number of points:           "
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


    if (Pstream::parRun())
    {
        const globalMeshData& parData = mesh.globalData();

        Info<< "Overall number of points:   " << parData.nTotalPoints() << nl
            << "                  faces:    " << parData.nTotalFaces() << nl
            << "                  cells:    " << parData.nTotalCells() << nl
            << endl;
    }


    Pout << nl << "Checking topology and geometry ..." << endl;

#include "checkBasicTopo.H"
    Pout << nl << "Basic topo ok ..." << endl;

#include "checkPatchTopo.H"

    Pout << nl << "Patch topo ok ..." << endl;

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

    Pout << "Topology check done.\n" << endl;
