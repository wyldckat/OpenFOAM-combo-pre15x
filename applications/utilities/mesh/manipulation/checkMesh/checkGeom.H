    boundBox bb(mesh.points());

    Pout<< "Domain bounding box: min = "
        << bb.min() << "    max = " << bb.max() << " meters. " << endl;

    if (Pstream::parRun())
    {
        const boundBox& bb = mesh.parallelData().bb();

        Info<< "Overall domain bounding box: min = "
            << bb.min() << "    max = " << bb.max() << " meters. " << endl;
    }

    Pout << "\nChecking geometry..." << endl;

#   include "checkBasicGeom.H"

    mesh.cellCentres().size();
    mesh.faceCentres().size();
    mesh.cellVolumes().size();
    mesh.faceAreas().size();

    Pout << "Geometry check done.\n" << endl;
