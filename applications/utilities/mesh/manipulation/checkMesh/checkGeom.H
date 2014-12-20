    boundBox bb(mesh.points());

    Info<< "Domain bounding box: min = "
        << bb.min() << "    max = " << bb.max() << " meters. " << endl;

    Info << "\nChecking geometry..." << endl;

#   include "checkBasicGeom.H"

    mesh.cellCentres().size();
    mesh.faceCentres().size();
    mesh.edgeVectors().size();
    mesh.cellVolumes().size();
    mesh.faceAreas().size();

    Info << "Geometry check done.\n" << endl;
