    OFstream file
    (
        runTime.timePath()/"Cs.dat"
    );

    const vectorField& C = aMesh.centres().internalField();

    forAll(C, faceI)
    {
        file<< C[faceI].x() << ' ' 
            << C[faceI].y() << ' '
            << ::acos(C[faceI].x()/mag(C[faceI]))*180.0/M_PI << ' '
            << Cs.internalField()[faceI] << ' '
            << mag(Us.internalField()[faceI]) << endl;
    }
