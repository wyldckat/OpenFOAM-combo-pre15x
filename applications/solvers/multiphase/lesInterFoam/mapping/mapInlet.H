if (planeMapping && targetPatchNamePresent && cut.cut())
{
    //SGS model mapping
    turbulence->k()().boundaryField()[targetPatchNumber] ==
        interPatch.faceInterpolate(cut.sample(turbulence->k()()));

    Ubar =
    (
        planeToPatchDict.lookup("Ubar")
    );

    U.boundaryField()[targetPatchNumber] ==
        interPatch.faceInterpolate(cut.sample(U));


    scalar Q = sum
    (
        mesh.Sf().boundaryField()[targetPatchNumber]
        & U.boundaryField()[targetPatchNumber]
    );

    scalar Qbar = sum
    (
        mesh.Sf().boundaryField()[targetPatchNumber] &
        Ubar.value()
    );

    pBar.value() =
    (
        sum(pd.boundaryField()[targetPatchNumber]
        *mesh.boundary()[targetPatchNumber].magSf())
        /sum(mesh.magSf().boundaryField()[targetPatchNumber])
    );

    Info << "Correction factor: " << Qbar/Q
         << ", Uncorrected velocity: " 
         << mag(Q/sum(mesh.magSf().boundaryField()[targetPatchNumber])) 
         << ", Target velocity: "
         << mag(Ubar).value()
         << ", Inlet pressure: " << pBar.value()
         << endl;

    if(Q != 0)
    {
        U.boundaryField()[targetPatchNumber] ==
            U.boundaryField()[targetPatchNumber] * Qbar/Q;
    }


    phi.boundaryField()[targetPatchNumber] ==
    (
        mesh.Sf().boundaryField()[targetPatchNumber] &
        U.boundaryField()[targetPatchNumber]
    );
        
    
    UbulkFile << runTime.value() << tab << mag(Ubar.value()) 
              << endl;
    pTraceFile << runTime.value() << tab 
               << pBar.value() << endl;
}
