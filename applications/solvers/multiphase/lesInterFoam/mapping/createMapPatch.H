IOdictionary planeToPatchDict
(
    IOobject
    (
        "planeToPatchDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Switch planeMapping
(
    planeToPatchDict.lookup("mappingPlane")
);

if (planeMapping)
{
    Info << "Plane mapping on." << endl;
}
else
{
    Info << "Plane mapping off." << endl;
}

cuttingPlane cut(mesh, plane(planeToPatchDict));

primitiveFacePatch planePatch
(   
    cut.faces(),
    cut.points()
);

word toPatchName
(
    planeToPatchDict.lookup("targetPatchName")
);


const fvPatchList& patches = mesh.boundary();
bool targetPatchNamePresent(false);
label targetPatchNumber = -1;

forAll(patches, patchI)
{
    if(patches[patchI].name() == toPatchName)
    {
        targetPatchNamePresent = true;
        targetPatchNumber = patchI;
    }
}   

if (!targetPatchNamePresent)
{
    targetPatchNumber = 0;

    if(planeMapping)
    {
        Serr<< "Proc " << Pstream::myProcNo()
            <<": Target patch not present. " << endl;
    }
}
else
{
    if (planeMapping)
    {
        Serr <<"Proc " << Pstream::myProcNo()
             <<": Target patch name: " 
             << patches[targetPatchNumber].name() 
             << endl;
    }
}


PatchToPatchInterpolation<primitiveFacePatch, primitivePatch> interPatch
(
    planePatch,
    patches[targetPatchNumber].patch()
);

//Initialise possible target values
dimensionedVector Ubar("Ubar", U.dimensions(), vector::zero);
dimensionedScalar pBar("pBar", pd.dimensions(), 0.0);

vectorField inletDir(0, vector::zero);

if (targetPatchNamePresent)
{
    inletDir.setSize(U.boundaryField()[targetPatchNumber].size());
    inletDir = U.boundaryField()[targetPatchNumber]
        /max(mag(U.boundaryField()[targetPatchNumber]), VSMALL);
}

mkDir(runTime.path()/runTime.timeName()/"data");

OFstream UbulkFile
(
    runTime.path()/runTime.timeName()/"data"/"Ubulk"
);
UbulkFile.precision(12);

OFstream pTraceFile
(
    runTime.path()/runTime.timeName()/"data"/"ptrace"
);
pTraceFile.precision(12);
