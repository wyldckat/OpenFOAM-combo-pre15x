Info<< "Reading alignPatchDict\n" << endl;

IOdictionary alignPatchDict
(
    IOobject
    (
        "alignPatchDict",
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

word patchName(alignPatchDict.lookup("patch"));
vector planeNormal(alignPatchDict.lookup("planeNormal"));
