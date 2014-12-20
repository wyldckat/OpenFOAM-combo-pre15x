    const dictionary& simple = mesh.solutionDict().subDict("SIMPLE");

    int nNonOrthCorr = 0;
    if (simple.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(simple.lookup("nNonOrthogonalCorrectors"));
    }

    bool momentumPredictor = true;
    if (simple.found("momentumPredictor"))
    {
        momentumPredictor = Switch(simple.lookup("momentumPredictor"));
    }

    bool fluxGradp = false;
    if (simple.found("fluxGradp"))
    {
        fluxGradp = Switch(simple.lookup("fluxGradp"));
    }

    scalar UphiCoeff = 0.9;
    if (simple.found("UphiCoeff"))
    {
        UphiCoeff = readScalar(simple.lookup("UphiCoeff"));
    }

    label pRefCell = 0;
    if (simple.found("pRefCell"))
    {
        pRefCell = readLabel(simple.lookup("pRefCell"));
    }

    label pRefCellNew = findRefCell(mesh, pRefCell);

    if (pRefCellNew != pRefCell)
    {
        pRefCell = pRefCellNew;

        ((dictionary&)simple).remove("pRefCell");
        ((dictionary&)simple).add("pRefCell", pRefCell);
    }

    scalar pRefValue = 0.0;
    if (simple.found("pRefValue"))
    {
        pRefValue = readScalar(simple.lookup("pRefValue"));
    }
