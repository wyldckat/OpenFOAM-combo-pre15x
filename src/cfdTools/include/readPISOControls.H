    const dictionary& piso = mesh.solutionDict().subDict("PISO");

    int nCorr(readInt(piso.lookup("nCorrectors")));

    int nNonOrthCorr = 0;
    if (piso.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(piso.lookup("nNonOrthogonalCorrectors"));
    }

    bool momentumPredictor = true;
    if (piso.found("momentumPredictor"))
    {
        momentumPredictor = Switch(piso.lookup("momentumPredictor"));
    }

    bool fluxGradp = false;
    if (piso.found("fluxGradp"))
    {
        fluxGradp = Switch(piso.lookup("fluxGradp"));
    }

    bool transSonic = false;
    if (piso.found("transSonic"))
    {
        transSonic = Switch(piso.lookup("transSonic"));
    }

    scalar UphiCoeff = 0.9;
    if (piso.found("UphiCoeff"))
    {
        UphiCoeff = readScalar(piso.lookup("UphiCoeff"));
    }

    int nOuterCorr = 1;
    if (piso.found("nOuterCorrectors"))
    {
        nOuterCorr = readInt(piso.lookup("nOuterCorrectors"));
    }

    label pRefCell = 0;
    if (piso.found("pRefCell"))
    {
        pRefCell = readLabel(piso.lookup("pRefCell"));
    }

    label pRefCellNew = findRefCell(mesh, pRefCell);

    if (pRefCellNew != pRefCell)
    {
        pRefCell = pRefCellNew;

        ((dictionary&)piso).remove("pRefCell");
        ((dictionary&)piso).add("pRefCell", pRefCell);
    }

    scalar pRefValue = 0.0;
    if (piso.found("pRefValue"))
    {
        pRefValue = readScalar(piso.lookup("pRefValue"));
    }
