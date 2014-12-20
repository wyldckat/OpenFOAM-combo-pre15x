    if (runTime.outputTime())
    {
        // Displacement gradient
        tetPointTensorField gradU = tetFec::grad(U);

        // Stress tensor
        tetPointTensorField sigma = mu.value()*(gradU + gradU.T())
          + lambda.value()*tr(gradU)*I;

        // sigmaXX
        tetPointScalarField sigmaXX
        (
            IOobject
            (
                "sigmaXX",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(tensor::XX)
        );

        // sigmaYY
        tetPointScalarField sigmaYY
        (
            IOobject
            (
                "sigmaYY",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(tensor::YY)
        );

        // sigmaXY
        tetPointScalarField sigmaXY
        (
            IOobject
            (
                "sigmaXY",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(tensor::XY)
        );

        runTime.write();
    }
