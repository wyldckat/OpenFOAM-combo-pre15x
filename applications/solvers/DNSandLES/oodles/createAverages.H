    label nAveragingSteps = 1;

    IFstream nAveragingStepsFile
    (
        runTime.path()/runTime.timeName()/"nAveragingSteps.raw"
    );


    volVectorField* UmeanPtr;
    volTensorField* RPtr;
    volTensorField* BmeanPtr;
    volScalarField* epsilonMeanPtr;

    volScalarField* pMeanPtr;
    volScalarField* pPrime2MeanPtr;

    if (nAveragingStepsFile.good())
    {
        nAveragingStepsFile >> nAveragingSteps;

        Info<< "Reading field Umean\n" << endl;
        UmeanPtr = new volVectorField
        (
            IOobject
            (
                "Umean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info<< "\nReading field R\n" << endl;
        RPtr = new volTensorField
        (
            IOobject
            (
                "R",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info<< "\nReading field Bmean\n" << endl;
        BmeanPtr = new volTensorField
        (
            IOobject
            (
                "Bmean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info<< "\nReading field epsilonMean\n" << endl;
        epsilonMeanPtr = new volScalarField
        (
            IOobject
            (
                "epsilonMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info<< "Reading field pMean\n" << endl;
        pMeanPtr = new volScalarField
        (
            IOobject
            (
                "pMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        Info<< "Reading field pPrime2Mean\n" << endl;
        pPrime2MeanPtr = new volScalarField
        (
            IOobject
            (
                "pPrime2Mean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
    }
    else
    {
        Info<< "Creating field Umean\n" << endl;
        UmeanPtr = new volVectorField
        (
            IOobject
            (
                "Umean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            U
        );

        Info<< "Creating field R\n" << endl;
        RPtr = new volTensorField
        (
            IOobject
            (
                "R",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqr(U) - sqr(*UmeanPtr)
        );

        Info<< "Creating field Bmean\n" << endl;
        BmeanPtr = new volTensorField
        (
            IOobject
            (
                "Bmean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sgsModel->B()
        );

        Info<< "Creating field epsilonMean\n" << endl;
        epsilonMeanPtr = new volScalarField
        (
            IOobject
            (
                "epsilonMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sgsModel->epsilon()
        );

        Info<< "Creating field pMean\n" << endl;
        pMeanPtr = new volScalarField
        (
            IOobject
            (
                "pMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            p
        );

        Info<< "Creating field pPrime2Mean\n" << endl;
        pPrime2MeanPtr = new volScalarField
        (
            IOobject
            (
                "pPrime2Mean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqr(p) - sqr(*pMeanPtr)
        );
    }

    volVectorField& Umean = *UmeanPtr;
    volTensorField& R = *RPtr;
    volTensorField& Bmean = *BmeanPtr;
    volScalarField& epsilonMean = *epsilonMeanPtr;

    volScalarField& pMean = *pMeanPtr;
    volScalarField& pPrime2Mean = *pPrime2MeanPtr;
