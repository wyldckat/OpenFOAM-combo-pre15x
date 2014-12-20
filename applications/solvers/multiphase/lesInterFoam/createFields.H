    label nGammaCorr
    (
        readLabel(piso.lookup("nGammaCorr"))
    );

    label nGammaSubCycles
    (
        readLabel(piso.lookup("nGammaSubCycles"))
    );

    Info<< "Reading field pd\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field gamma\n" << endl;
    volScalarField gamma
    (
        IOobject
        (
            "gamma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#   include "createPhi.H"

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );



    autoPtr<transportModel> twoPhaseLaminar
    (
        transportModel::New(U, phi,"twoPhase")
    );
    
    dimensionedScalar rho1
    (
        dictionary(transportProperties.subDict("phase1")).lookup("rho")
    );

    dimensionedScalar rho2
    (
        dictionary(transportProperties.subDict("phase2")).lookup("rho")
    );


    // Need to store rho for ddt(rho, U)
    volScalarField rho(gamma*rho1 + (1.0- gamma)*rho2);

    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // gamma solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );


    Info<< "Calculating field g.h\n" << endl;
    surfaceScalarField ghf = fvc::interpolate(g & mesh.C());


    // Construct interface from gamma distribution
    interfaceProperties interface(gamma, U, transportProperties);

    // Construct LES model
    autoPtr<LESmodel> turbulence
    (
        LESmodel::New(U, phi, twoPhaseLaminar())
    );
