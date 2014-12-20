    Info<< "Reading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField beta
    (
        IOobject
        (
            "beta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1.0 - alpha
        //,alpha.boundaryField().types()
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Ua\n" << endl;
    volVectorField Ua
    (
        IOobject
        (
            "Ua",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Ub\n" << endl;
    volVectorField Ub
    (
        IOobject
        (
            "Ub",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha*Ua + beta*Ub
    );

    Info<< "Reading field k\n" << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field epsilon\n" << endl;
    volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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

    dimensionedScalar rhoa
    (
        transportProperties.lookup("rhoa")
    );

    dimensionedScalar rhob
    (
        transportProperties.lookup("rhob")
    );

    dimensionedScalar nua
    (
        transportProperties.lookup("nua")
    );

    dimensionedScalar nub
    (
        transportProperties.lookup("nub")
    );

    dimensionedScalar da
    (
        transportProperties.lookup("da")
    );

    dimensionedScalar db
    (
        transportProperties.lookup("db")
    );

    dimensionedScalar Cvm
    (
        transportProperties.lookup("Cvm")
    );

    dimensionedScalar Cl
    (
        transportProperties.lookup("Cl")
    );

    dimensionedScalar Ct
    (
        transportProperties.lookup("Ct")
    );

#   include "createPhia.H"
#   include "createPhib.H"

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh
        ),
        fvc::interpolate(alpha)*phia
      + fvc::interpolate(beta)*phib
    );


    IOdictionary turbulenceProperties
    (
        IOobject
        (
            "turbulenceProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Switch turbulence
    (
        turbulenceProperties.lookup("turbulence")
    );

    dictionary kEpsilonCoeffs
    (
        turbulenceProperties.subDict("kEpsilonCoeffs")
    );


    dimensionedScalar Cmu
    (
        kEpsilonCoeffs.lookup("Cmu")
    );

    dimensionedScalar C1
    (
        kEpsilonCoeffs.lookup("C1")
    );

    dimensionedScalar C2
    (
        kEpsilonCoeffs.lookup("C2")
    );

    dimensionedScalar alphak
    (
        kEpsilonCoeffs.lookup("alphak")
    );

    dimensionedScalar alphaEps
    (
        kEpsilonCoeffs.lookup("alphaEps")
    );

    dictionary wallFunctionCoeffs
    (
        turbulenceProperties.subDict("wallFunctionCoeffs")
    );

    dimensionedScalar kappa
    (
        wallFunctionCoeffs.lookup("kappa")
    );

    dimensionedScalar E
    (
        wallFunctionCoeffs.lookup("E")
    );

    nearWallDist y(mesh);

    Info<< "Calculating field nutb\n" << endl;
    volScalarField nutb
    (
        IOobject
        (
            "nutb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Cmu*sqr(k)/epsilon
    );

    Info<< "Calculating field nuEffa\n" << endl;
    volScalarField nuEffa
    (
        IOobject
        (
            "nuEffa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(Ct)*nutb + nua
    );

    Info<< "Calculating field nuEffb\n" << endl;
    volScalarField nuEffb
    (
        IOobject
        (
            "nuEffb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nutb + nub
    );


    Info<< "Calculating field DDtUa and DDtUb\n" << endl;

    volVectorField DDtUa =
        fvc::ddt(Ua)
      + fvc::div(phia, Ua)
      - fvc::div(phia)*Ua;

    volVectorField DDtUb =
        fvc::ddt(Ub)
      + fvc::div(phib, Ub)
      - fvc::div(phib)*Ub;


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
