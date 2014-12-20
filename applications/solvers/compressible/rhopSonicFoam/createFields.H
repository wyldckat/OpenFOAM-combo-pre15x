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
    p.oldTime();

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        1.0/(R*T)
    );
    psi.oldTime();

    Info<< "Reading field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh 
    );

    rho == p*psi;

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

    volVectorField rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    rhoU == rho*U;

    volScalarField rhoE
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    rhoE == rho*Cv*T + 0.5*rho*magSqr(rhoU/rho);

#   include "createPhi.H"
    phi.oldTime();

    surfaceScalarField phiv
    (
        IOobject
        (
            "phiv",
            runTime.timeName(),
            mesh
        ),
        phi/fvc::interpolate(rho),
        phi.boundaryField().types()
    );
