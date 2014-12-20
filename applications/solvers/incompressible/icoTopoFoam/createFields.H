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

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
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

    volScalarField A
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("A", dimless/dimTime, 1),
        zeroGradientFvPatchScalarField::typeName
    );

    volVectorField H
    (
        IOobject
        (
            "H",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("H", A.dimensions()*U.dimensions(), vector::zero),
        zeroGradientFvPatchVectorField::typeName
    );

    surfaceScalarField phiField
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U) & mesh.Sf()
    );

    scalar totalVolume = sum(mesh.V());
