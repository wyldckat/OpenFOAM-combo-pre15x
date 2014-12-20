{

#   include "UEqn.H"

    wordList pcorrTypes(pd.boundaryField().types());

    for (label i=0; i<pd.boundaryField().size(); i++)
    {
        if(pd.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", pd.dimensions(), 0.0),
        pcorrTypes
    );

    surfaceScalarField rUAf("(rho|A(U))", fvc::interpolate(1.0/UEqn.A()));

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rUAf, pcorr) == fvc::div(phi)
        );

        fvScalarMatrix::reference pRef =
            pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();
        pcorrEqn.unsetReference(pRef);

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pcorrEqn.flux();
        }
    }
}

#include "continuityErrs.H"

