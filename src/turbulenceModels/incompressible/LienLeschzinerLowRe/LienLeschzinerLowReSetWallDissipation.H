// Set the dissipation in the near-wall cell to the value prescribed by the
// Lien-Leschziner low-Re model with Wolfstein length-scale prescription

    volScalarField boundaryCentral
    (
        IOobject
        (
            "boundaryCentral",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    );

    volScalarField boundarySource
    (
        IOobject
        (
            "boundarySource",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", epsilon_.dimensions()/dimTime, 0.0)
    );

{
    //- use constant Cmu for epsilon in the near-wall cell
    scalar Cmu75 = pow(Cmu.value(), 0.75);

    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (typeid(currPatch) == typeid(wallFvPatch))
        {
            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                boundaryCentral[faceCelli] += 1.0;
                boundarySource[faceCelli] +=
                    1.0e+10
                    *Cmu75*pow(k_[faceCelli], 1.5)
                    /(
                         kappa_.value()*y_[faceCelli]
                         *(1.0 - exp(-Aepsilon.value()*yStar_[faceCelli]))
                     )
                    *exp(-Amu.value()*sqr(yStar_[faceCelli]));

            }
        }
    }

    // perform the averaging

    scalarField& bcInt = boundaryCentral.internalField();

    forAll (bcInt, celli)
    {
        if (boundaryCentral[celli] > 0.0)
        {
            boundarySource[celli] /= boundaryCentral[celli];

            boundaryCentral[celli] = 1.0e+10;
        }
    }
}
