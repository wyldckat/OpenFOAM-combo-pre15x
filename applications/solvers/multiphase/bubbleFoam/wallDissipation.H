// Initialise the boundary condition fields for epsilon to zero
// They are used to accumulate the information for all the wall faces
// of a cell which is subsequently averaged.

    volScalarField boundaryCentral
    (
        IOobject
        (
            "boundaryCentral",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    );

    volScalarField boundarySource
    (
        IOobject
        (
            "boundarySource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", epsilon.dimensions()/dimTime, 0.0)
    );


{
    scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    scalar Cmu75 = ::pow(Cmu.value(), 0.75);

    const fvPatchList& patches = mesh.boundary();

    //- Initialise the near-wall P field to zero
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (typeid(currPatch) == typeid(wallFvPatch))
        {
            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                G[faceCelli] = 0.0;
            }
        }
    }

    //- Accumulate the wall face contributions to P and boundarySource
    //  Increment boundaryCentral for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (typeid(currPatch) == typeid(wallFvPatch))
        {
            const scalarField& nuw = nutb.boundaryField()[patchi];

            scalarField magFaceGradU =
                mag(((const volVectorField&)U).boundaryField()[patchi]
               .snGrad());

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                scalar yPlus =
                    Cmu25*y[patchi][facei]
                    *::sqrt(k[faceCelli])
                    /nub.value();

                // for corner cells (with two boundary or more faces),
                // boundarySource and P in the near - wall cell
                // are calculated as an average

                boundaryCentral[faceCelli] += 1.0;
                boundarySource[faceCelli] +=
                    1.0e+10
                    *Cmu75*::pow(k[faceCelli], 1.5)
                    /(kappa.value()*y[patchi][facei]);

                if (yPlus > 11.6)
                {
                    G[faceCelli] +=
                        nuw[facei]*magFaceGradU[facei]
                        *Cmu25*::sqrt(k[faceCelli])
                        /(kappa.value()*y[patchi][facei]);
                }
            }
        }
    }


    // perform the averaging

    scalarField& Gint = G.internalField();

    forAll (Gint, celli)
    {
        if (boundaryCentral[celli] > 0.0)
        {
            Gint[celli] /= boundaryCentral[celli];
            boundarySource[celli] /= boundaryCentral[celli];

            boundaryCentral[celli] = 1.0e+10;
        }
    }
}
