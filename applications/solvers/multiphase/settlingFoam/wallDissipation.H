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
        dimensionedScalar("zero", rho.dimensions()/dimTime, 0.0)
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
        dimensionedScalar
        (
            "zero",
            rho.dimensions()*epsilon.dimensions()/dimTime,
            0.0
        )
    );


{
    scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    scalar Cmu75 = ::pow(Cmu.value(), 0.75);

    const fvPatchList& patches = mesh.boundary();

    //- Initialise the near-wall P field to zero
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                G[faceCelli] = 0.0;
            }
        }
    }

    //- Accumulate the wall face contributions to P and boundarySource
    //  Increment boundaryCentral for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            const scalarField& rhow = rho.boundaryField()[patchi];

            const scalarField muw = mul.boundaryField()[patchi];
            const scalarField& mutw = mut.boundaryField()[patchi];

            scalarField magFaceGradU =
                mag(U.boundaryField()[patchi].snGrad());

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                scalar yPlus =
                    Cmu25*y[patchi][facei]*::sqrt(k[faceCelli])
                   /(muw[facei]/rhow[facei]);

                // For corner cells (with two boundary or more faces),
                // boundarySource and G in the near - wall cell
                // are calculated as an average

                boundaryCentral[faceCelli] += 1.0;

                boundarySource[faceCelli] +=
                    1.0e+10
                   *Cmu75*rho[faceCelli]*::pow(k[faceCelli], 1.5)
                   /(kappa.value()*y[patchi][facei]);

                if (yPlus > 11.6)
                {
                    G[faceCelli] +=
                        mutw[facei]*magFaceGradU[facei]
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
        if (boundaryCentral[celli] > 0.5)
        {
            Gint[celli] /= boundaryCentral[celli];
            boundarySource[celli] /= boundaryCentral[celli];

            boundaryCentral[celli] = 1.0e+10*rho[celli];
        }
    }
}
