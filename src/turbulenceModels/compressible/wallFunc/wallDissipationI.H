/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Global
    wallDissipation

Description
    Calculate wall dissipation from wall-functions.

\*---------------------------------------------------------------------------*/

// Initialise the boundary condition fields for epsilon to zero
// They are used to accumulate the information for all the wall faces
// of a cell which is subsequently averaged.

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
        dimensionedScalar("0", rho_.dimensions()/dimTime, 0.0)
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
        dimensionedScalar
        ("0", rho_.dimensions()*epsilon_.dimensions()/dimTime, 0.0)
    );


{
    scalar Cmu25 = pow(Cmu.value(), 0.25);
    scalar Cmu75 = pow(Cmu.value(), 0.75);

    const fvPatchList& patches = mesh_.boundary();

    //- Initialise the near-wall G field to zero
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (typeid(curPatch) == typeid(wallFvPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                G[faceCelli] = 0.0;
            }
        }
    }

    //- Accumulate the wall face contributions to G and boundarySource
    //  Increment boundaryCentral for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (typeid(curPatch) == typeid(wallFvPatch))
        {
            const scalarField& rhow = rho_.boundaryField()[patchi];

            const scalarField muw = mu().boundaryField()[patchi];
            const scalarField& mutw = mut_.boundaryField()[patchi];

            scalarField magFaceGradU =
                mag(U_.boundaryField()[patchi].snGrad());

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                scalar yPlus =
                    Cmu25*turbulenceModel::y_[patchi][facei]*sqrt(k_[faceCelli])
                    /(muw[facei]/rhow[facei]);

                // For corner cells (with two boundary or more faces),
                // boundarySource and G in the near - wall cell
                // are calculated as an average

                boundaryCentral[faceCelli] += 1.0;

                boundarySource[faceCelli] +=
                    1.0e+10
                    *Cmu75*rho_[faceCelli]*pow(k_[faceCelli], 1.5)
                    /(kappa_.value()*turbulenceModel::y_[patchi][facei]);

                if (yPlus > yPlusLam_)
                {
                    G[faceCelli] +=
                        mutw[facei]*magFaceGradU[facei]
                      //- *0.5*log(2.0*yPlus/yPlusLam_)    //  UMIST correction
                        *Cmu25*sqrt(k_[faceCelli])
                        /(kappa_.value()*turbulenceModel::y_[patchi][facei]);
                }
            }
        }
    }


    // Perform the averaging

    scalarField& Gint = G.internalField();

    forAll (Gint, celli)
    {
        if (boundaryCentral[celli] > 0.5)
        {
            Gint[celli] /= boundaryCentral[celli];
            boundarySource[celli] /= boundaryCentral[celli];

            boundaryCentral[celli] = 1.0e+10*rho_[celli];
        }
    }
}


// ************************************************************************* //
