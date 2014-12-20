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

Description
    Simple central-difference snGrad scheme with non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#include "correctedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
correctedSnGrad<Type>::~correctedSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
correctedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // construct GeometricField<Type, fvPatchField, surfaceMesh>
    tmp<GeometricField<Type, fvPatchField, surfaceMesh> > tssf
    (
        new GeometricField<Type, fvPatchField, surfaceMesh>
        (
            IOobject
            (
                "snGradCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, fvPatchField, surfaceMesh>& ssf = tssf();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        ssf.replace
        (
            cmpt,
            mesh.correctionVectors()
          & linear
            <
                typename 
                outerProduct<vector, typename pTraits<Type>::cmptType>::type
            >(mesh).interpolate
            (
                gradScheme<typename pTraits<Type>::cmptType>::New
                (
                    mesh,
                    mesh.gradScheme(ssf.name())
                )()
                //gaussGrad<typename pTraits<Type>::cmptType>(mesh)
               .grad(vf.component(cmpt))
            )
        );
    }

#   ifdef NEW_NON_ORTH
    // Set reference to difference factors array
    const scalarField& DeltaCoeffs = mesh.deltaCoeffs().internalField();

    // Owner/neighbour addressing
    const volVectorField& cellCentres = mesh.C();
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    const surfaceVectorField& areas = mesh.Sf();
    const surfaceScalarField& magAreas = mesh.magSf();
    const surfaceVectorField& corrVecs = mesh.correctionVectors();

    forAll(owner, faceI)
    {
        vector delta =
            cellCentres[neighbour[faceI]] - cellCentres[owner[faceI]];

        vector unitArea = areas[faceI]/magAreas[faceI];

        scalar orthCorrCoeff =
            (delta & (unitArea - corrVecs[faceI]))/mag(delta);

        Type orthCorr =
            (orthCorrCoeff/mag(delta) - DeltaCoeffs[faceI])
           *(vf[neighbour[faceI]] - vf[owner[faceI]]);

        ssf[faceI] += orthCorr;
    }
#   endif

    return tssf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
