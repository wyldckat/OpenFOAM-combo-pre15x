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
    
\*---------------------------------------------------------------------------*/

#include "gaussLaplacianScheme.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<fvMatrix<Type> >
gaussLaplacianScheme<Type>::fvmLaplacian
(
    const surfaceScalarField& gamma,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<surfaceScalarField> tdeltaCoeffs =
        this->tsnGradScheme_().deltaCoeffs(vf);
    const surfaceScalarField& deltaCoeffs = tdeltaCoeffs();

    surfaceScalarField gammaMagSf = gamma*this->mesh().magSf();

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    fvm.upper() = deltaCoeffs.internalField()*gammaMagSf.internalField();
    fvm.negSumDiag();

    forAll(fvm.psi().boundaryField(), patchI)
    {
        const fvPatchField<Type>& psf = fvm.psi().boundaryField()[patchI];
        const fvPatchScalarField& patchGamma =
            gammaMagSf.boundaryField()[patchI];

        fvm.internalCoeffs()[patchI] = patchGamma*psf.gradientInternalCoeffs();
        fvm.boundaryCoeffs()[patchI] = -patchGamma*psf.gradientBoundaryCoeffs();
    }

    if (this->tsnGradScheme_().corrected())
    {
        if (this->mesh().fluxRequired(vf.name()))
        {
            fvm.faceFluxCorrectionPtr() = new
            GeometricField<Type, fvPatchField, surfaceMesh>
            (
                gammaMagSf*this->tsnGradScheme_().correction(vf)
            );

            fvm.source() -=
                this->mesh().V()*
                fvc::div
                (
                    *fvm.faceFluxCorrectionPtr()
                )().internalField();
        }
        else
        {
            fvm.source() -=
                this->mesh().V()*
                fvc::div
                (
                    gammaMagSf*this->tsnGradScheme_().correction(vf)
                )().internalField();
        }
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
gaussLaplacianScheme<Type>::fvmLaplacian
(
    const surfaceTensorField& gamma,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvmLaplacian(tr(gamma), vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
gaussLaplacianScheme<Type>::fvcLaplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tLaplacian
    (
        fvc::div(this->tsnGradScheme_().snGrad(vf)*vf.mesh().magSf())
    );

    tLaplacian().rename("laplacian(" + vf.name() + ')');

    return tLaplacian;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
gaussLaplacianScheme<Type>::fvcLaplacian
(
    const surfaceScalarField& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tLaplacian
    (
        fvc::div(gamma*this->tsnGradScheme_().snGrad(vf)*vf.mesh().magSf())
    );

    tLaplacian().rename("laplacian(" + gamma.name() + ',' + vf.name() + ')');

    return tLaplacian;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
gaussLaplacianScheme<Type>::fvcLaplacian
(
    const surfaceTensorField& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    surfaceVectorField SfGamma = mesh.Sf() & gamma;

    surfaceVectorField Sn = mesh.Sf()/mesh.magSf();
    surfaceVectorField SfGammaCorr = SfGamma - (SfGamma & Sn)*Sn;


    tmp<GeometricField<Type, fvPatchField, volMesh> > tLaplacian
    (
        fvc::div
        (
            (SfGamma & Sn)*this->tsnGradScheme_().snGrad(vf)
            //+ (SfGammaCorr & grad(vf))
        )
    );

    tLaplacian().rename("laplacian(" + gamma.name() + ',' + vf.name() + ')');

    return tLaplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
