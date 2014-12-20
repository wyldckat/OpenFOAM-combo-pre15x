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

#include "gaussFaConvectionScheme.H"
#include "facEdgeIntegrate.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, edgeMesh> >
gaussConvectionScheme<Type>::flux
(
    const edgeScalarField& faceFlux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    return faceFlux*tinterpScheme_().interpolate(vf);
}


template<class Type>
tmp<faMatrix<Type> >
gaussConvectionScheme<Type>::famDiv
(
    const edgeScalarField& faceFlux,
    GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    tmp<edgeScalarField> tweights = tinterpScheme_().weights(vf);
    const edgeScalarField& weights = tweights();

    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam();

    fam.lower() = -weights.internalField()*faceFlux.internalField();
    fam.upper() = fam.lower() + faceFlux.internalField();
    fam.negSumDiag();

    forAll(fam.psi().boundaryField(), patchI)
    {
        const faPatchField<Type>& psf = fam.psi().boundaryField()[patchI];
        const faPatchScalarField& patchFlux = faceFlux.boundaryField()[patchI];
        const faPatchScalarField& pw = weights.boundaryField()[patchI];

        fam.internalCoeffs()[patchI] = patchFlux*psf.valueInternalCoeffs(pw);
        fam.boundaryCoeffs()[patchI] = -patchFlux*psf.valueBoundaryCoeffs(pw);
    }

    if (this->mesh().timeScheme() == faSchemes::CN)
    {
        fam *= 0.5;
        fam.source() = -(fam & vf.oldTime().internalField());
    }

    if (tinterpScheme_().corrected())
    {
        fam += fac::edgeIntegrate(faceFlux*tinterpScheme_().correction(vf));
    }

    return tfam;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
gaussConvectionScheme<Type>::facDiv
(
    const edgeScalarField& faceFlux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > tConvection
    (
        fac::edgeIntegrate(flux(faceFlux, vf))
    );

    tConvection().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );

    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
