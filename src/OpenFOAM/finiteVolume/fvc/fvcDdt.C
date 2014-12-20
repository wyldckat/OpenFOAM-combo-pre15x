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

#include "fvcDdt.H"
#include "fvMesh.H"
#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
ddt
(
    const dimensioned<Type> dt,
    const fvMesh& mesh
)
{
    return fv::ddtScheme<Type>::New
    (
        mesh,
        mesh.ddtScheme("ddt(" + dt.name() + ')')
    )().fvcDdt(dt);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
ddt0
(
    const dimensioned<Type> dt,
    const fvMesh& mesh
)
{
    return fv::ddtScheme<Type>::New
    (
        mesh,
        mesh.ddtScheme("ddt0(" + dt.name() + ')')
    )().fvcDdt0(dt);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
ddt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt(" + vf.name() + ')')
    )().fvcDdt(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
ddt0
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt0(" + vf.name() + ')')
    )().fvcDdt0(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
ddt0
(
    const GeometricField<Type, fvPatchField, surfaceMesh>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt0(" + vf.name() + ')')
    )().fvcDdt0(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
ddt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt(" + rho.name() + ',' + vf.name() + ')')
    )().fvcDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
ddt0
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt0(" + rho.name() + ',' + vf.name() + ')')
    )().fvcDdt0(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
ddt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt(" + rho.name() + ',' + vf.name() + ')')
    )().fvcDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
ddt0
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt0(" + rho.name() + ',' + vf.name() + ')')
    )().fvcDdt0(rho, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
