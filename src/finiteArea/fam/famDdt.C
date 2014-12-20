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

#include "areaFields.H"
#include "edgeFields.H"
#include "faMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type> >
ddt
(
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            vf.dimensions()*dimArea/dimTime
        )
    );

    faMatrix<Type>& fam = tfam();

    scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    if (mesh.timeScheme() == faSchemes::SS)
    {}
    else if
    (
        mesh.timeScheme() == faSchemes::EI
     || mesh.timeScheme() == faSchemes::CN
    )
    {
        fam.diag() = rDeltaT*mesh.S();

        if (mesh.moving())
        {
            fam.source() =
                rDeltaT*fam.psi().oldTime().internalField()*mesh.S0();
        }
        else
        {
            fam.source() =
                rDeltaT*fam.psi().oldTime().internalField()*mesh.S();
        }
    }
    else if (mesh.timeScheme() == faSchemes::BD)
    {
        scalar deltaT = mesh().time().deltaT().value();
        scalar deltaT0 = mesh().time().deltaT0().value();

        scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
        scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
        scalar coefft0  = coefft + coefft00;

        fam.diag() = (coefft*rDeltaT)*mesh.S();

        if (mesh.moving())
        {
            fam.source() = rDeltaT*
            (
                coefft0*fam.psi().oldTime().internalField()*mesh.S0()
              - coefft00*fam.psi().oldTime().oldTime().internalField()
               *mesh.S00()
            );
        }
        else
        {
            fam.source() = rDeltaT*mesh.S()*
            (
                coefft0*fam.psi().oldTime().internalField()
              - coefft00*fam.psi().oldTime().oldTime().internalField()
            );
        }
    }
    else
    {
        FatalErrorIn("ddt(GeometricField<Type, faPatchField, areaMesh>&)")
            << "Unsupported temporal differencing scheme : "
            << int(mesh.timeScheme())
            << abort(FatalError);
    }

    return tfam;
}


template<class Type>
tmp<faMatrix<Type> >
ddt
(
    const dimensionedScalar& rho,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );
    faMatrix<Type>& fam = tfam();

    scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    if (mesh.timeScheme() == faSchemes::SS)
    {}
    else if
    (
        mesh.timeScheme() == faSchemes::EI
     || mesh.timeScheme() == faSchemes::CN
    )
    {
        fam.diag() = rDeltaT*rho.value()*mesh.S();

        if (mesh.moving())
        {
            fam.source() = rDeltaT
               *rho.value()*fam.psi().oldTime().internalField()*mesh.S0();
        }
        else
        {
            fam.source() = rDeltaT
               *rho.value()*fam.psi().oldTime().internalField()*mesh.S();
        }
    }
    else if (mesh.timeScheme() == faSchemes::BD)
    {
        scalar deltaT = mesh().time().deltaT().value();
        scalar deltaT0 = mesh().time().deltaT0().value();

        scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
        scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
        scalar coefft0  = coefft + coefft00;

        fam.diag() = (coefft*rDeltaT*rho.value())*mesh.S();

        if (mesh.moving())
        {
            fam.source() = rDeltaT*rho.value()*
            (
                coefft0*fam.psi().oldTime().internalField()*mesh.S0()
              - coefft00*fam.psi().oldTime().oldTime().internalField()
               *mesh.S00()
            );
        }
        else
        {
            fam.source() = rDeltaT*mesh.S()*rho.value()*
            (
                coefft0*fam.psi().oldTime().internalField()
              - coefft00*fam.psi().oldTime().oldTime().internalField()
            );
        }
    }
    else
    {
        FatalErrorIn
        (
            "ddt(const dimesionedScalar&, "
            "GeometricField<Type, faPatchField, areaMesh>&)"
        )   << "Unsupported temporal differencing scheme : "
            << int(mesh.timeScheme())
            << abort(FatalError);
    }

    return tfam;
}


template<class Type>
tmp<faMatrix<Type> >
ddt
(
    const areaScalarField& rho,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );
    faMatrix<Type>& fam = tfam();

    scalar rDeltaT = 1.0/mesh().time().deltaT().value();

    if (mesh.timeScheme() == faSchemes::SS)
    {}
    else if
    (
        mesh.timeScheme() == faSchemes::EI
     || mesh.timeScheme() == faSchemes::CN
    )
    {
        fam.diag() = rDeltaT*rho.internalField()*mesh.S();

        if (mesh.moving())
        {
            fam.source() = rDeltaT
               *rho.oldTime().internalField()
               *fam.psi().oldTime().internalField()*mesh.S0();
        }
        else
        {
            fam.source() = rDeltaT
               *rho.oldTime().internalField()
                *fam.psi().oldTime().internalField()*mesh.S();
        }
    }
    else if (mesh.timeScheme() == faSchemes::BD)
    {
        scalar deltaT = mesh().time().deltaT().value();
        scalar deltaT0 = mesh().time().deltaT0().value();

        scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
        scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
        scalar coefft0  = coefft + coefft00;

        fam.diag() =
            (coefft*rDeltaT)*rho.internalField()*mesh.S();

        if (mesh.moving())
        {
            fam.source() = rDeltaT*
            (
                coefft0*rho.oldTime().internalField()
               *fam.psi().oldTime().internalField()*mesh.S0()
              - coefft00*rho.oldTime().oldTime().internalField()
               *fam.psi().oldTime().oldTime().internalField()*mesh.S00()
            );
        }
        else
        {
            fam.source() = rDeltaT*mesh.S()*
            (
                coefft0*rho.oldTime().internalField()
               *fam.psi().oldTime().internalField()
              - coefft00*rho.oldTime().oldTime().internalField()
               *fam.psi().oldTime().oldTime().internalField()
            );
        }
    }
    else
    {
        FatalErrorIn
        (
            "ddt(areaScalarField&, GeometricField<Type, "
            "faMesh, areaMesh>&)"
        )   << "Unsupported temporal differencing scheme : "
            << int(mesh.timeScheme())
            << abort(FatalError);
    }

    return tfam;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
