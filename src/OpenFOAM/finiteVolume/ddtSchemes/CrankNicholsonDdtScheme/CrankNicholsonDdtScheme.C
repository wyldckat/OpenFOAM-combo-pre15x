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
    
\*---------------------------------------------------------------------------*/

#include "CrankNicholsonDdtScheme.H"
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
template<class GeoMesh>
GeometricField<Type, fvPatchField, GeoMesh>&
CrankNicholsonDdtScheme<Type>::ddt0_
(
    const word& name,
    const dimensionSet& dims
)
{
    if 
    (
        !mesh().objectRegistry::foundObject
        <GeometricField<Type, fvPatchField, GeoMesh> >
        (
            name
        )
    )
    {
        scalar t = mesh().time().value();
        scalar t0 = t - mesh().time().deltaT().value();

        IOobject ddt0Header
        (
            name,
            mesh().time().timeName(t0),
            mesh(),
            IOobject::NO_READ
        );

        if (ddt0Header.headerOk())
        {
            label ti = mesh().time().timeIndex();
            ((Time&)mesh().time()).setTime(t0, ti - 1);

            new GeometricField<Type, fvPatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh()
            );

            ((Time&)mesh().time()).setTime(t, ti);
        }
        else
        {
            new GeometricField<Type, fvPatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dims/dimTime,
                    pTraits<Type>::zero
                )
            );
        }
    }

    return 
        (GeometricField<Type, fvPatchField, GeoMesh>&)
        mesh().objectRegistry::lookupObject
        <GeometricField<Type, fvPatchField, GeoMesh> >(name);
}


template<class Type>
template<class GeoMesh>
bool CrankNicholsonDdtScheme<Type>::evaluate
(
    const GeometricField<Type, fvPatchField, GeoMesh>& ddt0
) const
{
    return ddt0.timeIndex() != mesh().time().timeIndex();
}


template<class Type>
template<class GeoMesh>
dimensionedScalar CrankNicholsonDdtScheme<Type>::rDeltaT_
(
    const GeometricField<Type, fvPatchField, GeoMesh>& ddt0
) const
{
    if (mesh().time().timeIndex() == 1 && sumMag(ddt0.internalField()) < VSMALL)
    {
        return 1.0/mesh().time().deltaT();
    }
    else
    {
        return 2.0/mesh().time().deltaT();
    }
}


template<class Type>
template<class GeoMesh>
dimensionedScalar CrankNicholsonDdtScheme<Type>::rDeltaT0_
(
    const GeometricField<Type, fvPatchField, GeoMesh>& ddt0
) const
{
    if (mesh().time().timeIndex() == 2 && sumMag(ddt0.internalField()) < VSMALL)
    {
        return 1.0/mesh().time().deltaT0();
    }
    else
    {
        return 2.0/mesh().time().deltaT0();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicholsonDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type> dt
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 =
        ddt0_<volMesh>("ddt0(" + dt.name() + ')', dt.dimensions());

    IOobject ddtIOobject
    (
        "ddt(" + dt.name() + ')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    tmp<GeometricField<Type, fvPatchField, volMesh> > tdtdt
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            ddtIOobject,
            mesh(),
            dimensioned<Type>
            (
                "0",
                dt.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );

    scalar rDeltaT = rDeltaT_(ddt0).value();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDeltaT0 = rDeltaT0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDeltaT0*dt.value()*(mesh().V0() - mesh().V00())
              - mesh().V00()*ddt0.internalField()
            )/mesh().V0();
        }

        tdtdt().internalField() =
        (
            (rDeltaT*dt.value())*(mesh().V() - mesh().V0())
          - mesh().V0()*ddt0.internalField()
        )/mesh().V();
    }

    return tdtdt;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicholsonDdtScheme<Type>::fvcDdt0
(
    const dimensioned<Type> dt
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 =
        ddt0_<volMesh>("ddt0(" + dt.name() + ')', dt.dimensions());

    IOobject ddtIOobject
    (
        "ddt(" + dt.name() + ')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    dimensionedScalar rDeltaT
    (
        "rDeltaT", dimless/dimTime, rDeltaT_(ddt0).value()
    );

    tmp<GeometricField<Type, fvPatchField, volMesh> > tdtdt0
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            ddtIOobject,
            mesh(),
            -rDeltaT*dt
        )
    );

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDeltaT0 = rDeltaT0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDeltaT0*dt.value()*(mesh().V0() - mesh().V00())
              - mesh().V00()*ddt0.internalField()
            )/mesh().V0();
        }

        tdtdt0().internalField() -=
            mesh().V0()*ddt0.internalField()/mesh().V();
    }

    return tdtdt0;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicholsonDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 =
        ddt0_<volMesh>("ddt0(" + vf.name() + ')', vf.dimensions());

    IOobject ddtIOobject
    (
        "ddt(" + vf.name() + ')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    dimensionedScalar rDeltaT
    (
        "rDeltaT", dimless/dimTime, rDeltaT_(ddt0).value()
    );

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDeltaT0 = rDeltaT0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDeltaT0*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                ) - mesh().V00()*ddt0.internalField()
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDeltaT0*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                ) - ddt0.boundaryField()
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                (
                    rDeltaT.value()*
                    (
                        mesh().V()*vf.internalField()
                      - mesh().V0()*vf.oldTime().internalField()
                    ) - mesh().V0()*ddt0.internalField()
                )/mesh().V(),
                rDeltaT.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                ) - ddt0.boundaryField()
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDeltaT0_(ddt0)*(vf.oldTime() - vf.oldTime().oldTime())
                - ddt0;
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(vf - vf.oldTime()) - ddt0
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicholsonDdtScheme<Type>::fvcDdt0
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 =
        ddt0_<volMesh>("ddt0(" + vf.name() + ')', vf.dimensions());

    IOobject ddtIOobject
    (
        "ddtExp(" + vf.name() + ')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    dimensionedScalar rDeltaT
    (
        "rDeltaT", dimless/dimTime, rDeltaT_(ddt0).value()
    );

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDeltaT0 = rDeltaT0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDeltaT0*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                ) - mesh().V00()*ddt0.internalField()
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDeltaT0*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                ) - ddt0.boundaryField()
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                (
                    (-rDeltaT.value())*vf.oldTime().internalField()
                  - ddt0.internalField()
                )*mesh().V0()/mesh().V(),
                (-rDeltaT.value())*vf.oldTime().boundaryField()
                  - ddt0.boundaryField()
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDeltaT0_(ddt0)*(vf.oldTime() - vf.oldTime().oldTime())
                - ddt0;
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                (-rDeltaT)*vf.oldTime() - ddt0
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
CrankNicholsonDdtScheme<Type>::fvcDdt0
(
    const GeometricField<Type, fvPatchField, surfaceMesh>& vf
)
{
    GeometricField<Type, fvPatchField, surfaceMesh>& ddt0 =
        ddt0_<surfaceMesh>("ddt0(" + vf.name() + ')', vf.dimensions());

    IOobject ddtIOobject
    (
        "ddtExp(" + vf.name() + ')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    dimensionedScalar rDeltaT
    (
        "rDeltaT", dimless/dimTime, rDeltaT_(ddt0).value()
    );

    if (evaluate(ddt0))
    {
        ddt0 = rDeltaT0_(ddt0)*(vf.oldTime() - vf.oldTime().oldTime()) - ddt0;
    }

    return tmp<GeometricField<Type, fvPatchField, surfaceMesh> >
    (
        new GeometricField<Type, fvPatchField, surfaceMesh>
        (
            ddtIOobject,
            (-rDeltaT)*vf.oldTime() - ddt0
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicholsonDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 = ddt0_<volMesh>
    (
        "ddt0(" + rho.name() + ',' + vf.name() + ')',
        rho.dimensions()*vf.dimensions()
    );

    IOobject ddtIOobject
    (
        "ddt(" + rho.name() + ',' + vf.name() + ')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    dimensionedScalar rDeltaT
    (
        "rDeltaT", dimless/dimTime, rDeltaT_(ddt0).value()
    );

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDeltaT0 = rDeltaT0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDeltaT0*rho.value()*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                ) - mesh().V00()*ddt0.internalField()
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDeltaT0*rho.value()*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                ) - ddt0.boundaryField()
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    rDeltaT.value()*rho.value()*
                    (
                        mesh().V()*vf.internalField()
                      - mesh().V0()*vf.oldTime().internalField()
                    ) - mesh().V0()*ddt0.internalField()
                )/mesh().V(),
                rDeltaT.value()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                ) - ddt0.boundaryField()
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDeltaT0_(ddt0)*rho*(vf.oldTime() - vf.oldTime().oldTime())
                - ddt0;
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*(vf - vf.oldTime()) - ddt0
            )
        );
    }
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicholsonDdtScheme<Type>::fvcDdt0
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 = ddt0_<volMesh>
    (
        "ddt0(" + rho.name() + ',' + vf.name() + ')',
        rho.dimensions()*vf.dimensions()
    );

    IOobject ddtIOobject
    (
        "ddtExp(" + rho.name() + ',' + vf.name() + ')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    dimensionedScalar rDeltaT
    (
        "rDeltaT", dimless/dimTime, rDeltaT_(ddt0).value()
    );

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDeltaT0 = rDeltaT0_(ddt0).value();

            ddt0.internalField() = 
            (
                (-rDeltaT0*rho.value())
               *mesh().V00()*vf.oldTime().oldTime().internalField()
              - mesh().V00()*ddt0.internalField()
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                (-rDeltaT0*rho.value())*vf.oldTime().oldTime().boundaryField()
              - ddt0.boundaryField()
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    (-rDeltaT.value()*rho.value())*vf.oldTime().internalField()
                  - ddt0.internalField()
                )*mesh().V0()/mesh().V(),
                (-rDeltaT.value()*rho.value())*vf.oldTime().boundaryField()
                  - ddt0.boundaryField()
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDeltaT0_(ddt0)*rho*(vf.oldTime() - vf.oldTime().oldTime())
                - ddt0;
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                (-rDeltaT*rho)*vf.oldTime() - ddt0
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicholsonDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 = ddt0_<volMesh>
    (
        "ddt0(" + rho.name() + ',' + vf.name() + ')',
        rho.dimensions()*vf.dimensions()
    );

    IOobject ddtIOobject
    (
        "ddt(" + rho.name() + ',' + vf.name() + ')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    dimensionedScalar rDeltaT
    (
        "rDeltaT", dimless/dimTime, rDeltaT_(ddt0).value()
    );

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDeltaT0 = rDeltaT0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDeltaT0*
                (
                    mesh().V0()*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - mesh().V00()*rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                ) - mesh().V00()*ddt0.internalField()
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDeltaT0*
                (
                    rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                  - rho.oldTime().oldTime().boundaryField()
                   *vf.oldTime().oldTime().boundaryField()
                ) - ddt0.boundaryField()
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    rDeltaT.value()*
                    (
                        mesh().V()*rho.internalField()*vf.internalField()
                      - mesh().V0()*rho.oldTime().internalField()
                       *vf.oldTime().internalField()
                    ) - mesh().V00()*ddt0.internalField()
                )/mesh().V(),
                rDeltaT.value()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()*vf.oldTime().boundaryField()
                ) - ddt0.boundaryField()
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDeltaT0_(ddt0)*
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            ) - ddt0;
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime()) - ddt0
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicholsonDdtScheme<Type>::fvcDdt0
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 = ddt0_<volMesh>
    (
        "ddt0(" + rho.name() + ',' + vf.name() + ')',
        rho.dimensions()*vf.dimensions()
    );

    IOobject ddtIOobject
    (
        "ddtExp(" + rho.name() + ',' + vf.name() + ')',
        mesh().time().timeName(),
        mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    dimensionedScalar rDeltaT
    (
        "rDeltaT", dimless/dimTime, rDeltaT_(ddt0).value()
    );

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDeltaT0 = rDeltaT0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDeltaT0*
                (
                    mesh().V0()*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - mesh().V00()*rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                )
              - mesh().V00()*ddt0.internalField()
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDeltaT0*
                (
                    rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                  - rho.oldTime().oldTime().boundaryField()
                   *vf.oldTime().oldTime().boundaryField()
                )
              - ddt0.boundaryField()
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    (-rDeltaT.value())*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - ddt0.internalField()
                )*mesh().V0()/mesh().V(),
                (-rDeltaT.value())*rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                 - ddt0.boundaryField()
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDeltaT0_(ddt0)*
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            ) - ddt0;
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                (-rDeltaT)*rho.oldTime()*vf.oldTime() - ddt0
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type> >
CrankNicholsonDdtScheme<Type>::fvmDdt
(
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 =
        ddt0_<volMesh>("ddt0(" + vf.name() + ')', vf.dimensions());

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

    scalar rDeltaT = rDeltaT_(ddt0).value();

    fvm.diag() = rDeltaT*mesh().V();

    vf.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDeltaT0 = rDeltaT0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDeltaT0*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                )
              - mesh().V00()*ddt0.internalField()
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDeltaT0*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                )
              - ddt0.boundaryField()
            );
        }

        fvm.source() =
        (
            rDeltaT*vf.oldTime().internalField() + ddt0.internalField()
        )*mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDeltaT0_(ddt0)*(vf.oldTime() - vf.oldTime().oldTime())
                - ddt0;
        }

        fvm.source() = 
        (
            rDeltaT*vf.oldTime().internalField() + ddt0.internalField()
        )*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
CrankNicholsonDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 = ddt0_<volMesh>
    (
        "ddt0(" + rho.name() + ',' + vf.name() + ')',
        rho.dimensions()*vf.dimensions()
    );

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalar rDeltaT = rDeltaT_(ddt0).value();
    fvm.diag() = rDeltaT*rho.value()*mesh().V();

    vf.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDeltaT0 = rDeltaT0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDeltaT0*rho.value()*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                )
              - mesh().V00()*ddt0.internalField()
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDeltaT0*rho.value()*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                )
              - ddt0.boundaryField()
            );
        }

        fvm.source() =
        (
            rDeltaT*rho.value()*vf.oldTime().internalField()
          + ddt0.internalField()
        )*mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDeltaT0_(ddt0)*rho*(vf.oldTime() - vf.oldTime().oldTime())
                - ddt0;
        }

        fvm.source() = 
        (
            rDeltaT*rho.value()*vf.oldTime().internalField()
          + ddt0.internalField()
        )*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
CrankNicholsonDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 = ddt0_<volMesh>
    (
        "ddt0(" + rho.name() + ',' + vf.name() + ')',
        rho.dimensions()*vf.dimensions()
    );

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalar rDeltaT = rDeltaT_(ddt0).value();
    fvm.diag() = rDeltaT*rho.internalField()*mesh().V();

    vf.oldTime().oldTime();
    rho.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDeltaT0 = rDeltaT0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDeltaT0*
                (
                    mesh().V0()*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - mesh().V00()*rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                )
              - mesh().V00()*ddt0.internalField()
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDeltaT0*
                (
                    rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                  - rho.oldTime().oldTime().boundaryField()
                   *vf.oldTime().oldTime().boundaryField()
                )
              - ddt0.boundaryField()
            );
        }

        fvm.source() =
        (
            rDeltaT*rho.internalField()*vf.oldTime().internalField()
          + ddt0.internalField()
        )*mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDeltaT0_(ddt0)*
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            ) - ddt0;
        }

        fvm.source() = 
        (
            rDeltaT*rho.oldTime().internalField()*vf.oldTime().internalField()
          + ddt0.internalField()
        )*mesh().V();
    }

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
