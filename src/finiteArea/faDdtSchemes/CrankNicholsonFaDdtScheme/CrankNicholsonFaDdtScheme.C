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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
    
\*---------------------------------------------------------------------------*/

#include "CrankNicholsonFaDdtScheme.H"
#include "facDiv.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
template<class GeoMesh>
GeometricField<Type, faPatchField, GeoMesh>&
CrankNicholsonFaDdtScheme<Type>::ddt0_
(
    const word& name,
    const dimensionSet& dims
)
{
    if 
    (
        !mesh()().objectRegistry::foundObject
        <GeometricField<Type, faPatchField, GeoMesh> >
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
            mesh()().time().timeName(t0),
            mesh()(),
            IOobject::NO_READ
        );

        if (ddt0Header.headerOk())
        {
            label ti = mesh().time().timeIndex();
            const_cast<Time&>(mesh().time()).setTime(t0, ti - 1);

            new GeometricField<Type, faPatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    mesh()().time().timeName(),
                    mesh()(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh()
            );

            const_cast<Time&>(mesh().time()).setTime(t, ti);
        }
        else
        {
            new GeometricField<Type, faPatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    mesh()().time().timeName(),
                    mesh()(),
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

    return const_cast<GeometricField<Type, faPatchField, GeoMesh>&>
    (
        mesh()().objectRegistry::lookupObject
        <GeometricField<Type, faPatchField, GeoMesh> >(name)
    );
}


template<class Type>
template<class GeoMesh>
bool CrankNicholsonFaDdtScheme<Type>::evaluate
(
    const GeometricField<Type, faPatchField, GeoMesh>& ddt0
) const
{
    return ddt0.timeIndex() != mesh().time().timeIndex();
}


template<class Type>
template<class GeoMesh>
dimensionedScalar CrankNicholsonFaDdtScheme<Type>::rDeltaT_
(
    const GeometricField<Type, faPatchField, GeoMesh>& ddt0
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
dimensionedScalar CrankNicholsonFaDdtScheme<Type>::rDeltaT0_
(
    const GeometricField<Type, faPatchField, GeoMesh>& ddt0
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
tmp<GeometricField<Type, faPatchField, areaMesh> >
CrankNicholsonFaDdtScheme<Type>::facDdt
(
    const dimensioned<Type> dt
)
{
    GeometricField<Type, faPatchField, areaMesh>& ddt0 =
        ddt0_<areaMesh>("ddt0(" + dt.name() + ')', dt.dimensions());

    IOobject ddtIOobject
    (
        "ddt(" + dt.name() + ')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    tmp<GeometricField<Type, faPatchField, areaMesh> > tdtdt
    (
        new GeometricField<Type, faPatchField, areaMesh>
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
                rDeltaT0*dt.value()*(mesh().S0() - mesh().S00())
              - mesh().S00()*ddt0.internalField()
            )/mesh().S0();
        }

        tdtdt().internalField() =
        (
            (rDeltaT*dt.value())*(mesh().S() - mesh().S0())
          - mesh().S0()*ddt0.internalField()
        )/mesh().S();
    }

    return tdtdt;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
CrankNicholsonFaDdtScheme<Type>::facDdt0
(
    const dimensioned<Type> dt
)
{
    GeometricField<Type, faPatchField, areaMesh>& ddt0 =
        ddt0_<areaMesh>("ddt0(" + dt.name() + ')', dt.dimensions());

    IOobject ddtIOobject
    (
        "ddt(" + dt.name() + ')',
        mesh()().time().timeName(),
        mesh()(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    dimensionedScalar rDeltaT
    (
        "rDeltaT", dimless/dimTime, rDeltaT_(ddt0).value()
    );

    tmp<GeometricField<Type, faPatchField, areaMesh> > tdtdt0
    (
        new GeometricField<Type, faPatchField, areaMesh>
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
                rDeltaT0*dt.value()*(mesh().S0() - mesh().S00())
              - mesh().S00()*ddt0.internalField()
            )/mesh().S0();
        }

        tdtdt0().internalField() -=
            mesh().S0()*ddt0.internalField()/mesh().S();
    }

    return tdtdt0;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
CrankNicholsonFaDdtScheme<Type>::facDdt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    GeometricField<Type, faPatchField, areaMesh>& ddt0 =
        ddt0_<areaMesh>("ddt0(" + vf.name() + ')', vf.dimensions());

    IOobject ddtIOobject
    (
        "ddt(" + vf.name() + ')',
        mesh()().time().timeName(),
        mesh()(),
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
                    mesh().S0()*vf.oldTime().internalField()
                  - mesh().S00()*vf.oldTime().oldTime().internalField()
                ) - mesh().S00()*ddt0.internalField()
            )/mesh().S0();

            ddt0.boundaryField() = 
            (
                rDeltaT0*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                ) - ddt0.boundaryField()
            );
        }

        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                (
                    rDeltaT.value()*
                    (
                        mesh().S()*vf.internalField()
                      - mesh().S0()*vf.oldTime().internalField()
                    ) - mesh().S0()*ddt0.internalField()
                )/mesh().S(),
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

        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                rDeltaT*(vf - vf.oldTime()) - ddt0
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
CrankNicholsonFaDdtScheme<Type>::facDdt0
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    GeometricField<Type, faPatchField, areaMesh>& ddt0 =
        ddt0_<areaMesh>("ddt0(" + vf.name() + ')', vf.dimensions());

    IOobject ddtIOobject
    (
        "ddtExp(" + vf.name() + ')',
        mesh()().time().timeName(),
        mesh()(),
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
                    mesh().S0()*vf.oldTime().internalField()
                  - mesh().S00()*vf.oldTime().oldTime().internalField()
                ) - mesh().S00()*ddt0.internalField()
            )/mesh().S0();

            ddt0.boundaryField() = 
            (
                rDeltaT0*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                ) - ddt0.boundaryField()
            );
        }

        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                (
                    (-rDeltaT.value())*vf.oldTime().internalField()
                  - ddt0.internalField()
                )*mesh().S0()/mesh().S(),
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

        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                (-rDeltaT)*vf.oldTime() - ddt0
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, edgeMesh> >
CrankNicholsonFaDdtScheme<Type>::facDdt0
(
    const GeometricField<Type, faPatchField, edgeMesh>& vf
)
{
    GeometricField<Type, faPatchField, edgeMesh>& ddt0 =
        ddt0_<edgeMesh>("ddt0(" + vf.name() + ')', vf.dimensions());

    IOobject ddtIOobject
    (
        "ddtExp(" + vf.name() + ')',
        mesh()().time().timeName(),
        mesh()(),
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

    return tmp<GeometricField<Type, faPatchField, edgeMesh> >
    (
        new GeometricField<Type, faPatchField, edgeMesh>
        (
            ddtIOobject,
            (-rDeltaT)*vf.oldTime() - ddt0
        )
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
CrankNicholsonFaDdtScheme<Type>::facDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    GeometricField<Type, faPatchField, areaMesh>& ddt0 = ddt0_<areaMesh>
    (
        "ddt0(" + rho.name() + ',' + vf.name() + ')',
        rho.dimensions()*vf.dimensions()
    );

    IOobject ddtIOobject
    (
        "ddt(" + rho.name() + ',' + vf.name() + ')',
        mesh()().time().timeName(),
        mesh()(),
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
                    mesh().S0()*vf.oldTime().internalField()
                  - mesh().S00()*vf.oldTime().oldTime().internalField()
                ) - mesh().S00()*ddt0.internalField()
            )/mesh().S0();

            ddt0.boundaryField() = 
            (
                rDeltaT0*rho.value()*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                ) - ddt0.boundaryField()
            );
        }

        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    rDeltaT.value()*rho.value()*
                    (
                        mesh().S()*vf.internalField()
                      - mesh().S0()*vf.oldTime().internalField()
                    ) - mesh().S0()*ddt0.internalField()
                )/mesh().S(),
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

        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                rDeltaT*rho*(vf - vf.oldTime()) - ddt0
            )
        );
    }
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
CrankNicholsonFaDdtScheme<Type>::facDdt0
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    GeometricField<Type, faPatchField, areaMesh>& ddt0 = ddt0_<areaMesh>
    (
        "ddt0(" + rho.name() + ',' + vf.name() + ')',
        rho.dimensions()*vf.dimensions()
    );

    IOobject ddtIOobject
    (
        "ddtExp(" + rho.name() + ',' + vf.name() + ')',
        mesh()().time().timeName(),
        mesh()(),
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
               *mesh().S00()*vf.oldTime().oldTime().internalField()
              - mesh().S00()*ddt0.internalField()
            )/mesh().S0();

            ddt0.boundaryField() = 
            (
                (-rDeltaT0*rho.value())*vf.oldTime().oldTime().boundaryField()
              - ddt0.boundaryField()
            );
        }

        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    (-rDeltaT.value()*rho.value())*vf.oldTime().internalField()
                  - ddt0.internalField()
                )*mesh().S0()/mesh().S(),
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

        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                (-rDeltaT*rho)*vf.oldTime() - ddt0
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
CrankNicholsonFaDdtScheme<Type>::facDdt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    GeometricField<Type, faPatchField, areaMesh>& ddt0 = ddt0_<areaMesh>
    (
        "ddt0(" + rho.name() + ',' + vf.name() + ')',
        rho.dimensions()*vf.dimensions()
    );

    IOobject ddtIOobject
    (
        "ddt(" + rho.name() + ',' + vf.name() + ')',
        mesh()().time().timeName(),
        mesh()(),
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
                    mesh().S0()*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - mesh().S00()*rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                ) - mesh().S00()*ddt0.internalField()
            )/mesh().S0();

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

        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    rDeltaT.value()*
                    (
                        mesh().S()*rho.internalField()*vf.internalField()
                      - mesh().S0()*rho.oldTime().internalField()
                       *vf.oldTime().internalField()
                    ) - mesh().S00()*ddt0.internalField()
                )/mesh().S(),
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

        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime()) - ddt0
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
CrankNicholsonFaDdtScheme<Type>::facDdt0
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    GeometricField<Type, faPatchField, areaMesh>& ddt0 = ddt0_<areaMesh>
    (
        "ddt0(" + rho.name() + ',' + vf.name() + ')',
        rho.dimensions()*vf.dimensions()
    );

    IOobject ddtIOobject
    (
        "ddtExp(" + rho.name() + ',' + vf.name() + ')',
        mesh()().time().timeName(),
        mesh()(),
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
                    mesh().S0()*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - mesh().S00()*rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                )
              - mesh().S00()*ddt0.internalField()
            )/mesh().S0();

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

        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    (-rDeltaT.value())*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - ddt0.internalField()
                )*mesh().S0()/mesh().S(),
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

        return tmp<GeometricField<Type, faPatchField, areaMesh> >
        (
            new GeometricField<Type, faPatchField, areaMesh>
            (
                ddtIOobject,
                (-rDeltaT)*rho.oldTime()*vf.oldTime() - ddt0
            )
        );
    }
}


template<class Type>
tmp<faMatrix<Type> >
CrankNicholsonFaDdtScheme<Type>::famDdt
(
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    GeometricField<Type, faPatchField, areaMesh>& ddt0 =
        ddt0_<areaMesh>("ddt0(" + vf.name() + ')', vf.dimensions());

    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            vf.dimensions()*dimArea/dimTime
        )
    );

    faMatrix<Type>& fam = tfam();

    scalar rDeltaT = rDeltaT_(ddt0).value();

    fam.diag() = rDeltaT*mesh().S();

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
                    mesh().S0()*vf.oldTime().internalField()
                  - mesh().S00()*vf.oldTime().oldTime().internalField()
                )
              - mesh().S00()*ddt0.internalField()
            )/mesh().S0();

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

        fam.source() =
        (
            rDeltaT*vf.oldTime().internalField() + ddt0.internalField()
        )*mesh().S0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDeltaT0_(ddt0)*(vf.oldTime() - vf.oldTime().oldTime())
                - ddt0;
        }

        fam.source() = 
        (
            rDeltaT*vf.oldTime().internalField() + ddt0.internalField()
        )*mesh().S();
    }

    return tfam;
}


template<class Type>
tmp<faMatrix<Type> >
CrankNicholsonFaDdtScheme<Type>::famDdt
(
    const dimensionedScalar& rho,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    GeometricField<Type, faPatchField, areaMesh>& ddt0 = ddt0_<areaMesh>
    (
        "ddt0(" + rho.name() + ',' + vf.name() + ')',
        rho.dimensions()*vf.dimensions()
    );

    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );
    faMatrix<Type>& fam = tfam();

    scalar rDeltaT = rDeltaT_(ddt0).value();
    fam.diag() = rDeltaT*rho.value()*mesh().S();

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
                    mesh().S0()*vf.oldTime().internalField()
                  - mesh().S00()*vf.oldTime().oldTime().internalField()
                )
              - mesh().S00()*ddt0.internalField()
            )/mesh().S0();

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

        fam.source() =
        (
            rDeltaT*rho.value()*vf.oldTime().internalField()
          + ddt0.internalField()
        )*mesh().S0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDeltaT0_(ddt0)*rho*(vf.oldTime() - vf.oldTime().oldTime())
                - ddt0;
        }

        fam.source() = 
        (
            rDeltaT*rho.value()*vf.oldTime().internalField()
          + ddt0.internalField()
        )*mesh().S();
    }

    return tfam;
}


template<class Type>
tmp<faMatrix<Type> >
CrankNicholsonFaDdtScheme<Type>::famDdt
(
    const areaScalarField& rho,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    GeometricField<Type, faPatchField, areaMesh>& ddt0 = ddt0_<areaMesh>
    (
        "ddt0(" + rho.name() + ',' + vf.name() + ')',
        rho.dimensions()*vf.dimensions()
    );

    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );
    faMatrix<Type>& fam = tfam();

    scalar rDeltaT = rDeltaT_(ddt0).value();
    fam.diag() = rDeltaT*rho.internalField()*mesh().S();

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
                    mesh().S0()*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - mesh().S00()*rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                )
              - mesh().S00()*ddt0.internalField()
            )/mesh().S0();

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

        fam.source() =
        (
            rDeltaT*rho.internalField()*vf.oldTime().internalField()
          + ddt0.internalField()
        )*mesh().S0();
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

        fam.source() = 
        (
            rDeltaT*rho.oldTime().internalField()*vf.oldTime().internalField()
          + ddt0.internalField()
        )*mesh().S();
    }

    return tfam;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
