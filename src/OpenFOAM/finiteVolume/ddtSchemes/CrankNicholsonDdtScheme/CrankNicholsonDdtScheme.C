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

#include "CrankNicholsonDdtScheme.H"
#include "surfaceInterpolate.H"
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
template<class GeoField>
GeoField& CrankNicholsonDdtScheme<Type>::ddt0_
(
    const word& name,
    const dimensionSet& dims
)
{
    if (!mesh().objectRegistry::foundObject<GeoField>(name))
    {
        scalar t = mesh().time().value();
        scalar t0 = t - mesh().time().deltaT().value();

        IOobject ddt0Header
        (
            name,
            mesh().time().timeName(t0),
            mesh()
        );

        if (ddt0Header.headerOk())
        {
            label ti = mesh().time().timeIndex();
            const_cast<Time&>(mesh().time()).setTime(t0, ti - 1);

            new GeoField
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

            const_cast<Time&>(mesh().time()).setTime(t, ti);
        }
        else
        {
            new GeoField
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
                dimensioned<typename GeoField::value_type>
                (
                    "0",
                    dims/dimTime,
                    pTraits<typename GeoField::value_type>::zero
                )
            );
        }
    }

    return const_cast<GeoField&>
    (
        mesh().objectRegistry::lookupObject<GeoField>(name)
    );
}


template<class Type>
template<class GeoField>
bool CrankNicholsonDdtScheme<Type>::evaluate(const GeoField& ddt0) const
{
    return ddt0.timeIndex() != mesh().time().timeIndex();
}


template<class Type>
template<class GeoField>
dimensionedScalar CrankNicholsonDdtScheme<Type>::rDtCoef_
(
    const GeoField& ddt0
) const
{
    if (mesh().time().timeIndex() == 1 && sumMag(ddt0.internalField()) < VSMALL)
    {
        return 1.0/mesh().time().deltaT();
    }
    else
    {
        return (1 + ocCoeff_)/mesh().time().deltaT();
    }
}


template<class Type>
template<class GeoField>
dimensionedScalar CrankNicholsonDdtScheme<Type>::rDtCoef0_
(
    const GeoField& ddt0
) const
{
    if (mesh().time().timeIndex() == 2 && sumMag(ddt0.internalField()) < VSMALL)
    {
        return 1.0/mesh().time().deltaT0();
    }
    else
    {
        return (1 + ocCoeff_)/mesh().time().deltaT0();
    }
}


template<class Type>
template<class GeoField>
tmp<GeoField> CrankNicholsonDdtScheme<Type>::offCentre_
(
    const GeoField& ddt0
) const
{
    if (ocCoeff_ < 1.0)
    {
        return ocCoeff_*ddt0;
    }
    else
    {
        return ddt0;
    }
}


template<class Type>
const FieldField<fvPatchField, Type>& ff
(
    const FieldField<fvPatchField, Type>& bf
)
{
    return bf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicholsonDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + dt.name() + ')',
            dt.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + dt.name() + ')',
        mesh().time().timeName(),
        mesh()
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

    scalar rDtCoef = rDtCoef_(ddt0).value();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDtCoef0*dt.value()*(mesh().V0() - mesh().V00())
              - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();
        }

        tdtdt().internalField() =
        (
            (rDtCoef*dt.value())*(mesh().V() - mesh().V0())
          - mesh().V0()*offCentre_(ddt0.internalField())
        )/mesh().V();
    }

    return tdtdt;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicholsonDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + vf.name() + ')',
            vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + vf.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(ddt0);

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDtCoef0*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                ) - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDtCoef0*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDtCoef.dimensions()*vf.dimensions(),
                (
                    rDtCoef.value()*
                    (
                        mesh().V()*vf.internalField()
                      - mesh().V0()*vf.oldTime().internalField()
                    ) - mesh().V0()*offCentre_(ddt0.internalField())
                )/mesh().V(),
                rDtCoef.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*(vf.oldTime() - vf.oldTime().oldTime())
                - offCentre_(ddt0);
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDtCoef*(vf - vf.oldTime()) - offCentre_(ddt0)
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicholsonDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<Type, fvPatchField, volMesh>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + rho.name() + ',' + vf.name() + ')',
            rho.dimensions()*vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + rho.name() + ',' + vf.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(ddt0);

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDtCoef0*rho.value()*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                ) - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDtCoef0*rho.value()*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDtCoef.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    rDtCoef.value()*rho.value()*
                    (
                        mesh().V()*vf.internalField()
                      - mesh().V0()*vf.oldTime().internalField()
                    ) - mesh().V0()*offCentre_(ddt0.internalField())
                )/mesh().V(),
                rDtCoef.value()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*rho*(vf.oldTime() - vf.oldTime().oldTime())
                - offCentre_(ddt0);
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDtCoef*rho*(vf - vf.oldTime()) - offCentre_(ddt0)
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
    GeometricField<Type, fvPatchField, volMesh>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + rho.name() + ',' + vf.name() + ')',
            rho.dimensions()*vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + rho.name() + ',' + vf.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(ddt0);

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDtCoef0*
                (
                    mesh().V0()*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - mesh().V00()*rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                ) - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDtCoef0*
                (
                    rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                  - rho.oldTime().oldTime().boundaryField()
                   *vf.oldTime().oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDtCoef.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    rDtCoef.value()*
                    (
                        mesh().V()*rho.internalField()*vf.internalField()
                      - mesh().V0()*rho.oldTime().internalField()
                       *vf.oldTime().internalField()
                    ) - mesh().V00()*offCentre_(ddt0.internalField())
                )/mesh().V(),
                rDtCoef.value()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()*vf.oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            ) - offCentre_(ddt0);
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDtCoef*(rho*vf - rho.oldTime()*vf.oldTime()) - offCentre_(ddt0)
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
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + vf.name() + ')',
            vf.dimensions()
        );

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

    scalar rDtCoef = rDtCoef_(ddt0).value();

    fvm.diag() = rDtCoef*mesh().V();

    vf.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDtCoef0*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                )
                - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDtCoef0*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                )
                - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        fvm.source() =
        (
            rDtCoef*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*(vf.oldTime() - vf.oldTime().oldTime())
                - offCentre_(ddt0);
        }

        fvm.source() = 
        (
            rDtCoef*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
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
    GeometricField<Type, fvPatchField, volMesh>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
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

    scalar rDtCoef = rDtCoef_(ddt0).value();
    fvm.diag() = rDtCoef*rho.value()*mesh().V();

    vf.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDtCoef0*rho.value()*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                )
              - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDtCoef0*rho.value()*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                )
              - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        fvm.source() =
        (
            rDtCoef*rho.value()*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*rho*(vf.oldTime() - vf.oldTime().oldTime())
                - offCentre_(ddt0);
        }

        fvm.source() = 
        (
            rDtCoef*rho.value()*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
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
    GeometricField<Type, fvPatchField, volMesh>& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
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

    scalar rDtCoef = rDtCoef_(ddt0).value();
    fvm.diag() = rDtCoef*rho.internalField()*mesh().V();

    vf.oldTime().oldTime();
    rho.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() = 
            (
                rDtCoef0*
                (
                    mesh().V0()*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - mesh().V00()*rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                )
              - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() = 
            (
                rDtCoef0*
                (
                    rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                  - rho.oldTime().oldTime().boundaryField()
                   *vf.oldTime().oldTime().boundaryField()
                )
              - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        fvm.source() =
        (
            rDtCoef*rho.internalField()*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            ) - offCentre_(ddt0);
        }

        fvm.source() = 
        (
            rDtCoef*rho.oldTime().internalField()*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*mesh().V();
    }

    return tfvm;
}



template<class Type>
tmp<typename CrankNicholsonDdtScheme<Type>::fluxFieldType>
CrankNicholsonDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    GeometricField<Type, fvPatchField, volMesh>& dUdt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + U.name() + ')',
            U.dimensions()
        );

    fluxFieldType& dphidt0 =
        ddt0_<fluxFieldType>
        (
            "ddt0(" + phi.name() + ')',
            phi.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(dUdt0);

    if (mesh().moving())
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                mesh(),
                dimensioned<typename flux<Type>::type>
                (
                    "0",
                    rDtCoef.dimensions()*rA.dimensions()*phi.dimensions(),
                    pTraits<typename flux<Type>::type>::zero
                )
            )
        );
    }
    else
    {
        if (evaluate(dUdt0))
        {
            dUdt0 = 
                rDtCoef0_(dUdt0)*(U.oldTime() - U.oldTime().oldTime())
              - offCentre_(dUdt0);
        }

        if (evaluate(dphidt0))
        {
            dphidt0 = 
                rDtCoef0_(dphidt0)*(phi.oldTime() - phi.oldTime().oldTime())
              - offCentre_(dphidt0);
        }

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())*
                (
                    fvc::interpolate(rA)
                   *(rDtCoef*phi.oldTime() + offCentre_(dphidt0))
                  - (
                        fvc::interpolate
                        (
                            rA*(rDtCoef*U.oldTime() + offCentre_(dUdt0))
                        ) & mesh().Sf()
                    )
                )
            )
        );
    }
}


template<class Type>
tmp<typename CrankNicholsonDdtScheme<Type>::fluxFieldType>
CrankNicholsonDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rA,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    GeometricField<Type, fvPatchField, volMesh>& dUdt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + U.name() + ')',
            rho.dimensions()*U.dimensions()
        );

    fluxFieldType& dphidt0 =
        ddt0_<fluxFieldType>
        (
            "ddt0(" + phi.name() + ')',
            phi.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddtPhiCorr("
      + rA.name() + ',' + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(dUdt0);

    if (mesh().moving())
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                mesh(),
                dimensioned<typename flux<Type>::type>
                (
                    "0",
                    rA.dimensions()*rho.dimensions()*phi.dimensions()/dimTime,
                    pTraits<typename flux<Type>::type>::zero
                )
            )
        );
    }
    else
    {
        if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimVelocity*dimArea
        )
        {
            if (evaluate(dUdt0))
            {
                dUdt0 = rDtCoef0_(dUdt0)*
                (
                    rho.oldTime()*U.oldTime()
                  - rho.oldTime().oldTime()*U.oldTime().oldTime()
                ) - offCentre_(dUdt0);
            }

            if (evaluate(dphidt0))
            {
                dphidt0 = rDtCoef0_(dphidt0)*
                (
                    phi.oldTime()
                  - fvc::interpolate(rho.oldTime().oldTime()/rho.oldTime())
                   *phi.oldTime().oldTime()
                ) - offCentre_(dphidt0);
            }

            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())*
                    (
                        fvc::interpolate(rA*rho.oldTime())
                       *(rDtCoef*phi.oldTime() + offCentre_(dphidt0))
                      - (
                            fvc::interpolate
                            (
                                rA*
                                (
                                    rDtCoef*rho.oldTime()*U.oldTime()
                                  + offCentre_(dUdt0)
                                )
                            ) & mesh().Sf()
                        )
                    )
                )
            );
        }
        else if 
        (
            U.dimensions() == dimDensity*dimVelocity
         && phi.dimensions() == dimDensity*dimVelocity*dimArea
        )
        {
            if (evaluate(dUdt0))
            {
                dUdt0 = rDtCoef0_(dUdt0)*
                (
                    U.oldTime() - U.oldTime().oldTime()
                ) - offCentre_(dUdt0);
            }

            if (evaluate(dphidt0))
            {
                dphidt0 = rDtCoef0_(dphidt0)*
                (
                    phi.oldTime() - phi.oldTime().oldTime()
                ) - offCentre_(dphidt0);
            }

            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    fvcDdtPhiCoeff(rho.oldTime(), U.oldTime(), phi.oldTime())*
                    (
                        fvc::interpolate(rA)
                       *(rDtCoef*phi.oldTime() + offCentre_(dphidt0))
                      - (
                            fvc::interpolate
                            (
                                rA*(rDtCoef*U.oldTime() + offCentre_(dUdt0))
                            ) & mesh().Sf()
                        )
                    )
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "CrankNicholsonDdtScheme<Type>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return fluxFieldType::null();
        }
    }
}


template<class Type>
tmp<surfaceScalarField> CrankNicholsonDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if 
    (
        !mesh().objectRegistry::foundObject
        <surfaceScalarField>("meshPhiCN_0")
    )
    {
        scalar t = mesh().time().value();
        scalar t0 = t - mesh().time().deltaT().value();

        IOobject ddt0Header
        (
            "meshPhiCN_0",
            mesh().time().timeName(t0),
            mesh()
        );

        if (ddt0Header.headerOk())
        {
            label ti = mesh().time().timeIndex();
            const_cast<Time&>(mesh().time()).setTime(t0, ti - 1);

            new surfaceScalarField
            (
                IOobject
                (
                    "meshPhiCN_0",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh()
            );

            const_cast<Time&>(mesh().time()).setTime(t, ti);
        }
        else
        {
            new surfaceScalarField
            (
                IOobject
                (
                    "meshPhiCN_0",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar
                (
                    "0",
                    dimVolume/dimTime,
                    0.0
                )
            );
        }
    }

    surfaceScalarField& meshPhi0 = const_cast<surfaceScalarField&>
    (
        mesh().objectRegistry::lookupObject<surfaceScalarField>("meshPhiCN_0")
    );

    scalar sumMagMeshPhi0 = sumMag(meshPhi0.internalField());

    const fvPatchList& patches = mesh().boundary();

    forAll (patches, patchI)
    {
        sumMagMeshPhi0 += sumMag(meshPhi0.boundaryField()[patchI]);
    }

    scalar C, C0;

    if (mesh().time().timeIndex() == 1 && sumMagMeshPhi0 < VSMALL)
    {
        C = 1;
    }
    else
    {
        C = 1 + ocCoeff_;
    }

    if (mesh().time().timeIndex() == 2 && sumMagMeshPhi0 < VSMALL)
    {
        C0 = 1;
    }
    else
    {
        C0 = 1 + ocCoeff_;
    }

    if (meshPhi0.timeIndex() != mesh().time().timeIndex())
    {
        meshPhi0 = C0*mesh().phi().oldTime() - offCentre_(meshPhi0);
    }

    return C*mesh().phi() - offCentre_(meshPhi0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
