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
    Finite Volume matrix

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "calculatedFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class Type2>
void fvMatrix<Type>::addToInternalField
(
    const unallocLabelList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorIn
        (
            "fvMatrix<Type>::addToInternalField(const unallocLabelList&, "
            "const Field&, Field&)"
        )   << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] += pf[faceI];
    }
}


template<class Type>
template<class Type2>
void fvMatrix<Type>::addToInternalField
(
    const unallocLabelList& addr,
    const tmp<Field<Type2> >& tpf,
    Field<Type2>& intf
) const
{
    addToInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
template<class Type2>
void fvMatrix<Type>::subtractFromInternalField
(
    const unallocLabelList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorIn
        (
            "fvMatrix<Type>::addToInternalField(const unallocLabelList&, "
            "const Field&, Field&)"
        )   << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] -= pf[faceI];
    }
}


template<class Type>
template<class Type2>
void fvMatrix<Type>::subtractFromInternalField
(
    const unallocLabelList& addr,
    const tmp<Field<Type2> >& tpf,
    Field<Type2>& intf
) const
{
    subtractFromInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
void fvMatrix<Type>::addBoundaryDiag
(
    scalarField& diag,
    const direction solveCmpt
) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchI),
            internalCoeffs_[patchI].component(solveCmpt),
            diag
        );
    }
}


template<class Type>
void fvMatrix<Type>::addCmptAvBoundaryDiag(scalarField& diag) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchI),
            cmptAv(internalCoeffs_[patchI]),
            diag
        );
    }
}


template<class Type>
void fvMatrix<Type>::addBoundarySource(Field<Type>& source) const
{
    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];
        const Field<Type>& pbc = boundaryCoeffs_[patchI];

        if (ptf.coupled())
        {
            tmp<Field<Type> > tpnf = ptf.patchNeighbourField();
            const Field<Type>& pnf = tpnf();

            const unallocLabelList& addr = lduAddr().patchAddr(patchI);

            forAll(addr, facei)
            {
                source[addr[facei]] += scale(pbc[facei], pnf[facei]);
            }
        }
        else
        {
            addToInternalField(lduAddr().patchAddr(patchI), pbc, source);
        }
    }
}


template<class Type>
void fvMatrix<Type>::addBoundarySource
(
    scalarField& source,
    const direction solveCmpt
) const
{
    forAll(psi_.boundaryField(), patchI)
    {
        if (!psi_.boundaryField()[patchI].coupled())
        {
            addToInternalField
            (
                lduAddr().patchAddr(patchI),
                boundaryCoeffs_[patchI].component(solveCmpt),
                source
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fvMatrix<Type>::fvMatrix
(
    GeometricField<Type, fvPatchField, volMesh>& psi,
    const dimensionSet& ds
)
:
    lduMatrix(psi.mesh().ldu(), psi.mesh().parallelData().patchSchedule()),
    psi_(psi),
    dimensions_(ds),
    source_(psi.size(), pTraits<Type>::zero),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(NULL),
    solvingComponent(0)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>(GeometricField<Type, fvPatchField, volMesh>&,"
               " const dimensionSet&) : "
               "constructing fvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    // Initialise coupling coefficients
    forAll (psi.mesh().boundary(), patchI)
    {
        internalCoeffs_.hook
        (
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );

        boundaryCoeffs_.hook
        (
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );
    }

    psi_.boundaryField().updateCoeffs();
}


template<class Type>
fvMatrix<Type>::fvMatrix(const fvMatrix<Type>& fvm)
:
    refCount(),
    lduMatrix(fvm),
    psi_(fvm.psi_),
    dimensions_(fvm.dimensions_),
    source_(fvm.source_),
    internalCoeffs_(fvm.internalCoeffs_),
    boundaryCoeffs_(fvm.boundaryCoeffs_),
    faceFluxCorrectionPtr_(NULL),
    solvingComponent(fvm.solvingComponent)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>::fvMatrix(const fvMatrix<Type>&) : "
            << "copying fvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (fvm.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, fvPatchField, surfaceMesh>
        (
            *(fvm.faceFluxCorrectionPtr_)
        );
    }
}


template<class Type>
fvMatrix<Type>::fvMatrix
(
    GeometricField<Type, fvPatchField, volMesh>& psi,
    Istream& is
)
:
    lduMatrix(psi.mesh().ldu(), psi.mesh().parallelData().patchSchedule()),
    psi_(psi),
    dimensions_(is),
    source_(is),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(NULL),
    solvingComponent(0)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>(GeometricField<Type, fvPatchField, volMesh>&,"
               " Istream&) : "
               "constructing fvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    // Initialise coupling coefficients
    forAll (psi.mesh().boundary(), patchI)
    {
        internalCoeffs_.hook
        (
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );

        boundaryCoeffs_.hook
        (
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                pTraits<Type>::zero
            )
        );
    }

}


template<class Type>
fvMatrix<Type>::~fvMatrix()
{
    if (debug)
    {
        Info<< "fvMatrix<Type>::~fvMatrix<Type>() : "
            << "destroying fvMatrix<Type> for field " << psi_.name()
            << endl;
    }

    if (faceFluxCorrectionPtr_)
    {
        delete faceFluxCorrectionPtr_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set reference level for solution
template<class Type>
typename fvMatrix<Type>::reference fvMatrix<Type>::setReference
(
    const label cell,
    const Type& value
)
{
    if (psi_.needReference())
    {
        if (Pstream::master())
        {
            //source()[cell] += diag()[cell]*value;
            //diag()[cell] += diag()[cell];

            autoPtr<constraint<Type> > refPtr
            (
                new constraint<Type>(*this, cell, value)
            );

            refPtr->eliminateEquation(*this);
            refPtr->setSource(*this);

            return refPtr;
        }
    }

    return autoPtr<constraint<Type> >(NULL);
}


// Unset reference level for solution
template<class Type>
void fvMatrix<Type>::unsetReference(const reference& ref)
{
    if (psi_.needReference())
    {
        if (Pstream::master())
        {
            ref->reconstructMatrix(*this);
        }
    }
}


template<class Type>
void fvMatrix<Type>::relax(const scalar alpha)
{
    if (alpha <= 0)
    {
        return;
    }

    {
        FieldField<Field, scalar> bouCoeffsCmptAvgMag
        (
            cmptAv(cmptMag(boundaryCoeffs_))
        );

        lduCoupledInterfacePtrsList interfaces(psi_.boundaryField().size());

        forAll (bouCoeffsCmptAvgMag, patchI)
        {
            interfaces[patchI] = &psi_.boundaryField()[patchI];
        }

        scalarField oldDiag(diag());
        lduMatrix::relax(bouCoeffsCmptAvgMag, interfaces, alpha);
        source() += (diag() - oldDiag)*psi_.internalField();
    }

    forAll(psi_.boundaryField(), patchI)
    {
        const fvPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (!ptf.coupled())
        {
            Field<Type> oldBoundaryDiag = internalCoeffs_[patchI];
            internalCoeffs_[patchI] /= alpha;

            addToInternalField
            (
                lduAddr().patchAddr(patchI),
                scale
                (
                    internalCoeffs_[patchI] - oldBoundaryDiag,
                    ptf.patchInternalField()
                ),
                source()
            );
        }
    }
}


template<class Type>
void fvMatrix<Type>::relax()
{
    if (psi_.mesh().relax(psi_.name()))
    {
        relax(psi_.mesh().relaxationFactor(psi_.name()));
    }
}


template<class Type>
tmp<scalarField> fvMatrix<Type>::D() const
{
    tmp<scalarField> tdiag(new scalarField(diag()));
    addCmptAvBoundaryDiag(tdiag());
    return tdiag;
}


template<class Type>
tmp<volScalarField> fvMatrix<Type>::A() const
{
    tmp<volScalarField> tAphi
    (
        new volScalarField
        (
            IOobject
            (
                "A("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/psi_.dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tAphi().internalField() = D()/psi_.mesh().V();
    tAphi().correctBoundaryConditions();

    return tAphi;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > fvMatrix<Type>::H() const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tHphi
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        scalarField psiCmpt = psi_.internalField().component(cmpt);

        scalarField boundaryDiagCmpt(psi_.size(), 0.0);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);
        boundaryDiagCmpt.negate();
        addCmptAvBoundaryDiag(boundaryDiagCmpt);

        tHphi().internalField().replace(cmpt, boundaryDiagCmpt*psiCmpt);
    }

    tHphi().internalField() += lduMatrix::H(psi_.internalField()) + source_;
    addBoundarySource(tHphi().internalField());

    tHphi().internalField() /= psi_.mesh().V();
    tHphi().correctBoundaryConditions();

    return tHphi;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, surfaceMesh> > fvMatrix<Type>::
flux() const
{
    if (!psi_.mesh().fluxRequired(psi_.name()))
    {
        FatalErrorIn("fvMatrix<Type>::flux()")
            << "flux requested but " << psi_.name()
            << " not specified in the fluxRequired sub-dictionary of fvSchemes."
            << abort(FatalError);
    }

    // construct GeometricField<Type, fvPatchField, surfaceMesh>
    tmp<GeometricField<Type, fvPatchField, surfaceMesh> > tfieldFlux
    (
        new GeometricField<Type, fvPatchField, surfaceMesh>
        (
            IOobject
            (
                "flux("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions()
        )
    );
    GeometricField<Type, fvPatchField, surfaceMesh>& fieldFlux = tfieldFlux();

    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
    {
        fieldFlux.internalField().replace
        (
            cmpt,
            lduMatrix::faceH(psi_.internalField().component(cmpt))
        );
    }

    FieldField<Field, Type> InternalContrib = internalCoeffs_;

    forAll(InternalContrib, patchI)
    {
        InternalContrib[patchI] =
            scale
            (
                InternalContrib[patchI],
                psi_.boundaryField()[patchI].patchInternalField()
            );
    }

    FieldField<Field, Type> NeighbourContrib = boundaryCoeffs_;

    forAll(NeighbourContrib, patchI)
    {
        if (psi_.boundaryField()[patchI].coupled())
        {
            NeighbourContrib[patchI] =
                scale
                (
                    NeighbourContrib[patchI],
                    psi_.boundaryField()[patchI].patchNeighbourField()
                );
        }
    }

    forAll(fieldFlux.boundaryField(), patchI)
    {
        fieldFlux.boundaryField()[patchI] = 
            InternalContrib[patchI] - NeighbourContrib[patchI];
    }

    if (faceFluxCorrectionPtr_)
    {
        fieldFlux += *faceFluxCorrectionPtr_;
    }

    return tfieldFlux;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void fvMatrix<Type>::operator=(const fvMatrix<Type>& fvmv)
{
    if (this == &fvmv)
    {
        FatalErrorIn("fvMatrix<Type>::operator=(const fvMatrix<Type>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    if (&psi_ != &(fvmv.psi_))
    {
        FatalErrorIn("fvMatrix<Type>::operator=(const fvMatrix<Type>&)")
            << "different fields"
            << abort(FatalError);
    }

    lduMatrix::operator=(fvmv);
    source_ = fvmv.source_;
    internalCoeffs_ = fvmv.internalCoeffs_;
    boundaryCoeffs_ = fvmv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ = *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, fvPatchField, surfaceMesh>
        (*fvmv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void fvMatrix<Type>::operator=(const tmp<fvMatrix<Type> >& tfvmv)
{
    operator=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void fvMatrix<Type>::negate()
{
    lduMatrix::negate();
    source_.negate();
    internalCoeffs_.negate();
    boundaryCoeffs_.negate();

    if (faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_->negate();
    }
}


template<class Type>
void fvMatrix<Type>::operator+=(const fvMatrix<Type>& fvmv)
{
    checkMethod(*this, fvmv, "+=");

    dimensions_ += fvmv.dimensions_;
    lduMatrix::operator+=(fvmv);
    source_ += fvmv.source_;
    internalCoeffs_ += fvmv.internalCoeffs_;
    boundaryCoeffs_ += fvmv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ += *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, fvPatchField, surfaceMesh>
        (
            *fvmv.faceFluxCorrectionPtr_
        );
    }
}


template<class Type>
void fvMatrix<Type>::operator+=(const tmp<fvMatrix<Type> >& tfvmv)
{
    operator+=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void fvMatrix<Type>::operator-=(const fvMatrix<Type>& fvmv)
{
    checkMethod(*this, fvmv, "+=");

    dimensions_ -= fvmv.dimensions_;
    lduMatrix::operator-=(fvmv);
    source_ -= fvmv.source_;
    internalCoeffs_ -= fvmv.internalCoeffs_;
    boundaryCoeffs_ -= fvmv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && fvmv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ -= *fvmv.faceFluxCorrectionPtr_;
    }
    else if (fvmv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, fvPatchField, surfaceMesh>
        (-*fvmv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void fvMatrix<Type>::operator-=(const tmp<fvMatrix<Type> >& tfvmv)
{
    operator-=(tfvmv());
    tfvmv.clear();
}


template<class Type>
void fvMatrix<Type>::operator+=
(
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(*this, su, "+=");
    source() -= su.mesh().V()*su.internalField();
}


template<class Type>
void fvMatrix<Type>::operator+=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void fvMatrix<Type>::operator-=
(
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(*this, su, "-=");
    source() += su.mesh().V()*su.internalField();
}


template<class Type>
void fvMatrix<Type>::operator-=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void fvMatrix<Type>::operator+=
(
    const dimensioned<Type>& su
)
{
    source() -= psi().mesh().V()*su;
}


template<class Type>
void fvMatrix<Type>::operator-=
(
    const dimensioned<Type>& su
)
{
    source() += psi().mesh().V()*su;
}


template<class Type>
void fvMatrix<Type>::operator*=
(
    const volScalarField& vsf
)
{
    dimensions_ *= vsf.dimensions();
    lduMatrix::operator*=(vsf.internalField());
    source_ *= vsf.internalField();

    forAll(boundaryCoeffs_, patchI)
    {
        scalarField psf = vsf.boundaryField()[patchI].patchInternalField();
        internalCoeffs_[patchI] *= psf;
        boundaryCoeffs_[patchI] *= psf;
    }

    if (faceFluxCorrectionPtr_)
    {
        FatalErrorIn("fvMatrix<Type>::operator*=(const tmp<volScalarField>&)")
            << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}


template<class Type>
void fvMatrix<Type>::operator*=
(
    const tmp<volScalarField>& tvsf
)
{
    operator*=(tvsf());
    tvsf.clear();
}


template<class Type>
void fvMatrix<Type>::operator*=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ *= ds.dimensions();
    lduMatrix::operator*=(ds.value());
    source_ *= ds.value();
    internalCoeffs_ *= ds.value();
    boundaryCoeffs_ *= ds.value();

    if (faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ *= ds.value();
    }
}


template<class Type>
tmp<Field<Type> > fvMatrix<Type>::operator&
(
    const Field<Type>& psi
) const
{
    tmp<Field<Type> > tApsi(new Field<Type>(source_));
    Field<Type>& Apsi = tApsi();

    addBoundarySource(Apsi);

    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        scalarField psiCmpt = psi.component(cmpt);

        scalarField boundaryDiagCmpt(psi.size(), 0.0);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);

        scalarField ApsiCmpt(psi.size());

        // There is no boundaries on psi, so no interface coupling
        // 
        lduMatrix::Amul
        (
            ApsiCmpt,
            psiCmpt,
            FieldField<Field, scalar>(0),
            lduCoupledInterfacePtrsList(0),
            cmpt
        );

        Apsi.replace
        (
            cmpt, ApsiCmpt + boundaryDiagCmpt*psiCmpt - Apsi.component(cmpt)
        );
    }

    return tApsi;
}


template<class Type>
tmp<Field<Type> > fvMatrix<Type>::operator&
(
    const tmp<Field<Type> >& tpsi
) const
{
    tmp<Field<Type> > tHpsi = operator&(tpsi());
    tpsi.clear();
    return tHpsi;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const fvMatrix<Type>& fvm1,
    const fvMatrix<Type>& fvm2,
    const char* op
)
{
    if (&fvm1.psi() != &fvm2.psi())
    {
        FatalErrorIn
        (
            "checkMethod(const fvMatrix<Type>&, const fvMatrix<Type>&)"
        )   << "incompatible fields for operation "
            << endl << "    "
            << "[" << fvm1.psi().name() << "] "
            << op
            << " [" << fvm2.psi().name() << "]"
            << abort(FatalError);
    }

    if (dimensionSet::debug && fvm1.dimensions() != fvm2.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const fvMatrix<Type>&, const fvMatrix<Type>&)"
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm1.psi().name() << fvm1.dimensions()/dimVolume << " ] "
            << op
            << " [" << fvm2.psi().name() << fvm2.dimensions()/dimVolume << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void checkMethod
(
    const fvMatrix<Type>& fvm,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const char* op
)
{
    if (dimensionSet::debug && fvm.dimensions()/dimVolume != vf.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const fvMatrix<Type>&, const GeometricField<Type, "
            "fvPatchField, volMesh>&)"
        )   <<  "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm.psi().name() << fvm.dimensions()/dimVolume << " ] "
            << op
            << " [" << vf.name() << vf.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void checkMethod
(
    const fvMatrix<Type>& fvm,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if (dimensionSet::debug && fvm.dimensions()/dimVolume != dt.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const fvMatrix<Type>&, const dimensioned<Type>&)"
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fvm.psi().name() << fvm.dimensions()/dimVolume << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
lduMatrix::solverPerformance solve
(
    fvMatrix<Type>& fvm,
    Istream& solverControls
)
{
    return fvm.solve(solverControls);
}

template<class Type>
lduMatrix::solverPerformance solve
(
    const tmp<fvMatrix<Type> >& tfvm,
    Istream& solverControls
)
{
    lduMatrix::solverPerformance solverPerf = 
        ((fvMatrix<Type>&)tfvm()).solve(solverControls);

    tfvm.clear();
    return solverPerf;
}


template<class Type>
lduMatrix::solverPerformance solve(fvMatrix<Type>& fvm)
{
    return fvm.solve();
}

template<class Type>
lduMatrix::solverPerformance solve(const tmp<fvMatrix<Type> >& tfvm)
{
    lduMatrix::solverPerformance solverPerf =
        ((fvMatrix<Type>&)tfvm()).solve();

    tfvm.clear();
    return solverPerf;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() += B;
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() += B;
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<fvMatrix<Type> > tC(tB.ptr());
    tC() += A;
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() += tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() -= B;
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() -= B;
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<fvMatrix<Type> > tC(tB.ptr());
    tC() -= A;
    tC().negate();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() -= tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const fvMatrix<Type>& A,
    const fvMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}


template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const fvMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}


template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const fvMatrix<Type>& A,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}


template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<fvMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}


template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const fvMatrix<Type>& A,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const GeometricField<Type, fvPatchField, volMesh>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const GeometricField<Type, fvPatchField, volMesh>& su,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>& A,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const GeometricField<Type, fvPatchField, volMesh>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    tC().source() -= su.mesh().V()*su.internalField();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const GeometricField<Type, fvPatchField, volMesh>& su,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= su.value()*A.psi().mesh().V();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const tmp<fvMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const dimensioned<Type>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() -= su.value()*A.psi().mesh().V();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator+
(
    const dimensioned<Type>& su,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const tmp<fvMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const dimensioned<Type>& su,
    const fvMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().negate();
    tC().source() -= su.value()*A.psi().mesh().V();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator-
(
    const dimensioned<Type>& su,
    const tmp<fvMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC().source() -= su.value()*tC().psi().mesh().V();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const fvMatrix<Type>& A,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(A, su, "==");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += su.mesh().V()*su.internalField();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const fvMatrix<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tsu().mesh().V()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const fvMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC().source() += A.psi().mesh().V()*su.value();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator==
(
    const tmp<fvMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC().source() += tC().psi().mesh().V()*su.value();
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const volScalarField& vsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() *= vsf;
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const tmp<volScalarField>& tvsf,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() *= tvsf;
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const volScalarField& vsf,
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() *= vsf;
    return tC;
}

template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const tmp<volScalarField>& tvsf,
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() *= tvsf;
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const dimensioned<scalar>& ds,
    const fvMatrix<Type>& A
)
{
    tmp<fvMatrix<Type> > tC(new fvMatrix<Type>(A));
    tC() *= ds;
    return tC;
}


template<class Type>
tmp<fvMatrix<Type> > operator*
(
    const dimensioned<scalar>& ds,
    const tmp<fvMatrix<Type> >& tA
)
{
    tmp<fvMatrix<Type> > tC(tA.ptr());
    tC() *= ds;
    return tC;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<(Ostream& os, const fvMatrix<Type>& fvm)
{
    os  << (const lduMatrix&)(fvm) << nl
        << fvm.dimensions_ << nl
        << fvm.source_ << nl
        << fvm.internalCoeffs_ << nl
        << fvm.boundaryCoeffs_ << endl;

    os.check("Ostream& operator<<(Ostream&, fvMatrix<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

#include "fvMatrixSolve.C"

// ************************************************************************* //
