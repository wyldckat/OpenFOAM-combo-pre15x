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
    Tetrahedral Finite Element matrix member functions and operators

\*---------------------------------------------------------------------------*/

#include "PstreamReduceOps.H"

#include "tetFemMatrix.H"
#include "tetPointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const label tetFemMatrix<Type>::fixFillIn = 2;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
tetFemMatrix<Type>::tetFemMatrix
(
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& psi,
    const dimensionSet& ds
)
:
    lduMatrix(psi.mesh().ldu(), psi.mesh().parallelData().patchSchedule()),
    psi_(psi),
    dimensions_(ds),
    source_(psi.size(), pTraits<Type>::zero),
    boundaryConditionsSet_(false),
    fixedEqns_(psi.mesh().ldu().size()/fixFillIn),
    solvingComponent(0)
{
#   ifdef DEBUGtetFemMatrix
    if (debug)
    {
        Info<< "tetFemMatrix<Type>(GeometricField<Type, tetPolyPatchField, "
            << "tetPointMesh>&, const dimensionSet&) : "
            << "constructing tetFemMatrix<Type> for field " << psi_.name()
            << endl;
    }
#   endif
}


template<class Type>
tetFemMatrix<Type>::tetFemMatrix(const tetFemMatrix<Type>& tetFem)
:
    refCount(),
    lduMatrix(tetFem),
    psi_(tetFem.psi_),
    dimensions_(tetFem.dimensions_),
    source_(tetFem.source_),
    boundaryConditionsSet_(false),
    fixedEqns_(psi_.mesh().ldu().size()/fixFillIn),
    solvingComponent(0)
{
#   ifdef DEBUGtetFemMatrix
    if (debug)
    {
        Info<< "tetFemMatrix<Type>::tetFemMatrix(const tetFemMatrix<Type>&) : "
            << "copying tetFemMatrix<Type> for field " << psi_.name()
            << endl;
    }
#   endif
}


template<class Type>
tetFemMatrix<Type>::tetFemMatrix
(
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& psi,
    Istream& is
)
:
    lduMatrix(psi.mesh().ldu(), psi.mesh().parallelData().patchSchedule()),
    psi_(psi),
    dimensions_(is),
    source_(is),
    boundaryConditionsSet_(false),
    fixedEqns_(psi.mesh().ldu().size()/fixFillIn),
    solvingComponent(0)
{
#   ifdef DEBUGtetFemMatrix
    if (debug)
    {
        Info<< "tetFemMatrix<Type>(GeometricField<Type, tetPolyPatchField, "
            << "tetPointMesh>&, Istream&) : "
            << "constructing tetFemMatrix<Type> for field " << psi_.name()
            << endl;
    }
#   endif
}


template<class Type>
tetFemMatrix<Type>::~tetFemMatrix()
{
#   ifdef DEBUGtetFemMatrix
    if (debug)
    {
        Info<< "tetFemMatrix<Type>::~tetFemMatrix<Type>() : "
            << "destroying tetFemMatrix<Type> for field " << psi_.name()
            << endl;
    }
#   endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Does the matrix need a reference level for solution
// template<class Type>
// bool tetFemMatrix<Type>::needReference()
// {
//     // Search all boundary conditions, if any are
//     // fixed-value or mixed (Robin) do not set reference level for solution.

//     static bool searched = false;
//     static bool needRef = true;

//     if (!searched)
//     {
//         const BoundaryField<Type>& patchFields = psi_.boundaryField();

//         forAll(patchFields, patchi)
//         {
//             if (patchFields[patchi].fixesValue())
//             {
//                 needRef = false;
//             }
//         }

//         reduce(needRef, andOp<bool>());

//         searched = true;
//     }

//     return needRef;
// }


// Set reference level for solution
template<class Type>
void tetFemMatrix<Type>::addConstraint
(
    const label vertex,
    const Type& value
)
{
    constraint<Type> cp(vertex, value);

    if (!fixedEqns_.found(vertex))
    {
        fixedEqns_.insert(vertex, cp);
    }
    else
    {
        Warning
            << "void tetFemMatrix<Type>::addConstraint(const label vertex, "
            << "const Type& value) : " << nl
            << "Adding onstraint on an already constrained point."
            << "  Point: " << vertex
            << endl;

        fixedEqns_[vertex].combine(cp);
    }
}


template<class Type>
void tetFemMatrix<Type>::relax(const scalar alpha)
{
    scalarField oldDiag(diag());

    lduMatrix::relax(alpha);

    source() += (diag() - oldDiag)*psi_.internalField();
}


template<class Type>
void tetFemMatrix<Type>::relax()
{
    scalar alpha = 0;

    if (psi_.mesh().relax(psi_.name()))
    {
        alpha = psi_.mesh().relaxationFactor(psi_.name());
    }

    if (alpha > 0)
    {
        relax(alpha);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void tetFemMatrix<Type>::operator=(const tetFemMatrix<Type>& tetFem)
{
    if (this == &tetFem)
    {
        FatalErrorIn
        (
            "tetFemMatrix<Type>::operator=(const tetFemMatrix<Type>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    if (&psi_ != &(tetFem.psi_))
    {
        FatalErrorIn
        (
            "tetFemMatrix<Type>::operator=(const tetFemMatrix<Type>&)"
        )   << "different fields"
            << abort(FatalError);
    }

    lduMatrix::operator=(tetFem);
    source_ = tetFem.source_;
    boundaryConditionsSet_ = false;
    fixedEqns_.clear();
}


template<class Type>
void tetFemMatrix<Type>::operator=(const tmp<tetFemMatrix<Type> >& ttetFem)
{
    operator=(ttetFem());
    ttetFem.clear();
}


template<class Type>
void tetFemMatrix<Type>::negate()
{
    lduMatrix::negate();
    source_.negate();
}


template<class Type>
void tetFemMatrix<Type>::operator+=(const tetFemMatrix<Type>& tetFem)
{
    checkMethod(*this, tetFem, "+=");

    dimensions_ += tetFem.dimensions_;
    lduMatrix::operator+=(tetFem);
    source_ += tetFem.source_;
}


template<class Type>
void tetFemMatrix<Type>::operator+=(const tmp<tetFemMatrix<Type> >& ttetFem)
{
    operator+=(ttetFem());
    ttetFem.clear();
}


template<class Type>
void tetFemMatrix<Type>::operator+=
(
    const GeometricField<Type, tetPolyPatchField, elementMesh>& su
)
{
    checkMethod(*this, su, "+=");
    source() -= distributeField(su.internalField());
}


template<class Type>
void tetFemMatrix<Type>::operator+=
(
    const tmp<GeometricField<Type, tetPolyPatchField, elementMesh> >& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void tetFemMatrix<Type>::operator-=(const tetFemMatrix<Type>& tetFem)
{
    checkMethod(*this, tetFem, "+=");

    dimensions_ -= tetFem.dimensions_;
    lduMatrix::operator-=(tetFem);
    source_ -= tetFem.source_;
}


template<class Type>
void tetFemMatrix<Type>::operator-=(const tmp<tetFemMatrix<Type> >& ttetFem)
{
    operator-=(ttetFem());
    ttetFem.clear();
}


template<class Type>
void tetFemMatrix<Type>::operator-=
(
    const GeometricField<Type, tetPolyPatchField, elementMesh>& su
)
{
    checkMethod(*this, su, "-=");
    source() += distributeField(su.internalField());
}


template<class Type>
void tetFemMatrix<Type>::operator-=
(
    const tmp<GeometricField<Type, tetPolyPatchField, elementMesh> >& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void tetFemMatrix<Type>::operator+=
(
    const dimensioned<Type>& su
)
{
    checkMethod(*this, su, "+=");
    source() -= distributeField(Field<Type>(psi.mesh().nCells(), su.value()));
}


template<class Type>
void tetFemMatrix<Type>::operator-=
(
    const dimensioned<Type>& su
)
{
    checkMethod(*this, su, "-=");
    source() += distributeField(Field<Type>(psi.mesh().nCells(), su.value()));
}


template<class Type>
void tetFemMatrix<Type>::operator*=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ *= ds.dimensions();
    lduMatrix::operator*=(ds.value());
    source_ *= ds.value();
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const tetFemMatrix<Type>& tetFem1,
    const tetFemMatrix<Type>& tetFem2,
    const char* op
)
{
    if (&tetFem1.psi() != &tetFem2.psi())
    {
        FatalErrorIn
        (
            "checkMethod(const tetFemMatrix<Type>&, "
            "const tetFemMatrix<Type>&) : "
        )   << "incompatible fields for operation "
            << endl << "    "
            << "[" << tetFem1.psi().name() << "] "
            << op
            << " [" << tetFem1.psi().name() << "]"
            << abort(FatalError);
    }

    if (dimensionSet::debug && tetFem1.dimensions() != tetFem2.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const tetFemMatrix<Type>&, "
            "const tetFemMatrix<Type>&) : "
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << tetFem1.psi().name() << tetFem1.dimensions()/dimVolume
            << " ] "
            << op
            << " [" << tetFem1.psi().name() << tetFem2.dimensions()/dimVolume
            << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void checkMethod
(
    const tetFemMatrix<Type>& tetFem,
    const GeometricField<Type, tetPolyPatchField, elementMesh>& vf,
    const char* op
)
{
    if
    (
        dimensionSet::debug
     && tetFem.dimensions()/dimVolume != vf.dimensions()
    )
    {
        FatalErrorIn
        (
            "checkMethod(const tetFemMatrix<Type>&, "
            "const GeometricField<Type, tetPolyPatchField, "
            "elementMesh>&) : "
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << tetFem.psi().name() << tetFem.dimensions()/dimVolume
            << " ] "
            << op
            << " [" << tetFem.psi().name() << vf.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void checkMethod
(
    const tetFemMatrix<Type>& tetFem,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if
    (
        dimensionSet::debug
     && tetFem.dimensions()/dimVolume != dt.dimensions()
    )
    {
        FatalErrorIn
        (
            "checkMethod(const tetFemMatrix<Type>&, "
            "const dimensioned<Type>&) : "
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << tetFem.psi().name() << tetFem.dimensions()/dimVolume
            << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
lduMatrix::solverPerformance solve(tetFemMatrix<Type>& tetFem)
{
    return tetFem.solve();
}


template<class Type>
lduMatrix::solverPerformance solve(const tmp<tetFemMatrix<Type> >& ttetFem)
{
    lduMatrix::solverPerformance solverPerf =
        ((tetFemMatrix<Type>&)ttetFem()).solve();

    ttetFem.clear();
    return solverPerf;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>& A,
    const tetFemMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += B;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tetFemMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += B;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>& A,
    const tmp<tetFemMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<tetFemMatrix<Type> > tC(tB.ptr());
    tC() += A;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tmp<tetFemMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>& A
)
{
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC().negate();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >& tA
)
{
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC().negate();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>& A,
    const tetFemMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= B;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tetFemMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= B;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>& A,
    const tmp<tetFemMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<tetFemMatrix<Type> > tC(tB.ptr());
    tC() -= A;
    tC().negate();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tmp<tetFemMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>& A,
    const tetFemMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tetFemMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>& A,
    const tmp<tetFemMatrix<Type> >& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tmp<tetFemMatrix<Type> >& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>& A,
    const GeometricField<Type, tetPolyPatchField, elementMesh>& su
)
{
    checkMethod(A, su, "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >& tA,
    const GeometricField<Type, tetPolyPatchField, elementMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>& A,
    const tmp<GeometricField<Type, tetPolyPatchField, elementMesh> >& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= tsu();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tmp<GeometricField<Type, tetPolyPatchField, elementMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= tsu();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const GeometricField<Type, tetPolyPatchField, elementMesh>& su,
    const tetFemMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const GeometricField<Type, tetPolyPatchField, elementMesh>& su,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, tetPolyPatchField, elementMesh> >& tsu,
    const tetFemMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= tsu();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<GeometricField<Type, tetPolyPatchField, elementMesh> >& tsu,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= tsu();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>& A,
    const GeometricField<Type, tetPolyPatchField, elementMesh>& su
)
{
    checkMethod(A, su, "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >& tA,
    const GeometricField<Type, tetPolyPatchField, elementMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>& A,
    const tmp<GeometricField<Type, tetPolyPatchField, elementMesh> >& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += tsu();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tmp<GeometricField<Type, tetPolyPatchField, elementMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += tsu();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const GeometricField<Type, tetPolyPatchField, elementMesh>& su,
    const tetFemMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC().negate();
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const GeometricField<Type, tetPolyPatchField, elementMesh>& su,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC() -= su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, tetPolyPatchField, elementMesh> >& tsu,
    const tetFemMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC().negate();
    tC() -= tsu();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<GeometricField<Type, tetPolyPatchField, elementMesh> >& tsu,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC() -= tsu();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tetFemMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const tmp<tetFemMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const dimensioned<Type>& su,
    const tetFemMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator+
(
    const dimensioned<Type>& su,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tetFemMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const tmp<tetFemMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const dimensioned<Type>& su,
    const tetFemMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC().negate();
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator-
(
    const dimensioned<Type>& su,
    const tmp<tetFemMatrix<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC().negate();
    tC() -= su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>& A,
    const GeometricField<Type, tetPolyPatchField, elementMesh>& su
)
{
    checkMethod(A, su, "==");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >& tA,
    const GeometricField<Type, tetPolyPatchField, elementMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += su;
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>& A,
    const tmp<GeometricField<Type, tetPolyPatchField, elementMesh> >& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += tsu();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >& tA,
    const tmp<GeometricField<Type, tetPolyPatchField, elementMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += tsu();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tetFemMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() += su;
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator==
(
    const tmp<tetFemMatrix<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<tetFemMatrix<Type> > tC(tA.ptr());
    tC() += su.value();
    return tC;
}


template<class Type>
tmp<tetFemMatrix<Type> > operator*
(
    const dimensioned<scalar>& ds,
    const tetFemMatrix<Type>& A
)
{
    tmp<tetFemMatrix<Type> > tC(new tetFemMatrix<Type>(A));
    tC() *= ds;
    return tC;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<(Ostream& os, const tetFemMatrix<Type>& tetFem)
{
    os  << (const lduMatrix&)(tetFem) << nl
        << tetFem.dimensions_ << nl
        << tetFem.source_ << endl;

    os.check("Ostream& operator<<(Ostream&, tetFemMatrix<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
