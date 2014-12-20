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
    lduMatrix is a general matrix class in which the coefficients are
    stored as three arrays, one for the upper triangle, one for the
    lower triangle and a third for the diagonal.  Addressing arrays must
    be supplied for the upper and lower triangles.  Diagonal, ICCG and BCG
    solvers are supplied and selected automatically by the solve member
    function according the the particular matrix structure.

    It would be better if this class were organised as a hierachy starting
    from an empty matrix, then deriving diagonal, symmetric and asymmetric
    matrices.  A job for the future.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "lduMatrix.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// declare static members

defineTypeNameAndDebug(lduMatrix, 1);

const scalar lduMatrix::great_ = 1.0e+20;
const scalar lduMatrix::small_ = 1.0e-20;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

lduMatrix::lduMatrix
(
    const lduAddressing& ldu,
    const patchScheduleList& lduCoupledInterfaceSchedule
)
:
    lduAddr_(ldu),
    lduCoupledInterfaceSchedule_(lduCoupledInterfaceSchedule),
    lowerPtr_(NULL),
    diagPtr_(NULL),
    upperPtr_(NULL),
    quietOperation_(false)
{}


lduMatrix::lduMatrix(const lduMatrix& A)
:
    lduAddr_(A.lduAddr_),
    lduCoupledInterfaceSchedule_(A.lduCoupledInterfaceSchedule_),
    lowerPtr_(NULL),
    diagPtr_(NULL),
    upperPtr_(NULL),
    quietOperation_(false)
{
    if (A.lowerPtr_)
    {
        lowerPtr_ = new scalarField(*(A.lowerPtr_));
    }

    if (A.diagPtr_)
    {
        diagPtr_ = new scalarField(*(A.diagPtr_));
    }

    if (A.upperPtr_)
    {
        upperPtr_ = new scalarField(*(A.upperPtr_));
    }
}


lduMatrix::lduMatrix
(
    const lduAddressing& ldu,
    const patchScheduleList& lduCoupledInterfaceSchedule,
    Istream& is
)
:
    lduAddr_(ldu),
    lduCoupledInterfaceSchedule_(lduCoupledInterfaceSchedule),
    lowerPtr_(new scalarField(is)),
    diagPtr_(new scalarField(is)),
    upperPtr_(new scalarField(is)),
    quietOperation_(false)
{}


lduMatrix::~lduMatrix()
{
    if (lowerPtr_)
    {
        delete lowerPtr_;
    }

    if (diagPtr_)
    {
        delete diagPtr_;
    }

    if (upperPtr_)
    {
        delete upperPtr_;
    }
}


scalarField& lduMatrix::lower()
{
    if (!lowerPtr_)
    {
        if (upperPtr_)
        {
            lowerPtr_ = new scalarField(*upperPtr_);
        }
        else
        {
            lowerPtr_ = new scalarField(lduAddr_.lowerAddr().size(), 0.0);
        }
    }

    return *lowerPtr_;
}


scalarField& lduMatrix::diag()
{
    if (!diagPtr_)
    {
        diagPtr_ = new scalarField(lduAddr_.size(), 0.0);
    }

    return *diagPtr_;
}


scalarField& lduMatrix::upper()
{
    if (!upperPtr_)
    {
        if (lowerPtr_)
        {
            upperPtr_ = new scalarField(*lowerPtr_);
        }
        else
        {
            upperPtr_ = new scalarField(lduAddr_.lowerAddr().size(), 0.0);
        }
    }

    return *upperPtr_;
}


const scalarField& lduMatrix::lower() const
{
    if (!lowerPtr_ && !upperPtr_)
    {
        FatalErrorIn("lduMatrix::lower() const")
            << "lowerPtr_ or upperPtr_ unallocated"
            << abort(FatalError);
    }

    if (lowerPtr_)
    {
        return *lowerPtr_;
    }
    else
    {
        return *upperPtr_;
    }
}


const scalarField& lduMatrix::diag() const
{
    if (!diagPtr_)
    {
        FatalErrorIn("const scalarField& lduMatrix::diag() const")
            << "diagPtr_ unallocated"
            << abort(FatalError);
    }

    return *diagPtr_;
}


const scalarField& lduMatrix::upper() const
{
    if (!lowerPtr_ && !upperPtr_)
    {
        FatalErrorIn("lduMatrix::upper() const")
            << "lowerPtr_ or upperPtr_ unallocated"
            << abort(FatalError);
    }

    if (upperPtr_)
    {
        return *upperPtr_;
    }
    else
    {
        return *lowerPtr_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const lduMatrix& ldum)
{
    if (ldum.lowerPtr_)
    {
        os  << "Lower triangle = "
            << *ldum.lowerPtr_
            << endl << endl;
    }

    if (ldum.diagPtr_)
    {
        os  << "diagonal = "
            << *ldum.diagPtr_
            << endl << endl;
    }

    if (ldum.upperPtr_)
    {
        os  << "Upper triangle = "
            << *ldum.upperPtr_
            << endl << endl;
    }

    os.check("Ostream& operator<<(Ostream&, const lduMatrix&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
