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

#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fvMesh::makeSf() const
{
    if (debug)
    {
        Info<< "void fvMesh::makeSf() : "
            << "assembling face areas"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (SfPtr_)
    {
        FatalErrorIn("fvMesh::makeSf()")
            << "face areas already exist"
            << abort(FatalError);
    }

    SfPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "S",
            pointsInstance(),
            meshSubDir,
            *this
        ),
        *this,
        dimArea
    );

    surfaceVectorField& s = *SfPtr_;

    const vectorField& allFaceAreas = faceAreas();

    s.internalField() =
        vectorField::subField(allFaceAreas, nInternalFaces());

    const fvPatchList& patches = boundary();

    forAll (patches, patchI)
    {
        s.boundaryField()[patchI] = patches[patchI].patchSlice(allFaceAreas);
    }
}


void fvMesh::makeMagSf() const
{
    if (debug)
    {
        Info<< "void fvMesh::makeMagSf() : "
            << "assembling mag face areas"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (magSfPtr_)
    {
        FatalErrorIn("void fvMesh::makeMagSf()")
            << "mag face areas already exist"
            << abort(FatalError);
    }

    // Note: Added stabilisation for faces with exactly zero area.
    // These should be caught on mesh checking but at least this stops
    // the code from producing Nans.  
    magSfPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "magSf",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mag(Sf())
#       ifdef BAD_MESH_STABILISATION
      + dimensionedScalar("vs", dimArea, VSMALL)
#       endif
    );
}


void fvMesh::makeC() const
{
    if (debug)
    {
        Info<< "void fvMesh::makeC() : "
            << "assembling cell centres"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (CPtr_)
    {
        FatalErrorIn("fvMesh::makeC()")
            << "cell centres already exist"
            << abort(FatalError);
    }

    CPtr_ = new volVectorField
    (
        IOobject
        (
            "C",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        *this,
        dimLength
    );

    volVectorField& c = *CPtr_;

    const vectorField& allCellCentres = cellCentres();
    const vectorField& fcs = primitiveMesh::faceCentres();

    c.internalField() =
        vectorField::subField(allCellCentres, nCells());

    const fvPatchList& patches = boundary();

    forAll (patches, patchI)
    {
        c.boundaryField()[patchI] = patches[patchI].patchSlice(fcs);
    }
}


void fvMesh::makeCf() const
{
    if (debug)
    {
        Info<< "void fvMesh::makeCf() : "
            << "assembling face centres"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (CfPtr_)
    {
        FatalErrorIn("fvMesh::makeCf()")
            << "face centres already exist"
            << abort(FatalError);
    }

    CfPtr_ = new surfaceVectorField
    (
        IOobject
        (
            "Cf",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        *this,
        dimLength
    );

    surfaceVectorField& fc = *CfPtr_;

    const vectorField& fcs = primitiveMesh::faceCentres();

    fc.internalField() =
        vectorField::subField(fcs, nInternalFaces());

    const fvPatchList& patches = boundary();

    forAll (patches, patchI)
    {
        fc.boundaryField()[patchI] = patches[patchI].patchSlice(fcs);
    }
}


void fvMesh::makePhi() const
{
    if (debug)
    {
        Info<< "void fvMesh::makePhi() : "
            << "creating zero flux field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (phiPtr_)
    {
        FatalErrorIn("fvMesh::makePhi()")
            << "flux field already exists"
            << abort(FatalError);
    }

    phiPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "meshPhi",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        *this,
        dimensionedScalar("0", dimVolume/dimTime, 0.0)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const scalarField& fvMesh::V() const
{
    return cellVolumes();
}


const scalarField& fvMesh::V0() const
{
    if (!V0Ptr_)
    {
        FatalErrorIn("fvMesh::V0() const")
            << "V0 is not available"
            << abort(FatalError);
    }

    return *V0Ptr_;
}


scalarField& fvMesh::setV0()
{
    if (!V0Ptr_)
    {
        FatalErrorIn("fvMesh::setV0()")
            << "V0 is not available"
            << abort(FatalError);
    }

    return *V0Ptr_;
}


const scalarField& fvMesh::V00() const
{
    if (!V00Ptr_)
    {
        V00Ptr_ = new scalarIOField
        (
            IOobject
            (
                "V00",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            V0()
        );

        V0Ptr_->writeOpt() = IOobject::AUTO_WRITE;
    }

    return *V00Ptr_;
}


const surfaceVectorField& fvMesh::Sf() const
{
    if (!SfPtr_)
    {
        makeSf();
    }

    return *SfPtr_;
}


const surfaceScalarField& fvMesh::magSf() const
{
    if (!magSfPtr_)
    {
        makeMagSf();
    }

    return *magSfPtr_;
}


const volVectorField& fvMesh::C() const
{
    if (!CPtr_)
    {
        makeC();
    }

    return *CPtr_;
}


const surfaceVectorField& fvMesh::Cf() const
{
    if (!CfPtr_)
    {
        makeCf();
    }

    return *CfPtr_;
}


const surfaceScalarField& fvMesh::phi() const
{
    if (!phiPtr_)
    {
        makePhi();
    }

    return *phiPtr_;
}


surfaceScalarField& fvMesh::setPhi()
{
    if (!phiPtr_)
    {
        makePhi();
    }

    return *phiPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
