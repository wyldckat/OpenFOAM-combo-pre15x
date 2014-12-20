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

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "slicedVolFields.H"
#include "slicedSurfaceFields.H"
#include "SubField.H"
#include "demandDrivenData.H"
#include "lduAddressingFvMesh.H"
#include "emptyPolyPatch.H"
#include "mapPolyMesh.H"
#include "MapFvFields.H"
#include "fvMeshMapper.H"
#include "fvSurfaceMapper.H"
#include "mapClouds.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fvMesh, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fvMesh::clearGeomNotVol()
{
    deleteDemandDrivenData(SfPtr_);
    deleteDemandDrivenData(magSfPtr_);
    deleteDemandDrivenData(CPtr_);
    deleteDemandDrivenData(CfPtr_);
}

void fvMesh::clearGeom()
{
    clearGeomNotVol();
    deleteDemandDrivenData(V0Ptr_);
    deleteDemandDrivenData(V00Ptr_);

    // Mesh motion flux cannot be deleted here because the old-time flux
    // needs to be saved.
}

void fvMesh::clearAddressing()
{
    deleteDemandDrivenData(lduPtr_);
}

void fvMesh::clearOut()
{
    clearGeom();
    surfaceInterpolation::clearOut();

    clearAddressing();

    // Clear mesh motion flux
    deleteDemandDrivenData(phiPtr_);
}


Vector<label> fvMesh::makeDirections()
{
    Vector<label> dir;

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        dir[cmpt] = 1;
    }

    label nEmptyPatches = 0;

    vector dirVec = vector::zero;

    forAll(boundaryMesh(), patchi)
    {
        if
        (
            isA<emptyPolyPatch>(boundaryMesh()[patchi])
         && boundaryMesh()[patchi].size()
        )
        {
            nEmptyPatches++;
            dirVec += sum(cmptMag(boundaryMesh()[patchi].faceAreas()));
        }
    }

    if (nEmptyPatches)
    {
        reduce(dirVec, sumOp<vector>());

        dirVec /= mag(dirVec);

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (dirVec[cmpt] > SMALL)
            {
                dir[cmpt] = -1;
            }
            else
            {
                dir[cmpt] = 1;
            }
        }
    }

    return dir;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fvMesh::fvMesh(const IOobject& io)
:
    polyMesh(io),
    surfaceInterpolation(*this),
    boundary_(*this, boundaryMesh()),
    directions_(makeDirections()),
    lduPtr_(NULL),
    curTimeIndex_(time().timeIndex()),
    V0Ptr_(NULL),
    V00Ptr_(NULL),
    SfPtr_(NULL),
    magSfPtr_(NULL),
    CPtr_(NULL),
    CfPtr_(NULL),
    phiPtr_(NULL)
{
    if (debug)
    {
        Info<< "Constructing fvMesh from IOobject"
            << endl;
    }

    if (file(time().timePath()/"V0"))
    {
        V0Ptr_ = new scalarIOField
        (
            IOobject
            (
                "V0",
                time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        );
    }
}


// Construct from components
fvMesh::fvMesh
(
    const IOobject& io,
    const pointField& points,
    const faceList& faces,
    const cellList& cells
)
:
    polyMesh(io, points, faces, cells),
    surfaceInterpolation(*this),
    boundary_(*this),
    directions_(makeDirections()),
    lduPtr_(NULL),
    curTimeIndex_(time().timeIndex()),
    V0Ptr_(NULL),
    V00Ptr_(NULL),
    SfPtr_(NULL),
    magSfPtr_(NULL),
    CPtr_(NULL),
    CfPtr_(NULL),
    phiPtr_(NULL)
{
    if (debug)
    {
        Info<< "Constructing fvMesh from components"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fvMesh::~fvMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Helper function for construction from pieces
void fvMesh::addFvPatches(const List<polyPatch*> & p)
{
    if (boundary().size() > 0)
    {
        FatalErrorIn("fvMesh::addFvPatches(const List<polyPatch*>& p)")
            << " boundary already exists"
            << abort(FatalError);
    }

    // first add polyPatches
    addPatches(p);
    boundary_.addPatches(boundaryMesh());
}


void fvMesh::removeFvBoundary()
{
    if (debug)
    {
        Info<< "void fvMesh::removeFvBoundary(): "
            << "Removing boundary patches."
            << endl;
    }

    // Remove fvBoundaryMesh data first.
    boundary_.clear();
    boundary_.setSize(0);
    polyMesh::removeBoundary();

    clearOut();
}


polyMesh::readUpdateState fvMesh::readUpdate()
{
    if (debug)
    {
        Info<< "polyMesh::readUpdateState fvMesh::readUpdate() : "
            << "Updating fvMesh.  ";
    }

    polyMesh::readUpdateState state = polyMesh::readUpdate();

    if (state == polyMesh::TOPO_PATCH_CHANGE)
    {
        if (debug)
        {
            Info << "Boundary and topological update" << endl;
        }

        boundary_.readUpdate(boundaryMesh());

        clearOut();
        
    }
    else if (state == polyMesh::TOPO_CHANGE)
    {
        if (debug)
        {
            Info << "Topological update" << endl;
        }

        clearOut();
    }
    else if (state == polyMesh::POINTS_MOVED)
    {
        if (debug)
        {
            Info << "Point motion update" << endl;
        }

        clearGeom();
    }
    else
    {
        if (debug)
        {
            Info << "No update" << endl;
        }
    }

    return state;
}


const fvBoundaryMesh& fvMesh::boundary() const
{
    return boundary_;
}


const lduAddressing& fvMesh::ldu() const
{
    if (!lduPtr_)
    {
        lduPtr_ = new lduAddressingFvMesh(*this);
    }

    return *lduPtr_;
}


void fvMesh::constructAndClear() const
{
    // Construct all necessary data and clear all storage

    static bool done = false;

    if (!done)
    {
        if (debug)
        {
            Info<< "void fvMesh::constructAndClear()"
                << endl;
        }

        // Force creation of geometric data
        Cf();
        Sf();
        magSf();

        C();
        V();

        // Force creation of interpolation weights
        weights();
        deltaCoeffs();

        // Clear ldu addressing - it holds references to deleted data.
        // It will be created on demand
        deleteDemandDrivenData(lduPtr_);

        //const_cast<fvMesh&>(*this).clearPrimitives();
        //const_cast<fvMesh&>(*this).polyMesh::clearGeom();

        done = true;
    }
}


void fvMesh::mapFields()
{
    const mapPolyMesh& meshMap = morphMap();

    // Create a mapper
    const fvMeshMapper mapper(*this);

    // Map all the volFields in the objectRegistry
    MapGeometricFields<scalar, fvPatchField, fvMeshMapper, volMesh>(mapper);
    MapGeometricFields<vector, fvPatchField, fvMeshMapper, volMesh>(mapper);
    MapGeometricFields<tensor, fvPatchField, fvMeshMapper, volMesh>(mapper);

    // Map all the surfaceFields in the objectRegistry
    MapGeometricFields<scalar, fvPatchField, fvMeshMapper, surfaceMesh>(mapper);
    MapGeometricFields<vector, fvPatchField, fvMeshMapper, surfaceMesh>(mapper);
    MapGeometricFields<tensor, fvPatchField, fvMeshMapper, surfaceMesh>(mapper);

    // Map all the clouds in the objectRegistry
    mapClouds(*this, meshMap);


    if (V0Ptr_)
    {
        scalarField& V0 = *V0Ptr_;

        scalarField V00(V0);
        V0.setSize(nCells());

        const labelList& cellMap = meshMap.cellMap();

        forAll(V0, i)
        {
            if (cellMap[i] > -1)
            {
                V0[i] = V00[cellMap[i]];
            }
            else
            {
                V0[i] = 0.0;
            }
        }
    }

    if (V00Ptr_)
    {
        MapInternalField<scalar, fvMeshMapper, volMesh>()(*V00Ptr_, mapper);
    }
}


void fvMesh::handleMorph()
{
    surfaceInterpolation::clearOut();
    clearGeomNotVol();

    // Map all fields
    mapFields();

    clearAddressing();

    // updateTopology() should also clear out the surfaceInterpolation.
    // This is a temporary solution
    surfaceInterpolation::movePoints();
}


void fvMesh::updateTopology()
{
    if (morphing())
    {
        // Have polyMesh do all topo changes
        polyMesh::updateTopology();

        handleMorph();
    }
}


tmp<scalarField> fvMesh::movePoints(const pointField& p)
{
    // Grab old time volumes if the time has been incremented
    if (curTimeIndex_ < time().timeIndex())
    {
        if (V00Ptr_ && V0Ptr_)
        {
            *V00Ptr_ = *V0Ptr_;
        }

        if (V0Ptr_)
        {
            *V0Ptr_ = V();
        }
        else
        {
            V0Ptr_ = new scalarIOField
            (
                IOobject
                (
                    "V0",
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                V()
            );
        }

        curTimeIndex_ = time().timeIndex();
    }


    // delete out of date geometrical information
    clearGeomNotVol();


    if (!phiPtr_)
    {
        makePhi();
    }

    surfaceScalarField& phi = *phiPtr_;

    // Grab old time mesh motion fluxes if the time has been incremented
    if (phi.timeIndex() != time().timeIndex())
    {
        phi.oldTime();
    }


    // move the polyMesh and grab mesh motion fluxes

    scalar rDeltaT = 1.0/time().deltaT().value();

    tmp<scalarField> tsweptVols = polyMesh::movePoints(p);
    scalarField& sweptVols = tsweptVols();

    phi.internalField() = scalarField::subField(sweptVols, nInternalFaces());
    phi.internalField() *= rDeltaT;

    const fvPatchList& patches = boundary();

    forAll (patches, patchI)
    {
        phi.boundaryField()[patchI] = patches[patchI].patchSlice(sweptVols);
        phi.boundaryField()[patchI] *= rDeltaT;
    }

    boundary_.movePoints();
    surfaceInterpolation::movePoints();

    return tsweptVols;
}


//- Write the underlying polyMesh and other data
bool fvMesh::write
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    return polyMesh::write(fmt, ver, cmp);
}


//- Write mesh using IO settings from the time
bool fvMesh::write() const
{
    return polyMesh::write();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool fvMesh::operator!=(const fvMesh& bm) const
{
    return &bm != this;
}


bool fvMesh::operator==(const fvMesh& bm) const
{
    return &bm == this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
