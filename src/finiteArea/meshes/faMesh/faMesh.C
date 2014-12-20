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

Description

\*---------------------------------------------------------------------------*/

#include "faMesh.H"
#include "Time.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "IndirectList.H"
#include "lduAddressingFaMesh.H"
#include "areaFields.H"
#include "edgeFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(faMesh, 0);

word faMesh::meshSubDir = "faMesh";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void faMesh::setPrimitiveMeshData()
{
    if (debug)
    {
        Info<< "void faMesh::setPrimitiveMeshData() const : "
            << "Setting primitive data" << endl;
    }

    const indirectPrimitivePatch& bp = patch();

    // Set faMesh edges
    edges_.setSize(bp.nEdges());

    label edgeI = -1;


    label nIntEdges = bp.nInternalEdges();

    for (label curEdge = 0; curEdge < nIntEdges; curEdge++)
    {
        edges_[++edgeI] =
            bp.edges()[curEdge];
    }

    forAll (boundary(), patchI)
    {
        forAll (boundary()[patchI], eI_)
        {
            edges_[++edgeI] =
                bp.edges()[boundary()[patchI][eI_]];
        }
    }

    nEdges_ = edges_.size();
    nInternalEdges_ = nIntEdges;


    // Set edge owner and neighbour
    edgeOwner_.setSize(nEdges());
    edgeNeighbour_.setSize(nInternalEdges());

    edgeI = -1;

    for (label curEdge = 0; curEdge < nIntEdges; curEdge++)
    {
        edgeOwner_[++edgeI] =
            bp.edgeFaces()[curEdge][0];

        edgeNeighbour_[edgeI] =
            bp.edgeFaces()[curEdge][1];
    }

    forAll (boundary(), patchI)
    {
        forAll (boundary()[patchI], eI)
        {
            edgeOwner_[++edgeI] =
                bp.edgeFaces()[boundary()[patchI][eI]][0];
        }
    }

    // Set number of faces
    nFaces_ = bp.size();

    // Set number of points
    nPoints_ = bp.nPoints();
}


void faMesh::clearGeomNotAreas() const
{
    if (debug)
    {
        Info<< "void faMesh::clearGeomNotAreas() const : "
            << "Clearing geometry" << endl;
    }

    deleteDemandDrivenData(patchPtr_);
    deleteDemandDrivenData(patchStartsPtr_);
    deleteDemandDrivenData(LePtr_);
    deleteDemandDrivenData(magLePtr_);
    deleteDemandDrivenData(centresPtr_);
    deleteDemandDrivenData(edgeCentresPtr_);
    deleteDemandDrivenData(faceAreaNormalsPtr_);
    deleteDemandDrivenData(edgeAreaNormalsPtr_);
    deleteDemandDrivenData(pointAreaNormalsPtr_);
    deleteDemandDrivenData(faceCurvaturesPtr_);
    deleteDemandDrivenData(edgeTransformTensorsPtr_);
    deleteDemandDrivenData(parallelDataPtr_);
    deleteDemandDrivenData(SPtr_);
}


void faMesh::clearGeom() const
{
    if (debug)
    {
        Info<< "void faMesh::clearGeom() const : "
            << "Clearing geometry" << endl;
    }

    clearGeomNotAreas();
    deleteDemandDrivenData(S0Ptr_);
    deleteDemandDrivenData(S00Ptr_);
}


void faMesh::clearAddressing() const
{
    if (debug)
    {
        Info<< "void faMesh::clearAddressing() const : "
            << "Clearing addressing" << endl;
    }

    deleteDemandDrivenData(lduPtr_);
}


void faMesh::clearOut() const
{
    clearGeom();
    clearAddressing();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from objectRegistry
faMesh::faMesh
(
    const polyMesh& m,
    IOobject::readOption r,
    IOobject::writeOption w
)
:
    edgeInterpolation(*this, m),
    mesh_(m),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            m.time().findInstance(meshDir(), "faceLabels"),
            meshSubDir,
            m,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    boundary_
    (
        IOobject
        (
            "boundary",
            m.time().findInstance(meshDir(), "boundary"),
            meshSubDir,
            m,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    patchPtr_(NULL),
    lduPtr_(NULL),
    SPtr_(NULL),
    S0Ptr_(NULL),
    S00Ptr_(NULL),
    patchStartsPtr_(NULL),
    LePtr_(NULL),
    magLePtr_(NULL),
    centresPtr_(NULL),
    edgeCentresPtr_(NULL),
    faceAreaNormalsPtr_(NULL),
    edgeAreaNormalsPtr_(NULL),
    pointAreaNormalsPtr_(NULL),
    faceCurvaturesPtr_(NULL),
    edgeTransformTensorsPtr_(NULL),
    parallelDataPtr_(NULL),
    moving_(false),
    curMotionTimeIndex_(m.time().timeIndex())
{
    if (debug)
    {
        Info<< "faMesh::faMesh(...) : "
            << "Creating faMesh from objectRegistry" << endl;
    }

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    setPrimitiveMeshData();

    if (file(m.time().timePath()/"S0"))
    {
        S0Ptr_ = new scalarIOField
        (
            IOobject
            (
                "S0",
                m.time().timeName(),
                meshSubDir,
                m.db(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        );
    }
}


// Construct from components without boundary.
faMesh::faMesh
(
    const polyMesh& m,
    const labelList& faceLabels,
    IOobject::readOption r,
    IOobject::writeOption w
)
:
    edgeInterpolation(*this, m),
    mesh_(m),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            m.instance(),
            meshSubDir,
            m,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        faceLabels
    ),
    boundary_
    (
        IOobject
        (
            "boundary",
            m.instance(),
            meshSubDir,
            m,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    patchPtr_(NULL),
    lduPtr_(NULL),
    SPtr_(NULL),
    S0Ptr_(NULL),
    S00Ptr_(NULL),
    patchStartsPtr_(NULL),
    LePtr_(NULL),
    magLePtr_(NULL),
    centresPtr_(NULL),
    edgeCentresPtr_(NULL),
    faceAreaNormalsPtr_(NULL),
    edgeAreaNormalsPtr_(NULL),
    pointAreaNormalsPtr_(NULL),
    faceCurvaturesPtr_(NULL),
    edgeTransformTensorsPtr_(NULL),
    parallelDataPtr_(NULL),
    moving_(false),
    curMotionTimeIndex_(m.time().timeIndex())
{
    if (debug)
    {
        Info<< "faMesh::faMesh(...) : "
            << "Creating faMesh from components" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

faMesh::~faMesh()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

fileName faMesh::meshDir() const
{
    return mesh_.dbDir()/meshSubDir;
}


const Time& faMesh::time() const
{
    return mesh_.time();
}


const indirectPrimitivePatch& faMesh::patch() const
{
    if (!patchPtr_)
    {
        patchPtr_ = new indirectPrimitivePatch
        (
            IndirectList<face>
            (
                mesh_.faces(),
                faceLabels_
            ),
            mesh_.points()
        );
    }
        
    return *patchPtr_;
}


indirectPrimitivePatch& faMesh::patch()
{
    if (!patchPtr_)
    {
        patchPtr_ = new indirectPrimitivePatch
        (
            IndirectList<face>
            (
                mesh_.faces(), 
                faceLabels_
            ),
            mesh_.points()
        );
    }
        
    return *patchPtr_;
}


const pointField& faMesh::points() const
{
    return patch().localPoints();
}


const edgeList& faMesh::edges() const
{
    return edges_;
}


const faceList& faMesh::faces() const
{
    return patch().localFaces();
}


void faMesh::addFaPatches(const List<faPatch*>& p)
{
    if (debug)
    {
        Info<< "void faMesh::addFaPatches(const List<faPatch*>& p) : "
            << "Adding patches to faMesh" << endl;
    }

    if (boundary().size() > 0)
    {
        FatalErrorIn("void faMesh::addPatches(const List<faPatch*>& p)")
            << "boundary already exists"
            << abort(FatalError);
    }

    boundary_.setSize(p.size());

    forAll(p, patchI)
    {
        boundary_.hook(p[patchI]);
    }

    setPrimitiveMeshData();

    boundary_.checkDefinition();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


const objectRegistry& faMesh::db() const
{
    return mesh_.db();
}


const faBoundaryMesh& faMesh::boundary() const
{
    return boundary_;
}


faBoundaryMesh& faMesh::boundary()
{
    return boundary_;
}


const lduAddressing& faMesh::ldu() const
{
    if (!lduPtr_)
    {
        calcLduAddressing();
    }

    return *lduPtr_;
}


const labelList& faMesh::patchStarts() const
{
    if (!patchStartsPtr_)
    {
        calcPatchStarts();
    }

    return *patchStartsPtr_;
}


const edgeVectorField& faMesh::Le() const
{
    if (!LePtr_)
    {
        calcLe();
    }

    return *LePtr_;
}


const edgeScalarField& faMesh::magLe() const
{
    if (!magLePtr_)
    {
        calcMagLe();
    }

    return *magLePtr_;
}


const areaVectorField& faMesh::centres() const
{
    if (!centresPtr_)
    {
        calcCentres();
    }

    return *centresPtr_;
}


const edgeVectorField& faMesh::edgeCentres() const
{
    if (!edgeCentresPtr_)
    {
        calcEdgeCentres();
    }

    return *edgeCentresPtr_;
}


const scalarField& faMesh::S() const
{
    if (!SPtr_)
    {
        calcS();
    }

    return *SPtr_;
}


const scalarField& faMesh::S0() const
{
    if (!S0Ptr_)
    {
        FatalErrorIn("faMesh::S0() const")
            << "S0 is not available"
            << abort(FatalError);
    }

    return *S0Ptr_;
}


const scalarField& faMesh::S00() const
{
    if (!S00Ptr_)
    {
        S00Ptr_ = new scalarIOField
        (
            IOobject
            (
                "S00",
                operator()().pointsInstance(),
                operator()(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            S0()
        );

        S0Ptr_->writeOpt() = IOobject::AUTO_WRITE;
    }

    return *S00Ptr_;
}


const areaVectorField& faMesh::faceAreaNormals() const
{
    if (!faceAreaNormalsPtr_)
    {
        calcFaceAreaNormals();
    }

    return *faceAreaNormalsPtr_;
}


const edgeVectorField& faMesh::edgeAreaNormals() const
{
    if (!edgeAreaNormalsPtr_)
    {
        calcEdgeAreaNormals();
    }

    return *edgeAreaNormalsPtr_;
}


const vectorField& faMesh::pointAreaNormals() const
{
    if (!pointAreaNormalsPtr_)
    {
        calcPointAreaNormals();
    }

    return *pointAreaNormalsPtr_;
}


const areaScalarField& faMesh::faceCurvatures() const
{
    if (!faceCurvaturesPtr_)
    {
        calcFaceCurvatures();
    }

    return *faceCurvaturesPtr_;
}


const FieldField<Field, tensor>& faMesh::edgeTransformTensors() const
{
    if (!edgeTransformTensorsPtr_)
    {
        calcEdgeTransformTensors();
    }

    return *edgeTransformTensorsPtr_;
}


// Return parallel info
const faProcTopology& faMesh::parallelData() const
{
    if (!parallelDataPtr_)
    {
        parallelDataPtr_ =
            new faProcTopology(boundary());
    }

    return *parallelDataPtr_;
}


tmp<scalarField> faMesh::movePoints(const vectorField& newPoints)
{
    moving_ = true;
    
    // Grab old time areas if the time has been incremented
    if (curMotionTimeIndex_ < operator()().time().timeIndex())
    {
        if (S00Ptr_ && S0Ptr_)
        {
            *S00Ptr_ = *S0Ptr_;
        }

        if (S0Ptr_)
        {
            *S0Ptr_ = S();
        }
        else
        {
            S0Ptr_ = new scalarIOField
            (
                IOobject
                (
                    "S0",
                    operator()().pointsInstance(),
                    operator()(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                S()
            );
        }

        curMotionTimeIndex_ = operator()().time().timeIndex();
    }

    clearGeomNotAreas();

    patch().movePoints(newPoints);
    boundary().movePoints(newPoints);
    edgeInterpolation::movePoints();

    tmp<scalarField> tresult
    (
        new scalarField
        (
            nEdges(),
            0.0
        )
    );
    
    return tresult;
}


bool faMesh::write
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    faceLabels_.write(fmt, ver, cmp);
    boundary_.write(fmt, ver, cmp);

    return false;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool faMesh::operator!=(const faMesh& m) const
{
    return &m != this;
}

bool faMesh::operator==(const faMesh& m) const
{
    return &m == this;
}


// * * * * * * * * * * * * * *  IOstream Operators * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const faMesh& m)
{
    os << "faceLabels " << m.faceLabels_ << endl;

    os << endl << m.boundary() << endl;

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
