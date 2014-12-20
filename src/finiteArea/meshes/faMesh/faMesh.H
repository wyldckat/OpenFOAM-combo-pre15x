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

Class
    faMesh

Description
    Finite area mesh.  Used for 2-D non-Euclidian finite area method.

SourceFiles
    faMesh.C
    faMeshDemandDrivenData.C

\*---------------------------------------------------------------------------*/

#ifndef faMesh_H
#define faMesh_H

#include "faBoundaryMesh.H"
#include "faProcTopology.H"
#include "edgeList.H"
#include "faceList.H"
#include "primitiveFieldsFwd.H"
#include "lduAddressing.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "indirectPrimitivePatch.H"
#include "edgeInterpolation.H"
#include "labelIOList.H"
#include "scalarIOField.H"
#include "FieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;
class lduAddressingFaMesh;

/*---------------------------------------------------------------------------*\
                           Class faMesh Declaration
\*---------------------------------------------------------------------------*/

class faMesh      
:
    public edgeInterpolation
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Face labels
        labelIOList faceLabels_;

        //- Boundary mesh
        faBoundaryMesh boundary_;


        // Primitive mesh data

            //- Edges, addressing into local point list
            edgeList edges_;

            //- Edge owner
            labelList edgeOwner_;

            //- Edge neighbour
            labelList edgeNeighbour_;


        // Primitive size data

            //- Number of points
            mutable label nPoints_;

            //- Number of edges
            mutable label nEdges_;

            //- Number of internal edges
            mutable label nInternalEdges_;

            //- Number of faces
            mutable label nFaces_;


    // Demand-driven data

        //- Primitive patch
        mutable indirectPrimitivePatch* patchPtr_;

        //- Ldu addressing data
        mutable lduAddressingFaMesh* lduPtr_;

        //- Face areas
        mutable scalarField* SPtr_;

        //- Face areas old time level
        mutable scalarIOField* S0Ptr_;

        //- Face areas old-old time level
        mutable scalarIOField* S00Ptr_;

        //- Patch starts in the edge list
        mutable labelList* patchStartsPtr_;

        //- Edge length vectors
        mutable edgeVectorField* LePtr_;

        //- Mag edge length vectors
        mutable edgeScalarField* magLePtr_;

        //- Face centres
        mutable areaVectorField* centresPtr_;

        //- Edge centres
        mutable edgeVectorField* edgeCentresPtr_;

        //- Face area normals
        mutable areaVectorField* faceAreaNormalsPtr_;

        //- Edge area normals
        mutable edgeVectorField* edgeAreaNormalsPtr_;

        //- Edge area normals
        mutable vectorField* pointAreaNormalsPtr_;

        //- Face curvatures
        mutable areaScalarField* faceCurvaturesPtr_;

        //- Edge transformation tensors
        mutable FieldField<Field, tensor>* edgeTransformTensorsPtr_;


        // Other mesh related data

            //- Parallel info
            mutable faProcTopology* parallelDataPtr_;

            // Mesh motion

                //- Is the mesh moving
                bool moving_;

                //- Current time index for mesh motion
                mutable label curMotionTimeIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faMesh(const faMesh&);

        //- Disallow default bitwise assignment
        void operator=(const faMesh&);

        //- Set primitive mesh data
        void setPrimitiveMeshData();

        // Private member functions to calculate demand driven data

            //- Calculate ldu addressing
            void calcLduAddressing() const;

            //- Calculate patch starts in the edge list
            void calcPatchStarts() const;

            //- Calculate edge lengths
            void calcLe() const;

            //- Calculate mag edge lengths
            void calcMagLe() const;

            //- Calculate face centres
            void calcCentres() const;

            //- Calculate edge centres
            void calcEdgeCentres() const;

            //- Calculate face areas
            void calcS() const;

            //- Calculate face area normals
            void calcFaceAreaNormals() const;

            //- Calculate edge area normals
            void calcEdgeAreaNormals() const;

            //- Calculate point area normals
            void calcPointAreaNormals() const;

            //- Calculate face curvatures
            void calcFaceCurvatures() const;

            //- Calculate edge transformation tensors
            void calcEdgeTransformTensors() const;

            //- Clear geometry but not the face areas
            void clearGeomNotAreas() const;

            //- Clear geometry
            void clearGeom() const;

            //- Clear addressing
            void clearAddressing() const;

            //- Clear demand-driven data
            void clearOut() const;

public:

    // Public typedefs

        typedef faMesh Mesh;
        typedef faBoundaryMesh BoundaryMesh;


    // Static data members

    // Declare name of the class and it's debug switch
    TypeName("faMesh");

    //- Return the mesh sub-directory name (usually "faMesh")
    static word meshSubDir;


    // Constructors

        //- Construct from objectRegistry, polyMesh reference
        //  and read/write options
        faMesh
        (
            const polyMesh& m,
            IOobject::readOption r=IOobject::MUST_READ,
            IOobject::writeOption w=IOobject::AUTO_WRITE
        );

        //- Construct from components without boundary.
        //  Boundary is added using addFaPatches() member function
        faMesh
        (
            const polyMesh& m,
            const labelList& faceLabels,
            IOobject::readOption r=IOobject::NO_READ,
            IOobject::writeOption w=IOobject::AUTO_WRITE
        );


    // Destructor

        virtual ~faMesh();


    // Member Functions

        // Database

            //- Return the local mesh directory (dbDir()/meshSubDir)
            fileName meshDir() const;

            //- Return reference to time
            const Time& time() const;


            //- Mesh size parameters

                inline label nPoints() const
                {
                    return nPoints_;
                }

                inline label nEdges() const
                {
                    return nEdges_;
                }

                inline label nInternalEdges() const
                {
                    return nInternalEdges_;
                }

                inline label nFaces() const
                {
                    return nFaces_;
                }

            // Primitive mesh data

                //- Return mesh points
                const pointField& points() const;

                //- Return edges
                const edgeList& edges() const;

                //- Return faces
                const faceList& faces() const;

                //- Edge owner addresing
                inline const labelList& edgeOwner() const
                {
                    return edgeOwner_;
                }

                //- Edge neighbour addressing
                inline const labelList& edgeNeighbour() const
                {
                    return edgeNeighbour_;
                }



        //- Add boundary patches. Constructor helper
        void addFaPatches(const List<faPatch*> &);

        //- Return reference to the mesh database
        const objectRegistry& db() const;

        //- Return constant reference to boundary mesh
        const faBoundaryMesh& boundary() const;

        //- Return reference to boundary mesh
        faBoundaryMesh& boundary();

        //- Return faMesh face labels
        const labelList& faceLabels() const
        {
            return faceLabels_;
        }


        //- Return parallel info
        const faProcTopology& parallelData() const;


        // Mesh motion

            //- Is mesh moving
            bool moving() const
            {
                return moving_;
            }

            //- Move points, returns volumes swept by faces in motion
            virtual tmp<scalarField> movePoints(const vectorField&);


        // Demand-driven data

            //- Return constant reference to primitive patch
            const indirectPrimitivePatch& patch() const;

            //- Return reference to primitive patch
            indirectPrimitivePatch& patch();

            //- Return ldu addressing
            const lduAddressing& ldu() const;

            //- Return patch starts
            const labelList& patchStarts() const;

            //- Return edge length vectors
            const edgeVectorField& Le() const;

            //- Return edge length magnitudes
            const edgeScalarField& magLe() const;

            //- Return face centres as areaVectorField
            const areaVectorField& centres() const;

            //- Return edge centres as edgeVectorField
            const edgeVectorField& edgeCentres() const;

            //- Return face areas
            const scalarField& S() const;

            //- Return old-time face areas
            const scalarField& S0() const;

            //- Return old-old-time face areas
            const scalarField& S00() const;

            //- Return face area normals
            const areaVectorField& faceAreaNormals() const;

            //- Return edge area normals
            const edgeVectorField& edgeAreaNormals() const;

            //- Return point area normals
            const vectorField& pointAreaNormals() const;

            //- Return face curvatures
            const areaScalarField& faceCurvatures() const;

            //- Return edge transformation tensors
            const FieldField<Field, tensor>& edgeTransformTensors() const;

            //- Return internal point labels
            labelList internalPoints() const;

            //- Return boundary point labels
            labelList boundaryPoints() const;

            //- Return edge length correction
            tmp<edgeScalarField> edgeLengthCorrection() const;

            // For consistency with older nomenclature, provide
            // wrapping for owner/neighbour addressing from faMesh

            //- Internal face owner
            const unallocLabelList& owner() const
            {
                return ldu().lowerAddr();
            }

            //- Internal face neighbour
            const unallocLabelList& neighbour() const
            {
                return ldu().upperAddr();
            }


        // Edit

            //- Write mesh
            bool write
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp
            ) const;


    // Member Operators

        //- Return reference to polyMesh
        const polyMesh& operator()() const
        {
            return mesh_;
        }

        bool operator!=(const faMesh& m) const;

        bool operator==(const faMesh& m) const;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

    friend Ostream& operator<<(Ostream&, const faMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
