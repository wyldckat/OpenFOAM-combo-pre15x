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

Class
    fvMesh

Description
    Mesh data needed to do the Finite Volume discretisation.

    NOTE ON USAGE:
    fvMesh contains all the topological and geometric information
    related to the mesh.  It is also responsible for keeping the data
    up-to-date.  This is done by deleting the cell volume, face area,
    cell/face centre, addressing and other derived information as
    required and recalculating it as necessary.  The fvMesh therefore
    reserves the right to delete the derived information upon every
    topological (mesh refinement/morphing) or geometric change (mesh
    motion).  It is therefore unsafe to keep local references to the
    derived data outside of the time loop.

SourceFiles
    fvMesh.C
    makeFvMeshData.C

\*---------------------------------------------------------------------------*/

#ifndef fvMesh_H
#define fvMesh_H

#include "polyMesh.H"
#include "primitiveMesh.H"
#include "fvBoundaryMesh.H"
#include "surfaceInterpolation.H"
#include "scalarIOField.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "lduAddressing.H"
#include "className.H"

#define BAD_MESH_STABILISATION

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduAddressingFvMesh;

/*---------------------------------------------------------------------------*\
                           Class fvMesh Declaration
\*---------------------------------------------------------------------------*/

class fvMesh
:
    public polyMesh,
    public surfaceInterpolation
{
    // Private data

        //- Boundary mesh
        fvBoundaryMesh boundary_;


    // Demand-driven data

        mutable lduAddressingFvMesh* lduPtr_;

        //- Current time index for cell volumes
        //  Note.  The whole mechanism will be replaced once the
        //  dimensionedField is created and the dimensionedField
        //  will take care of the old-time levels.  
        mutable label curTimeIndex_;

        //- Cell volumes old time level
        mutable scalarIOField* V0Ptr_;

        //- Cell volumes old-old time level
        mutable scalarIOField* V00Ptr_;

        //- Face area vectors
        mutable surfaceVectorField* SfPtr_;

        //- Mag face area vectors
        mutable surfaceScalarField* magSfPtr_;

        //- Cell centres
        mutable volVectorField* CPtr_;

        //- Face centres
        mutable surfaceVectorField* CfPtr_;

        //- Face motion fluxes
        mutable surfaceScalarField* phiPtr_;


    // Private Member Functions

        // Storage management

            //- Clear geometry but not the cell volumes
            void clearGeomNotVol();

            //- Clear geometry
            void clearGeom();

            //- Clear addressing
            void clearAddressing();

            //- Clear all geometry and addressing
            void clearOut();


       // Make geometric data

            void makeSf() const;
            void makeMagSf() const;
            void makePhi() const;

            void makeC() const;
            void makeCf() const;


        //- Disallow construct as copy
        fvMesh(const fvMesh&);

        //- Disallow assignment
        void operator=(const fvMesh&);

protected:

        //- Map all fields in time using current polyMesh::morphMap.
        virtual void mapFields();

        //- Update mesh info after topo change
        virtual void handleMorph();

        //- Remove boundary patches. Warning: fvPatchFields hold ref to
        //  these fvPatches.
        void removeFvBoundary();

        //- Return cell face motion fluxes
        surfaceScalarField& setPhi();

        //- Return old-time cell volumes
        scalarField& setV0();

public:

    // Declare name of the class and it's debug switch
    ClassName("fvMesh");


    // Constructors

        //- Construct from IOobject
        explicit fvMesh(const IOobject& io);

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        fvMesh
        (
            const IOobject& io,
            const pointField& points,
            const faceList& faces,
            const cellList& cells
        );

        //- Construct from cell shapes
        fvMesh
        (
            const IOobject& io,
            const pointField& points,
            const cellShapeList& shapes,
            const faceListList& boundaryFaces,
            const wordList& boundaryPatchNames,
            const wordList& boundaryPatchTypes,
            const word& defaultBoundaryPatchType,
            const wordList& boundaryPatchPhysicalTypes
        );


    // Destructor

        virtual ~fvMesh();


    // Member Functions

        // Helpers

            //- Add boundary patches. Constructor helper
            void addFvPatches(const List<polyPatch* > &);

            //- Update the mesh based on the mesh files saved in time
            //  directories
            virtual readUpdateState readUpdate();

        // Access

            //- Return reference to DB
            const Time& time() const
            {
                return polyMesh::time();
            }

            //- Return reference to name 
            //  Note: name() is currently ambiguous due to derivation from
            //  surfaceInterpolation
            const word& name() const
            {
                return polyMesh::name();
            }

            //- Return reference to boundary mesh
            const fvBoundaryMesh& boundary() const;

            //- Return cell volumes
            const scalarField& V() const;

            //- Return old-time cell volumes
            const scalarField& V0() const;

            //- Return old-old-time cell volumes
            const scalarField& V00() const;

            //- Return cell face area vectors
            const surfaceVectorField& Sf() const;

            //- Return cell face area magnitudes
            const surfaceScalarField& magSf() const;

            //- Return cell face motion fluxes
            const surfaceScalarField& phi() const;

            //- Return cell centres as volVectorField
            const volVectorField& C() const;

            //- Return face centres as surfaceVectorField
            const surfaceVectorField& Cf() const;

            //- Return ldu addressing
            const lduAddressing& ldu() const;

            // For consistency with older nomenclature, provide
            // wrapping for owner/neighbour addressing from fvMesh
            // 

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

            //- Construct all optional data and clear basic info from basicMesh
            void constructAndClear() const;

            //- Update mesh topology using the morph engine
            virtual void updateTopology();

            //- Move points, returns volumes swept by faces in motion
            virtual tmp<scalarField> movePoints(const vectorField&);


        // Write

            //- Write the underlying polyMesh and other data
            virtual bool write
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp
            ) const;

            //- Write mesh using IO settings from time
            virtual bool write() const;


    // Member Operators

        bool operator!=(const fvMesh&) const;
        bool operator==(const fvMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
