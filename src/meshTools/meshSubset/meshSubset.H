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
    meshSubset

Description
    Post-processing mesh subset tool.  Given the original mesh and the
    list of selected cells, it creates the mesh consisting only of the
    desired cells, with the mapping list for points, faces, and cells.

    Puts all exposed internal faces into either
    - a user supplied patch
    - a newly created patch "oldInternalFaces"

SourceFiles
    meshSubset.C

\*---------------------------------------------------------------------------*/

#ifndef meshSubset_H
#define meshSubset_H

#include "fvMesh.H"
#include "fvPatchFieldMapper.H"
#include "GeometricField.H"
#include "emptyFvPatchFields.H"
#include "labelHashSet.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class meshSubset Declaration
\*---------------------------------------------------------------------------*/

class meshSubset
:
    public fvMesh
{

public:

    //- Patch-field subset interpolation class
    class patchFieldSubset
    :
        public fvPatchFieldMapper
    {
        labelField directAddressing_;

    public:

        // Constructors

            //- Construct given addressing
            patchFieldSubset(const meshSubset& ms, const label patchNo)
            :
                directAddressing_
                (
                    (const labelField&)
                    (
                        labelField::subField
                        (
                            ms.faceMap(),
                            ms.subMesh().boundary()[patchNo].size(),
                            ms.subMesh().boundary()[patchNo].patch().start()
                        )
                    )
                  - ms.boundary()[ms.patchMap()[patchNo]].patch().start()
                )
            {}


        // Destructor

            ~patchFieldSubset()
            {}


        // Member Functions

            label size() const
            {
                return directAddressing_.size();
            }

            virtual bool resizeOnly() const
            {
                return false;
            }

            bool direct() const
            {
                return true;
            }

            const unallocLabelList& directAddressing() const
            {
                return directAddressing_;
            }
    };


private:

    // Private data

        //- Subset mesh pointer
        fvMesh* meshSubsetPtr_;

        //- Point mapping array
        labelList pointMap_;

        //- Face mapping array
        labelList faceMap_;

        //- Cell mapping array
        labelList cellMap_;

        //- Patch mapping array
        labelList patchMap_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshSubset(const meshSubset&);

        //- Disallow default bitwise assignment
        void operator=(const meshSubset&);

        //- Check if subset has been performed
        bool checkCellSubset() const;

        //- Mark points in Map
        static void markPoints(const labelList&, Map<label>&); 


public:

    // Constructors

        //- Construct given a mesh and the map of cells to subset
        explicit meshSubset(const IOobject& io);


    // Destructor

        ~meshSubset();


    // Member Functions

        // Edit

            //- Set the subset. Create "oldInternalFaces" patch for exposed
            //  internal faces (patchID==-1) or use supplied patch.
            void setCellSubset
            (
                const labelHashSet& globalCellMap,
                const label patchID = -1
            );


        // Access

            //- Return reference to subset mesh
            const fvMesh& subMesh() const;

            //- Return point map
            const labelList& pointMap() const;

            //- Return face map
            const labelList& faceMap() const;

            //- Return cell map
            const labelList& cellMap() const;

            //- Return patch map
            const labelList& patchMap() const;


        // Field mapping

            //- Map volume field
            template<class Type>
            tmp<GeometricField<Type, fvPatchField, volMesh> > interpolate
            (
                const GeometricField<Type, fvPatchField, volMesh>&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "meshSubsetInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
