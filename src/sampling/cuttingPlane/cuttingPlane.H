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
    cuttingPlane

Description

SourceFiles
    cuttingPlane.C

\*---------------------------------------------------------------------------*/

#ifndef cuttingPlane_H
#define cuttingPlane_H

#include <libgen.h>
#include "plane.H"
#include "polyMesh.H"
#include "labelList.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cuttingPlane Declaration
\*---------------------------------------------------------------------------*/

class cuttingPlane
:
    public plane
{
    // Private data

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Point intersections
        vectorField cuttingPoints_;

        //- Cut faces
        faceList cutFaces_;

        //- List of cells cut by the plane
        labelList cutCells_;

        //- Cell cutting points
        labelListList cutCellCuttingPoints_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cuttingPlane(const cuttingPlane&);

        //- Disallow default bitwise assignment
        void operator=(const cuttingPlane&);


        // Helper member functions

            //- Compare scalar sizes
            //- Find cut cells
            void cutCells
            (
                const scalarField& dotProducts, 
                const cellList& CellFaces,
                const faceList& Faces
            );

            //- Find cut edges and points for each cut cell
            //  0 = edges, 1 = points
            labelListList cuttingPoints
            (
                const scalarField& dotProducts,
                const labelListList& CellEdges,
                const pointField& Points
            );

            //- Assign cutting points to cut cells and delete zero faces
            void cutCellCuttingPoints
            (
                const labelListList intersections,
                const labelListList& CellEdges,
                const cellList& CellFaces,
                const faceList& Faces
            );

            //- Order cutting point labels for each face
            void orderFacePoints();

            //- Display point list
            void showPoints();


            //- Container class for sorting
            struct sl
            {
                scalar s_;
                label  l_;

                friend bool operator<(const sl& one, const sl& two)
                {
                    return one.s_ < two.s_;
                }
            };

            typedef List<sl> slList; 


public:

    // Constructors

        //- Construct from components: Mesh reference and plane
        cuttingPlane(const polyMesh& mesh, const plane&);


    // Member Functions

        //- Return vectorField of cutting points
        const pointField& points() const;

        //- Return faceList of points in cells
        const faceList& faces() const;

        //- Return labelList of cut cells
        const labelList& cells() const;

        //- Return true or false to question: have any cells been cut?
        bool cut();

        //- Sample the cell field
        template<class Type>
        tmp<Field<Type> > sample(const Field<Type>&) const;

        template<class Type>
        tmp<Field<Type> > sample(const tmp<Field<Type> >&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cuttingPlaneSample.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
