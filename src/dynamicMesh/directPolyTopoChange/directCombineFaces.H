/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    directCombineFaces

Description
    Combines boundary faces into single face. The faces get the patch
    of the first face ('the master')

SourceFiles
    directCombineFaces.C

\*---------------------------------------------------------------------------*/

#ifndef directCombineFaces_H
#define directCombineFaces_H

#include "typeInfo.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;
class directPolyTopoChange;
class mapPolyMesh;
class face;

/*---------------------------------------------------------------------------*\
                           Class directCombineFaces Declaration
\*---------------------------------------------------------------------------*/

class directCombineFaces
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;


    // Private Member Functions

        //- Mark all faces connected through non-feature edges to faceI 
        //  with the current region
        static void regionFaces
        (
            const indirectPrimitivePatch&,
            const scalar featureCos,
            const label regionI,            // current region
            const label faceI,
            labelList& region               // mark array
        );

        //- Test if face is convex. Allow slight concavity through concaveSin.
        static bool convexFace
        (
            const scalar concaveSin,
            const pointField&,
            const face&
        );

        //- Gets outside of patch as a face (in mesh point labels)
        static face getOutsideFace(const indirectPrimitivePatch&);

        //- Disallow default bitwise copy construct
        directCombineFaces(const directCombineFaces&);

        //- Disallow default bitwise assignment
        void operator=(const directCombineFaces&);

public:

    //- Runtime type information
    ClassName("directCombineFaces");


    // Constructors

        //- Construct from mesh
        directCombineFaces(const polyMesh& mesh);


    // Member Functions

        // Helper functions

            //- Extract lists of connected faces (might be only one face)
            //  on a single cell.
            //  The input faces have to be connected through non-feature edges,
            //  be on the same cell and on the same patch.
            //  featureCos : cos of angle between connected faces. Any face
            //               with cos > featureCos will be included.
            //  concaveSin : max allowed angle between consecutive edges for
            //               resulting face to be considered convex (and valid).
            //               Set to < -1 to disable.
            labelListList getFaceRegions
            (
                const scalar featureCos,
                const scalar concaveSin,
                const labelList& faceLabels
            ) const;


            //- Extract lists of all (non-coupled) boundary faces that can
            //  be merged. Uses getFaceRegions.
            labelListList getMergeSets
            (
                const scalar featureCos,
                const scalar concaveSin
            ) const;


        // Topology changes

            //- Play commands into polyTopoChange to combine faces. Gets
            //  labelListList of sets of faces to combine. Does no check
            //  for whether resulting face is legal.
            void setRefinement(const labelListList&, directPolyTopoChange&);

            //- Force recalculation of locally stored data on topological change
            void updateMesh(const mapPolyMesh&)
            {}

        // Access


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
