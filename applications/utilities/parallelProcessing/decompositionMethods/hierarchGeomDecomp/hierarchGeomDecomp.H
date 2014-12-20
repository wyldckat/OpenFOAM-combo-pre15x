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
    hierarchGeomDecomp

Description
    Does hierarchical decomposition of points. Works by first sorting the
    points in x direction into equal sized bins, then in y direction and
    finally in z direction.

    Uses single array to hold decomposition which is indexed as if it is a
    3 dimensional array:

        finalDecomp[i,j,k] is indexed as

        i*n[0]*n[1] + j*n[1] + k

    E.g. if we're sorting 'xyz': the first sort (over the x-component)
    determines in which x-domain the point goes. Then for each of the x-domains
    the points are sorted in y direction and each individual x-domain gets
    split into three y-domains. And similar for the z-direction.

    Since the domains are of equal size the maximum difference in size is
    n[0]*n[1] (or n[1]*n[2]?) (small anyway)


SourceFiles
    hierarchGeomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef hierarchGeomDecomp_H
#define hierarchGeomDecomp_H

#include "geomDecomp.H"
#include "FixedList.H"
#include "direction.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hierarchGeomDecomp Declaration
\*---------------------------------------------------------------------------*/

class hierarchGeomDecomp
:
    public geomDecomp
{
    // Private data

        //- Decomposition order in terms of components.
        FixedList<direction, 3> decompOrder_;


    // Private Member Functions

        //- Convert ordering string ("xyz") into list of components.
        void setDecompOrder();

        //- Recursively sort in x,y,z (or rather acc. to decompOrder_)
        void sortComponent
        (
            const pointField&,
            const labelList& slice,         // slice of points to decompose
            const direction componentIndex, // index in decompOrder_
            const label prevMult,           // multiplication factor
            labelList& finalDecomp          // overall decomposition
        );


        //- Disallow default bitwise copy construct and assignment
        void operator=(const hierarchGeomDecomp&);
        hierarchGeomDecomp(const hierarchGeomDecomp&);


public:

    //- Runtime type information
    TypeName("hierarchical");


    // Constructors

        //- Construct given the decomposition dictionary
        hierarchGeomDecomp(const dictionary& decompositionDict);

        //- Construct given the decomposition dictionary and mesh
        hierarchGeomDecomp
        (
            const dictionary& decompositionDict,
            const primitiveMesh& mesh
        );


    // Destructor

        ~hierarchGeomDecomp()
        {}


    // Member Functions

        labelList decompose(const pointField&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
