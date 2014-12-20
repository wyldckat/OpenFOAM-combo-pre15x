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
    constantPlane

Description

SourceFiles
    constantPlane.C

\*---------------------------------------------------------------------------*/

#ifndef constantPlane_H
#define constantPlane_H

#include "surface.H"
#include "cuttingPlane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class meshSearch;
class dictionary;
class plane;
class volPointInterpolation;
template<class Type> class fieldsCache;

/*---------------------------------------------------------------------------*\
                           Class constantPlane Declaration
\*---------------------------------------------------------------------------*/

class constantPlane
:
    public surface
{
    // Private data

        //- Cutting plane description
        const plane planeDesc_;

        //- Make triangles or keep faces
        const bool triangulate_;

        //- Generated points
        pointField points_;

        //- Generated faces
        faceList faces_;

        //- Fro every face original cell in mesh
        labelList meshCells_;


    // Private Member Functions

        //- Cut & triangulate
        void makeTriangles(const cuttingPlane&);

        //- Copy cut; do not triangulate
        void copyFaces(const cuttingPlane&);

        //- Do all to create geometry. Cut&triangulate if nessecary
        void createGeometry();

        //- Interpolate field onto faces
        template <class T>
        Field<T> doInterpolate
        (
            const word& fieldName,
            const fieldsCache<T>& cache
        ) const;

public:

    //- Runtime type information
    TypeName("constantPlane");

    // Constructors

        //- Construct from components
        constantPlane
        (
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const word& name,
            const plane& planeDesc,
            const bool triangulate = true
        );

        //- Construct from dictionary
        constantPlane
        (
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const dictionary& dict          
        );


    // Destructor

        virtual ~constantPlane();


    // Member Functions

        const plane& planeDesc() const
        {
            return planeDesc_;
        }

        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- Points of surface
        virtual const pointField& points() const
        {
            return points_;
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return faces_;
        }

        //- Correct for mesh movement and/or field changes
        virtual void correct
        (
            const bool meshChanged,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes,
            const fieldsCache<scalar>&,
            const fieldsCache<vector>&,
            const fieldsCache<tensor>&
        );

        //- interpolate field to surface
        virtual scalarField interpolate
        (
            const word&,
            const fieldsCache<scalar>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;

        //- interpolate field to surface
        virtual vectorField interpolate
        (
            const word&,
            const fieldsCache<vector>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;

        //- interpolate field to surface
        virtual tensorField interpolate
        (
            const word&,
            const fieldsCache<tensor>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
