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
    volPointInterpolation

Description

SourceFiles
    volPointInterpolation.C
    volPointInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef volPointInterpolation_H
#define volPointInterpolation_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "scalarList.H"
#include "tmp.H"
#include "className.H"
#include "FieldFields.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class pointMesh;

/*---------------------------------------------------------------------------*\
                       Class volPointInterpolation Declaration
\*---------------------------------------------------------------------------*/

class volPointInterpolation
{
    // Private data

        const fvMesh& fvMesh_;
        const pointMesh& pointMesh_;

        //- Interpolation scheme weighting factor array.
        mutable FieldField<Field, scalar>* pointWeightsPtr_;

        //- Primitive patch interpolators
        mutable PtrList<primitivePatchInterpolation>* patchInterpolatorsPtr_;

        //- List of points that will be corrected
        mutable labelList* boundaryPointsPtr_;

        //- Weights
        mutable FieldField<Field, scalar>* pointBoundaryWeightsPtr_;


    // Private member functions

        //- Return patch interpolators
        const PtrList<primitivePatchInterpolation>& patchInterpolators() const;

        //- Construct point weighting factors
        void makeWeights() const;

        //- Construct boundary point addressing
        void makeBoundaryAddressing() const;

        //- Construct boundary point weighting factors
        void makeBoundaryWeights() const;

        //- Clear addressing
        void clearAddressing() const;

        //- Clear geometry
        void clearGeom() const;


protected:

        const fvMesh& vMesh() const
        {
            return fvMesh_;
        }

        const pointMesh& pMesh() const
        {
            return pointMesh_;
        }


public:

    // Declare name of the class and it's debug switch
    ClassName("volPointInterpolation");


    // Constructors

        //- Constructor given fvMesh and pointMesh.
        volPointInterpolation(const fvMesh&, const pointMesh&);


    // Destructor

        ~volPointInterpolation();


    // Member functions

        // Access

            //- Return reference to weights arrays.
            //  This also constructs the weighting factors if neccessary.
            const FieldField<Field, scalar>& pointWeights() const;

            //- Return list of points to correct
            const labelList& boundaryPoints() const;

            //- Return reference to boundary weights arrays.
            const FieldField<Field, scalar>& pointBoundaryWeights() const;


        // Edit

            //- Update mesh topology using the morph engine
            void updateTopology();

            //- Correct weighting factors for moving mesh.
            bool movePoints();


    // Interpolation functions

        //- Interpolate from volField to pointField
        //  using inverse distance weighting
        template<class Type>
        void interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Interpolate volField returning pointField
        //  using inverse distance weighting
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Interpolate tmp<volField> returning pointField
        //  using inverse distance weighting
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "volPointInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
