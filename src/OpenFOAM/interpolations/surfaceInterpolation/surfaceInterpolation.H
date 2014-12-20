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
    surfaceInterpolation

Description
    Cell to surface interpolation scheme. Included in fvMesh.

SourceFiles
    surfaceInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceInterpolation_H
#define surfaceInterpolation_H

#include "tmp.H"
#include "scalar.H"
#include "fvSchemes.H"
#include "fvSolution.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class surfaceInterpolation Declaration
\*---------------------------------------------------------------------------*/

class surfaceInterpolation
:
    public fvSchemes,
    public fvSolution
{
    // Private data

        // Reference to fvMesh
        const fvMesh& mesh_;

        // Demand-driven data

            //- Central-differencing weighting factors
            mutable surfaceScalarField* weightingFactors_;

            //- Face-gradient difference factors
            mutable surfaceScalarField* differenceFactors_;

            //- Is mesh orthogonal
            mutable bool orthogonal_;

            //- Non-orthogonality correction vectors
            mutable surfaceVectorField* correctionVectors_;

            //- Is mesh skew
            mutable bool skew_;

            //- Skew correction vectors
            mutable surfaceVectorField* skewCorrectionVectors_;

            //- Least-squares gradient vectors
            mutable surfaceVectorField* leastSquarePvectors_;
            mutable surfaceVectorField* leastSquareNvectors_;


    // Private member functions

        //- Construct central-differencing weighting factors
        void makeWeights() const;

        //- Construct face-gradient difference factors
        void makeDeltaCoeffs() const;

        //- Construct non-orthogonality correction vectors
        void makeCorrectionVectors() const;

        //- Construct skewness correction vectors
        void makeSkewCorrectionVectors() const;

        //- Construct Least-squares gradient vectors
        void makeLeastSquareVectors() const;


protected:

    // Protected member functions

        // Storage management

            //- Clear all geometry and addressing
            void clearOut();


public:

    // Declare name of the class and it's debug switch
    ClassName("surfaceInterpolation");


    // Constructors

        //- Construct given an fvMesh
        explicit surfaceInterpolation(const fvMesh&);


    // Destructor

        ~surfaceInterpolation();


    // Member functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to weighting factors array
        const surfaceScalarField& weights() const;

        //- Return reference to difference factors array
        const surfaceScalarField& deltaCoeffs() const;

        //- Return whether mesh is orthogonal or not
        bool orthogonal() const;

        //- Return reference to non-orthogonality correction vectors array
        const surfaceVectorField& correctionVectors() const;

        //- Return whether mesh is skew or not
        bool skew() const;

        //- Return reference to skew vectors array
        const surfaceVectorField& skewCorrectionVectors() const;

        //- Return reference to owner least square vectors
        const surfaceVectorField& leastSquarePvectors() const;

        //- Return reference to neighbour least square vectors
        const surfaceVectorField& leastSquareNvectors() const;

        //- Do what is neccessary if the mesh has moved
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
