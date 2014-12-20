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
    laplaceMotionSolver

Description
    Mesh motion solver for a polyMesh.  Based on solving the
    vertex-based Laplace motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable motionU.

SourceFiles
    laplaceMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef laplaceMotionSolver_H
#define laplaceMotionSolver_H

#include "motionSolver.H"
#include "tetPolyMesh.H"
#include "tetPointFields.H"
#include "elementFields.H"
#include "lduMatrix.H"
#include "motionDiff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiff;

/*---------------------------------------------------------------------------*\
                           Class laplaceMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class laplaceMotionSolver
:
    public motionSolver
{
    // Private data

        //- Diffusion
        motionDiff* diffusionPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        laplaceMotionSolver(const laplaceMotionSolver&);

        //- Disallow default bitwise assignment
        void operator=(const laplaceMotionSolver&);

protected:

        //- Is diffusion field frozen
        bool frozen_;

        //- First motion
        bool firstMotion_;

        //- Return pointer to motion diffusivity
        motionDiff* diffusion()
        {
            return diffusionPtr_;
        }

        // Holds recent solver performance
        lduMatrix::solverPerformance solverPerf_;

public:

    //- Runtime type information
    TypeName("laplace");


    // Constructors

        laplaceMotionSolver(const polyMesh& mesh);


    // Destructor

        ~laplaceMotionSolver();


    // Member Functions

        //- Solve for motion
        virtual void solve();

        //- Return recent solver performance
        const lduMatrix::solverPerformance& solverPerformance() const
        {
            return solverPerf_;
        }
    
        //- Return pointer to motion diffusivity
        const elementScalarField& gamma() const
        {
            return diffusionPtr_->gamma();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
