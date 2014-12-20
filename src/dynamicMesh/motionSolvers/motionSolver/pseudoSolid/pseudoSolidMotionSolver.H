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
    pseudoSolidMotionSolver

Description
    Mesh motion solver for a polyMesh.  Based on solving the
    vertex-based pseudo-solid motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable motionU.

    Note: Diffusion field from laplaceMotionSolver is used 
          as Lame's constant \mu

SourceFiles
    pseudoSolidMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef pseudoSolidMotionSolver_H
#define pseudoSolidMotionSolver_H

#include "laplaceMotionSolver.H"

#include "tetPolyMesh.H"
#include "tetPointFields.H"
#include "elementFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pseudoSolidMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class pseudoSolidMotionSolver
:
    public laplaceMotionSolver
{
    // Private data

        //- Poisson's ration
        scalar nu_;

        //- Number of correctors in solution procedure of pseudo-solid eqn
        label nCorrectors_;

        //- Convergence tolerance for solution of pseudo-solid eqn
        scalar convergenceTolerance_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pseudoSolidMotionSolver(const pseudoSolidMotionSolver&);

        //- Disallow default bitwise assignment
        void operator=(const pseudoSolidMotionSolver&);


public:

    //- Runtime type information
    TypeName("pseudoSolid");


    // Constructors

        pseudoSolidMotionSolver(const polyMesh& mesh);


    // Destructor

        ~pseudoSolidMotionSolver();


    // Member Functions

        //- Solve for motion
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
