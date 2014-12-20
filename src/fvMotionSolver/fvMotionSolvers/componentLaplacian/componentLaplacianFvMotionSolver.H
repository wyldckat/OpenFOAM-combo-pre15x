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
    componentLaplacianFvMotionSolver

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    Laplacian motion equation for the given component.

SourceFiles
    componentLaplacianFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef componentLaplacianFvMotionSolver_H
#define componentLaplacianFvMotionSolver_H

#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiffusivity;

/*---------------------------------------------------------------------------*\
                           Class componentLaplacianFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class componentLaplacianFvMotionSolver
:
    public fvMotionSolver
{
    // Private data

        //- The component name to solve for
        word cmptName_;

        //- The component to solve for
        direction cmpt_;

        //- Cell-centre motion field
        mutable volScalarField cellMotionU_;

        //- Point motion field
        mutable pointScalarField pointMotionU_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        componentLaplacianFvMotionSolver
        (
            const componentLaplacianFvMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const componentLaplacianFvMotionSolver&);


public:

    //- Runtime type information
    TypeName("componentLaplacian");


    // Constructors

        componentLaplacianFvMotionSolver(const polyMesh& mesh, Istream& msData);


    // Destructor

        ~componentLaplacianFvMotionSolver();


    // Member Functions

        //- Non-const access to the cellMotionU in order to allow changes to the
        //  boundary motion
        volScalarField& cellMotionU()
        {
            return cellMotionU_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
