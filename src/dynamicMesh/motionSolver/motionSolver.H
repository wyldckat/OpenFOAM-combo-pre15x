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
    motionSolver

Description
    Mesh motion solver for a polyMesh.  Based on solving the
    vertex-based motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable motionU.

    For simple situations, call newPoints() member directly (this
    internally solves the equations).  curPoints() and solve() members
    are added for additional manipulation of the solution.

SourceFiles
    motionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef motionSolver_H
#define motionSolver_H

#include "polyMesh.H"
#include "tetPolyMesh.H"
#include "tetPointFields.H"
#include "elementFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class twoDPointCorrector;

/*---------------------------------------------------------------------------*\
                           Class motionSolver Declaration
\*---------------------------------------------------------------------------*/

class motionSolver
:
    public IOdictionary
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Motion mesh
        tetPolyMesh tetMesh_;

        //- Motion field
        tetPointVectorField motionU_;

        //- Frozen diffusion field
        elementScalarField* gammaPtr_;

        //- 2-D motion corrector pointer
        twoDPointCorrector* correct2DPtr_;

        //- First motion
        bool firstMotion_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        motionSolver(const motionSolver&);

        //- Disallow default bitwise assignment
        void operator=(const motionSolver&);


public:

    // Static data members

    // Constructors

        //- Construct from components
        motionSolver(const polyMesh&);


    // Destructor

        ~motionSolver();


    // Member Functions

        //- Return reference to motion field
        tetPointVectorField& motionU()
        {
            return motionU_;
        }

        //- Return reference to motion diffusivity
        elementScalarField& gamma();

        //- Freeze diffusivity.  Return true if freezing and false if
        //  already frozen
        bool freezeGamma();

        //- Release diffusivity.  Return true if releasing and false
        //  if already released
        bool releaseGamma();

        //- Provide new points for motion.  Solves for motion
        tmp<pointField> newPoints();

        //- Provide current points for motion.  Uses current motion field
        tmp<pointField> curPoints() const;

        //- Solve for motion
        void solve();


        // Mesh quality measures

            //- Return mean cell non-orthogonality
            tmp<scalarField> nonOrthogonality() const;

            //- Return mean cell skewness
            tmp<Foam::scalarField> cellSkewness() const;


        //- Update topology
        void updateTopology();

        //- Read motionProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
