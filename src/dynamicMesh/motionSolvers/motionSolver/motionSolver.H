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
    motionSolver

Description
    Virtual base class for mesh motion solver.

SourceFiles
    motionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef motionSolver_H
#define motionSolver_H

#include "polyMesh.H"
#include "tetPolyMesh.H"
#include "tetPointFields.H"
#include "elementFields.H"
#include "tetFec.H"

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
private:

    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Motion mesh
        tetPolyMesh tetMesh_;

        //- Motion field
        tetPointVectorField motionU_;

        //- Total displacement field
        mutable tetPointVectorField* totDisplacementPtr_;

        //- 2-D motion corrector pointer
        twoDPointCorrector* correct2DPtr_;

public:

    //- Runtime type information
    TypeName("motionSolver");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            motionSolver,
            dictionary,
            (const polyMesh& mesh),
            (mesh)
        );


    // Selectors

        //- Select null constructed
        static autoPtr<motionSolver> New(const polyMesh& mesh);


    // Constructors

        //- Construct from components

        motionSolver(const polyMesh& mesh);


    // Destructor

        virtual ~motionSolver();


    // Member Functions
        
        //- Return reference to mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to motion mesh
        const tetPolyMesh& tetMesh() const
        {
            return tetMesh_;
        }

        //- Return reference to motion field
        tetPointVectorField& motionU()
        {
            return motionU_;
        }

        //- 
        void storeTotDisplacement() const
        {
            if(!totDisplacementPtr_)
            {
                totDisplacementPtr_ = 
                    new tetPointVectorField
                    (
                        IOobject
                        (
                            "totalMotionU",
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        tetMesh(),
                        dimensionedVector("zero", dimLength, vector::zero)
                    );
            }
        }

        bool needTotDisplacement() const
        {
            return totDisplacementPtr_;
        }


        //- Return reference to total displacement field
        tetPointVectorField& totDisplacement()
        {
            return *totDisplacementPtr_;
        }

        const tetPointVectorField& totDisplacement() const
        {
            return *totDisplacementPtr_;
        }

        //-
        const twoDPointCorrector* correct2DPtr() const
        {
            return correct2DPtr_;
        }

        //- Provide new points for motion.  Solves for motion
        virtual tmp<pointField> newPoints();

        //- Provide current points for motion.  Uses current motion field
        tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve() = 0;


        //- Update topology
        virtual void updateTopology();


        //- Distortion energy density
        tmp<elementScalarField> distortionEnergy() const;

        //- Deformation energy density
        tmp<elementScalarField> deformationEnergy() const;

        //- Total distortion energy density
        tmp<elementScalarField> totDistortionEnergy() const;

        //- Total deformation energy density
        tmp<elementScalarField> totDeformationEnergy() const;

        //- Return mean cell non-orthogonality
        tmp<scalarField> nonOrthogonality() const;

        //- Return mean cell skewness
        tmp<Foam::scalarField> cellSkewness() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
