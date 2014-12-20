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
    processorFaPatch

Description
    Processor patch.

SourceFiles
    processorFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorFaPatch_H
#define processorFaPatch_H

#include "coupledFaPatch.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class processorFaPatch Declaration
\*---------------------------------------------------------------------------*/

class processorFaPatch
:
    public coupledFaPatch
{
    // Private data

        int myProcNo_;
        int neighbProcNo_;

        //- Processor-neighbbour patch edge centres
        vectorField neighbEdgeCentres_;

        //- Processor-neighbbour patch edge lengths
        vectorField neighbEdgeLengths_;

        //- Processor-neighbbour patch neighbour face centres
        vectorField neighbEdgeFaceCentres_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

        //- Make patch face - neighbour cell distances
        void makeDeltaCoeffs(scalarField&) const;


protected:

    // Protected Member functions

        //- Initialise the calculation of the patch geometry
        void initGeometry();

        //- Calculate the patch geometry
        void calcGeometry();

        //- Initialise the patches for moving points
        void initMovePoints(const pointField&);

        //- Correct patches after moving points
        void movePoints(const pointField&);


public:

    //- Runtime type information
    TypeName(processorPolyPatch::typeName_());


    // Constructors

        //- Construct from dictionary
        processorFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        )
        :
            coupledFaPatch(name, dict, index, bm),
            myProcNo_(readLabel(dict.lookup("myProcNo"))),
            neighbProcNo_(readLabel(dict.lookup("neighbProcNo")))

        {}


    // Member functions

        //- Return processor number
        int myProcNo() const
        {
            return myProcNo_;
        }

        //- Return neigbour processor number
        int neighbProcNo() const
        {
            return neighbProcNo_;
        }

        //- Return true if running parallel
        virtual bool coupled() const
        {
            if (Pstream::parRun())
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        //- Return processor-neighbbour patch edge centres
        const vectorField& neighbEdgeCentres() const
        {
            return neighbEdgeCentres_;
        }

        //- Return processor-neighbbour patch edge lengths
        const vectorField& neighbEdgeLengths() const
        {
            return neighbEdgeLengths_;
        }

        //- Return processor-neighbbour patch neighbour face centres
        const vectorField& neighbEdgeFaceCentres() const
        {
            return neighbEdgeFaceCentres_;
        }


        //- Return delta (P to N) vectors across coupled patch
        tmp<vectorField> delta() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
