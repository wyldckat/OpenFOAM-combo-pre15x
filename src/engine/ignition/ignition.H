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
    ignition

Description

SourceFiles
    ignition.C
    ignitionIO.C

\*---------------------------------------------------------------------------*/

#ifndef ignition_H
#define ignition_H

#include "ignitionSite.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class ignition Declaration
\*---------------------------------------------------------------------------*/

class ignition
{
    // Private data

        const fvMesh& mesh_;

        Switch ignite_;

        ptrList<ignitionSite> ignSites_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ignition(const ignition&);

        //- Disallow default bitwise assignment
        void operator=(const ignition&);


public:

    // Constructors

        //- Construct from Istream and database
        ignition(const dictionary&, const Time&, const fvMesh&);

        //- Construct from Istream and engineTime
        ignition(const dictionary&, const engineTime&, const fvMesh&);


    // Member Functions

        //- Return the set of ignition sites updated for moving meshes
        const ptrList<ignitionSite>& sites();


        // Check

            bool ignite() const
            {
                return ignite_;
            }

            //- Are any of the ignition site currently igniting
            bool igniting() const;

            //- Has the mixture been ignited?
            bool ignited() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
