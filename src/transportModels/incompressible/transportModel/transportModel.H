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
    transportModel

Description

SourceFiles
    transportModel.C
    newTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef transportModel_H
#define transportModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class transportModel Declaration
\*---------------------------------------------------------------------------*/

class transportModel
:
    public IOdictionary
{

protected:

    // Protected data

        const volVectorField& U_;
        const surfaceScalarField& phi_;


private:

        word phaseName_;


protected:

        dictionary phaseTransportProperties_;


    // Private Member Functions

        //- Disallow copy construct
        transportModel(const transportModel&);

        //- Disallow default bitwise assignment
        void operator=(const transportModel&);

        //- Lookup and return the phase transport properties dictionary
        static const dictionary& lookupPhaseTransportProperties
        (
            const dictionary& transportProperties,
            const word& phaseName = ""
        );


public:

    //- Runtime type information
    TypeName("transportModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            transportModel,
            dictionary,
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
                const word& phaseName
            ),
            (U, phi, phaseName)
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<transportModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& phaseName = ""
        );


    // Constructors

        //- Construct from components
        transportModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& phaseName = ""
        );


    // Destructor

        virtual ~transportModel()
        {}


    // Member Functions

        //- Return the name of this phase
        const word& phaseName() const
        {
            return phaseName_;
        }

        //- Return the phase transport properties dictionary
        const dictionary& phaseTransportProperties() const
        {
            return phaseTransportProperties_;
        }

        //- Return the strain rate mag(grad(U))
        tmp<volScalarField> strainRate() const;

        //- Return the laminar viscosity
        virtual const volScalarField& nu() const = 0;

        //- Correct the laminar viscosity
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
