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

Application
    yPlusLES

Description
    Calculates the yPlus of the near-wall cells for an LES.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/LESmodel/LESmodel.H"
#include "nearWallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

#       include "createFields.H"
        volScalarField nuEff = sgsModel->nuEff();

        volScalarField yPlus
        (
            IOobject
            (
                "yPlus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
        );

        const fvPatchList& patches = mesh.boundary();    
    
        forAll(patches, patchi)
        {
            const fvPatch& currPatch = patches[patchi];

            if (typeid(currPatch) == typeid(wallFvPatch))
            {
                yPlus.boundaryField()[patchi] =
                    d[patchi]
                   *sqrt
                    (
                        nuEff.boundaryField()[patchi]
                       *mag(U.boundaryField()[patchi].snGrad())
                    )
                   /sgsModel->nu().boundaryField()[patchi];
            }
        }

        yPlus.write();
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
