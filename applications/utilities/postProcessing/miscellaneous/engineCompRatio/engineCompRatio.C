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
    engineCompRatioFoam

Description
    Calculate the geometric compression ratio.
    Note that if you have valves and/or extra volumes it will not work,
    since it calculates the volume at BDC and TCD.

\*---------------------------------------------------------------------------*/

#include "engineTime.H"
#include "fvCFD.H"
#include "tetFem.H"
#include "fixedValueTetPolyPatchFields.H"
#include "slipTetPolyPatchFields.H"
#include "symmetryFvPatch.H"
#include "wedgeFvPatch.H"
#include "emptyFvPatch.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createEngineTime.H"
#   include "createMeshNoClear.H"
#   include "createEngineMovingMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalar eps = 1.0e-10;
    scalar fullCycle = 360.0;

    scalar ca0 = -180.0;
    scalar ca1 = 0.0;
    scalar ca = runTime.theta();

    while (ca > ca0)
    {
        ca0 += fullCycle;
        ca1 += fullCycle;
    }

    if (mag(ca - ca0) > eps)
    {
        while(mag(ca - ca0) > eps)
        {
            ca = runTime.theta();
            scalar t0 = runTime.userTimeToTime(ca0 - ca);
            runTime.setDeltaT(t0);
            runTime++;
            Info << "CA = " << ca << endl;
#           include "movePiston.H"
        }
    }

    scalar Vmax = sum(mesh.V());

    while(mag(ca-ca1) > eps)
    {
        ca = runTime.theta();
        scalar t1 = runTime.userTimeToTime(ca1-ca);
        runTime.setDeltaT(t1);
        runTime++;
        Info << "CA = " << runTime.theta() << endl;
#       include "movePiston.H"
    }

    scalar Vmin = sum(mesh.V());

    Info << "\nVmax = " << Vmax;
    Info << ", Vmin = " << Vmin << endl;
    Info << "Vmax/Vmin = " << Vmax/Vmin << endl;
    Info << "\n end\n";

    return(0);
}


// ************************************************************************* //
