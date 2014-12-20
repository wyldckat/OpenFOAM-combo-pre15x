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

Description

\*---------------------------------------------------------------------------*/

#include "ignitionSite.H"
#include "engineTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ignitionSite::ignitionSite(Istream& is, const Time& db, const fvMesh& mesh)
:
    db_(db),
    ignitionSiteDict_(is),
    location_(ignitionSiteDict_.lookup("location")),
    diameter_(readScalar(ignitionSiteDict_.lookup("diameter"))),
    time_
    (
        db_.userTimeToTime
        (
            readScalar(ignitionSiteDict_.lookup("start"))
        )
    ),
    duration_
    (
        db_.userTimeToTime
        (
            readScalar(ignitionSiteDict_.lookup("duration"))
        )
    ),
    strength_(readScalar(ignitionSiteDict_.lookup("strength")))
{
    // Check state of Istream
    is.check("ignitionSite::ignitionSite(Istream&)");

    findIgnitionCells(mesh);
}


ignitionSite::ignitionSite
(
    Istream& is,
    const engineTime& edb,
    const fvMesh& mesh
)
:
    db_(edb),
    ignitionSiteDict_(is),
    location_(ignitionSiteDict_.lookup("location")),
    diameter_(readScalar(ignitionSiteDict_.lookup("diameter"))),
    time_
    (
        db_.userTimeToTime
        (
            edb.degToTime(readScalar(ignitionSiteDict_.lookup("start")))
        )
    ),
    duration_
    (
        db_.userTimeToTime
        (
            edb.degToTime(readScalar(ignitionSiteDict_.lookup("duration")))
        )
    ),
    strength_(readScalar(ignitionSiteDict_.lookup("strength")))
{
    // Check state of Istream
    is.check("ignitionSite::ignitionSite(Istream&)");

    findIgnitionCells(mesh);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
