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
    combineSampleSets

Description
    Various utility functions for collection of samples on processors.

SourceFiles
    combineSampleSets.C

\*---------------------------------------------------------------------------*/

#ifndef combineSampleSets_H
#define combineSampleSets_H

#include "labelList.H"
#include "ptrList.H"
#include "coordSet.H"
#include "sampleSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Combine points from all processors. Sort by curveDist and produce
// index list. Valid result only on master processor.
void combineSampleSets
(
    const Foam::ptrList<Foam::sampleSet>& sampleSets,
    Foam::ptrList<Foam::coordSet>& masterSampleSets,
    Foam::labelListList& indexSets
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
