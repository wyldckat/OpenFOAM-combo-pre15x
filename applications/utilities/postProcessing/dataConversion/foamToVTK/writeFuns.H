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
    writeFunctions

Description
    Various functions for collecting and writing binary data.

SourceFiles
    writeFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef writeFunctions_H
#define writeFunctions_H

#include "labelList.H"
#include "floatScalar.H"
#include "OFstream.H"
#include "DynamicList.H"
#include "point.H"
#include "volFields.H"
#include "pointFields.H"
#include "vtkMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class writeFuns Declaration
\*---------------------------------------------------------------------------*/

class writeFuns
{
    //- Swap halves of word.
    static void swapWord(label& word32);
    static void swapWords(const label nWords, label* words32);

public:

    //- Write ascii or binary. If binary optionally in-place swaps argument
    static void write(std::ostream&, const bool, DynamicList<floatScalar>&);
    static void write(std::ostream&, const bool, DynamicList<label>&);

    //- Write ascii or binary. If binary optionally in-place swaps argument
    static void write(std::ostream&, const bool, List<floatScalar>&);
    static void write(std::ostream&, const bool, labelList&);

    //- Append point to DynamicList
    static void insert(const point&, DynamicList<floatScalar>& dest);

    //- Append elements of labelList to DynamicList
    static void insert(const labelList&, DynamicList<label>&);
    //- Append elements of scalarList to DynamicList
    static void insert(const List<scalar>&, DynamicList<floatScalar>&);

    //- Append elements of scalarList to DynamicList using map
    static void insert
    (
        const labelList& map,
        const List<scalar>& source,
        DynamicList<floatScalar>&
    );

    //- Append points to DynamicList of floats
    static void insert(const List<point>& source, DynamicList<floatScalar>&);
    //- As above but using map
    static void insert
    (
        const labelList& map,
        const List<point>& source,
        DynamicList<floatScalar>&
    );

    //- Write volScalarField with cell values (including decomposed cells)
    static void writeVSF
    (
        std::ostream& os,
        const bool binary,
        const volScalarField& vsf,
        const vtkMesh& vMesh
    );
    static void writeVVF
    (
        std::ostream& os,
        const bool binary,
        const volVectorField& vvf,
        const vtkMesh& vMesh
    );

    //- Write pointScalarField on all mesh points. Interpolate to cell centre
    //  for decomposed cell centres.
    static void writePSF
    (
        std::ostream& os,
        const bool binary,
        const pointScalarField& psf,
        const vtkMesh& vMesh
    );
    static void writePVF
    (
        std::ostream& os,
        const bool binary,
        const pointVectorField& pvf,
        const vtkMesh& vMesh
    );

    //- Write interpolated field on points and original cell values on
    //  decomposed cell centres.
    static void writeIntVSF
    (
        std::ostream& os,
        const bool binary,
        const volScalarField& vsf,
        const pointScalarField& psf,
        const vtkMesh& vMesh
    );

    static void writeIntVVF
    (
        std::ostream& os,
        const bool binary,
        const volVectorField& vvf,
        const pointVectorField& pvf,
        const vtkMesh& vMesh
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
