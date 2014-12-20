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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 information
inline bool Foam::regionInfo::update
(
    const regionInfo& w2,
    const label thisFaceI,
    const label thisCellI,
    const label neighbourFaceI,
    const label neighbourCellI
)
{
    if (w2.region() == -1)
    {
        FatalErrorIn
        (
            "regionInfo::update(const regionInfo&, const label,"
            "const label, const label)"
        ) << "Problem: w2 is -1" << abort(FatalError);

        return false;
    }
    else if (region() != -1)
    {
        // our already set. Stop any transfer
        return false;
    }
    else
    {
        region_ = w2.region();

        return true;
    }
}
    

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::regionInfo::regionInfo()
:
    region_(-1)
{}


// Construct from components
inline Foam::regionInfo::regionInfo(const label region)
:
    region_(region)
{}


// Construct as copy
inline Foam::regionInfo::regionInfo(const regionInfo& w2)
:
    region_(w2.region())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::regionInfo::valid() const
{
    return region_ != -1;
}


// No geometric data so never any problem on cyclics
inline bool Foam::regionInfo::sameGeometry
(
    const regionInfo&,
    const scalar
) const
{
    return true;
}


// No geometric data.
inline void Foam::regionInfo::leaveDomain
(
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{}


// No geometric data.
inline void Foam::regionInfo::transform(const tensor& rotTensor)
{}


// No geometric data.
inline void Foam::regionInfo::enterDomain
(
    const polyPatch& patch,
    const label patchFaceI,
    const point& faceCentre
)
{}


// Update this with neighbour information
inline bool Foam::regionInfo::updateCell
(
    const polyMesh& mesh,
    const label thisCellI,
    const label neighbourFaceI,
    const regionInfo& neighbourInfo,
    const scalar tol
)
{
    return update
    (
        neighbourInfo,
        -1,
        thisCellI,
        neighbourFaceI,
        -1
    );
}    


// Update this with neighbour information
inline bool Foam::regionInfo::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label neighbourCellI,
    const regionInfo& neighbourInfo,
    const scalar tol
)
{
    return update
    (
        neighbourInfo,
        thisFaceI,
        -1,
        -1,
        neighbourCellI
    );
}    

// Update this with neighbour information
inline bool Foam::regionInfo::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const regionInfo& neighbourInfo,
    const scalar tol
)
{
    return update
    (
        neighbourInfo,
        thisFaceI,
        -1,
        -1,
        -1
    );
}    


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::regionInfo::operator==(const Foam::regionInfo& rhs) const
{
    return region() == rhs.region();
}


inline bool Foam::regionInfo::operator!=(const Foam::regionInfo& rhs) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
