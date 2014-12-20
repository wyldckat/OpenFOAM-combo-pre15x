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

#include "PrimitivePatchInterpolation.H"
#include "faceList.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Patch>
const scalarField&
PrimitivePatchInterpolation<Patch>::faceToEdgeWeights() const
{
    if (!faceToEdgeWeightsPtr_)
    {
        makeFaceToEdgeWeights();
    }

    return *faceToEdgeWeightsPtr_;
}


template<class Patch>
void PrimitivePatchInterpolation<Patch>::makeFaceToEdgeWeights() const
{
    if (faceToEdgeWeightsPtr_)
    {
        FatalErrorIn
        (
            "PrimitivePatchInterpolation<Patch>::makeFaceToEdgeWeights() const"
        )   << "Face-to-edge weights already calculated"
            << abort(FatalError);
    }

    const pointField& points = patch_.localPoints();
    const faceList& faces = patch_.localFaces();
    const edgeList& edges = patch_.edges();
    const labelList& internalEdges = patch_.internalEdges();
    const labelListList& edgeFaces = patch_.edgeFaces();

    faceToEdgeWeightsPtr_ = new scalarField(internalEdges.size());
    scalarField& weights = *faceToEdgeWeightsPtr_;

    label curEdge, ownerFace, neighbourFace;

    forAll (internalEdges, edgeI)
    {
        curEdge = internalEdges[edgeI];

        ownerFace = edgeFaces[curEdge][0];
        neighbourFace = edgeFaces[curEdge][1];

        vector P = faces[ownerFace].centre(points);
        vector N = faces[neighbourFace].centre(points);
        vector S = points[edges[curEdge].start()];
        vector e = edges[curEdge].vec(points);

        scalar alpha = - ( ( (N - P)^(S - P) )&( (N - P)^e ) )/
            ( ( (N - P)^e )&( (N - P)^e ) );

        vector E = S + alpha*e;

        weights[edgeI] = mag(N - E)/(mag(N - E) + mag(E - P));
    }
}


template<class Patch>
void PrimitivePatchInterpolation<Patch>::clearWeights()
{
    deleteDemandDrivenData(faceToEdgeWeightsPtr_);
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class Patch>
PrimitivePatchInterpolation<Patch>::~PrimitivePatchInterpolation()
{
    clearWeights();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Patch>
template<class Type>
tmp<Field<Type> > PrimitivePatchInterpolation<Patch>::faceToPointInterpolate
(
    const Field<Type>& ff
) const
{
    // Check size of the given field
    if (ff.size() != patch_.size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > PrimitivePatchInterpolation::"
            "faceToPointInterpolate(const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << patch_.size() << " field size: " << ff.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            patch_.nPoints(), pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    // get reference to addressing
    const labelListList& pointFaces = patch_.pointFaces();

    forAll (pointFaces, pointI)
    {
        const labelList& curFaces = pointFaces[pointI];

        forAll (curFaces, faceI)
        {
            result[pointI] += ff[curFaces[faceI]];
        }

        result[pointI] /= curFaces.size();
    }

    return tresult;
}


template<class Patch>
template<class Type>
tmp<Field<Type> > PrimitivePatchInterpolation<Patch>::faceToPointInterpolate
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = faceToPointInterpolate(tff());
    tff.clear();
    return tint;
}


template<class Patch>
template<class Type>
tmp<Field<Type> > PrimitivePatchInterpolation<Patch>::pointToFaceInterpolate
(
    const Field<Type>& pf
) const
{
    if (pf.size() != patch_.nPoints())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > PrimitivePatchInterpolation::"
            "pointToFaceInterpolate(const Field<Type> pf)"
        )   << "given field does not correspond to patch. Patch size: "
            << patch_.nPoints() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            patch_.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    // get reference to addressing
    const faceList& localFaces = patch_.localFaces();

    forAll (result, faceI)
    {
        const labelList& curPoints = localFaces[faceI];

        forAll (curPoints, pointI)
        {
            result[faceI] += pf[curPoints[pointI]];
        }

        result[faceI] /= curPoints.size();
    }

    return tresult;
}


template<class Patch>
template<class Type>
tmp<Field<Type> > PrimitivePatchInterpolation<Patch>::pointToFaceInterpolate
(
    const tmp<Field<Type> >& tpf
) const
{
    tmp<Field<Type> > tint = pointToFaceInterpolate(tpf());
    tpf.clear();
    return tint;
}


template<class Patch>
template<class Type>
tmp<Field<Type> > PrimitivePatchInterpolation<Patch>::faceToEdgeInterpolate
(
    const Field<Type>& pf
) const
{
    // Check size of the given field
    if (pf.size() != patch_.size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > PrimitivePatchInterpolation::"
            "faceToEdgeInterpolate(const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << patch_.size() << " field size: " << pf.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            patch_.nEdges(), pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();
	
    const pointField& points = patch_.localPoints();
    const faceList& faces = patch_.localFaces();
    const edgeList& edges = patch_.edges();
    const labelList& internalEdges = patch_.internalEdges();
    const labelList& boundaryEdges = patch_.boundaryEdges();
    const labelListList& edgeFaces = patch_.edgeFaces();

    const scalarField& weights = faceToEdgeWeights();

    label curEdge, ownerFace, neighbourFace;

    forAll (internalEdges, edgeI)
    {
        curEdge = internalEdges[edgeI];

        ownerFace = edgeFaces[curEdge][0];
        neighbourFace = edgeFaces[curEdge][1];

        result[curEdge] = 
            weights[edgeI]*pf[ownerFace]
          + (1.0 - weights[edgeI])*pf[neighbourFace];
    }

    forAll (boundaryEdges, edgeI)
    {
        const label curEdge = boundaryEdges[edgeI];

        const label ownerFace = edgeFaces[curEdge][0];

        result[curEdge] = pf[ownerFace];
    }

    return tresult;
}


template<class Patch>
template<class Type>
tmp<Field<Type> > PrimitivePatchInterpolation<Patch>::faceToEdgeInterpolate
(
    const tmp<Field<Type> >& tpf
) const
{
    tmp<Field<Type> > tint = faceToEdgeInterpolate(tpf());
    tpf.clear();
    return tint;
}


template<class Patch>
bool PrimitivePatchInterpolation<Patch>::movePoints()
{
    clearWeights();

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
