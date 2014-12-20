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

#include "particle.H"
#include "Cloud.H"
#include "processorPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wedgePolyPatch.H"
#include "cyclicPolyPatch.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class particleType>
labelList particle<particleType>::findFaces
(
    const vector& position
)
{
    const polyMesh& mesh = cloud_.polyMesh_;
    const labelList& faces = mesh.cells()[celli_];
    const vector& C = mesh.cellCentres()[celli_];

    labelList faceList(0);
    forAll(faces, i)
    {
        label facei = faces[i];
        scalar lam = lambda(C, position, facei);

        if ((lam > 0) && (lam < 1.0))
        {
            label n = faceList.size();
            faceList.setSize(n+1);
            faceList[n] = facei;
        }
    }
    
    return faceList;
}


template<class particleType>
labelList particle<particleType>::findFaces
(
    const vector& position,
    const label celli,
    const scalar fraction
)
{
    const polyMesh& mesh = cloud_.polyMesh_;
    const labelList& faces = mesh.cells()[celli];
    const vector& C = mesh.cellCentres()[celli];

    labelList faceList(0);
    forAll(faces, i)
    {
        label facei = faces[i];
        scalar lam = lambda(C, position, facei, fraction);

        if ((lam > 0) && (lam < 1.0))
        {
            label n = faceList.size();
            faceList.setSize(n+1);
            faceList[n] = facei;
        }
    }
    
    return faceList;
}


template<class particleType>
void particle<particleType>::prepareForParallelTransfer
(
    const label patchi,
    const label facei
)
{
    celli_ = patchFace(patchi, facei);
}


template<class particleType>
void particle<particleType>::correctAfterParallelTransfer
(
    const label patchi,
    const label facei
)
{
    const processorPolyPatch& ppp = 
        refCast<const processorPolyPatch>
        (cloud_.pMesh().boundaryMesh()[patchi]);

    celli_ = ppp.faceCells()[facei];

    if (!ppp.parallel())
    {
        if (ppp.forwardT().size() == 1)
        {
            const tensor& T = ppp.forwardT()[0];
            transformPosition(T);
            transformProperties(T);
        }
        else
        {
            const tensor& T = ppp.forwardT()[facei];
            transformPosition(T);
            transformProperties(T);
        }
    }
    else if (ppp.separated())
    {
        if (ppp.separation().size() == 1)
        {
            position_ += ppp.separation()[0];
        }
        else
        {
            position_ += ppp.separation()[facei];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class particleType>
particle<particleType>::particle
(
    const Cloud<particleType>& cloud,
    const vector& position,
    const label celli
)
:
    cloud_(cloud),
    position_(position),
    celli_(celli),
    onBoundary_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class particleType>
label particle<particleType>::track
(
    const vector& endPosition,
    scalar& fraction
)
{
    scalar fractionLeft = 1.0;
    fraction = 0.0;
    scalar f = 0.0;
    bool onB = false;
    label facei = -1;
    // using the onBoundary() function in the while-statement
    // will not work for a particle traveling from boundary and away.
    while(!onB && (fractionLeft > SMALL))
    {
        f = fraction;
        facei = trackToFace(endPosition, f);
        fraction += fractionLeft*f;
        fractionLeft = 1.0 - fraction;
        if (facei > -1)
        {
            if (!cloud_.softAlgorithm_)
            {
                onB = onBoundary();
            }
        }
    }

    return facei;
}


template<class particleType>
label particle<particleType>::trackToFace
(
    const vector& endPosition,
    scalar& fraction
)
{
    const polyMesh& mesh = cloud_.polyMesh_;

    labelList faces = findFaces(endPosition);

    scalar fracOut = 0.0;
    label facei = -1;

    if (faces.size() == 0) // inside cell
    {
        position_ = endPosition;
        fracOut = 1.0;
        onBoundary_ = false;
    }
    else // hit face
    {
        scalar lambdaMin = GREAT;

        if (faces.size() == 1)
        {
            lambdaMin = lambda(position_, endPosition, faces[0], fraction);
            facei = faces[0];
        }
        else
        {
            // if the particle have to cross more than one cell to reach the
            // endPosition, we check which way to go.
            // if one of the faces is a boundary face and the particle is
            // outside we choose the boundary face.
            // the particle is outside if one of the lambda's is > 1 or < 0
            forAll(faces, i)
            {
                scalar lam = lambda(position_, endPosition, faces[i], fraction);

                if (lam < lambdaMin)
                {
                    lambdaMin = lam;
                    facei = faces[i];
                }
            }
        }

        bool internalFace = cloud_.internalFace(facei);

        // NN.
        // for warped faces the particle can be 'outside' the cell
        // this will yield a lambda larger than 1, or smaller than 0
        // for values < 0, the particle travels away from the cell
        // and we don't move the particle, only change cell.
        // for values larger than 1, we move the particle to endPosition only.
        if (lambdaMin > 0.0)
        {
            if (lambdaMin <= 1.0)
            {
                position_ += lambdaMin*(endPosition - position_);
                fracOut = lambdaMin;
            }
            else
            {
                fracOut = 1.0;
                position_ = endPosition;
            }
        }
        else
        {
            // soft-sphere particles can travel outside the domain
            // but we don't use lambda since this the particle
            // is going away from face
            if (cloud_.softAlgorithm_)
            {
                fracOut = 1.0;
                position_ = endPosition;
            }
        }

        // change cell
        if (internalFace) // Internal face
        {
            onBoundary_ = false;

            if (celli_ == cloud_.owner_[facei])
            {
                celli_ = cloud_.neighbour_[facei];
            }
            else if (celli_ == cloud_.neighbour_[facei])
            {
                celli_ = cloud_.owner_[facei];
            }
            else
            {
                FatalErrorIn
                (
                    "particle::trackToFace"
                    "(const vector& endPosition, scalar& fraction)"
                )   << "addressing failure"
                    << abort(FatalError);
            }
        }
        else
        {
            // On boundary, can't go any further
            onBoundary_ = true;

            // soft-sphere algorithm ignores the boundary
            if (cloud_.softAlgorithm_)
            {
                position_ = endPosition;
                fracOut = 1.0;
            }

            label patchi = patch(facei);
            const polyPatch& patch = mesh.boundaryMesh()[patchi];

            if (typeid(patch) == typeid(wedgePolyPatch))
            {
                //transformProperties(I - 2.0*n_*n_);
            }
            else if (typeid(patch) == typeid(symmetryPolyPatch))
            {
                /*
                label facei = patchFace(patchi);
                
                vector nf = patch.faceAreas()[facei];
                nf /= mag(nf);

                transformProperties(I - 2.0*nf*nf);
                */
            }
            else if (typeid(patch) == typeid(cyclicPolyPatch))
            {
                const cyclicPolyPatch& cpp =
                    refCast<const cyclicPolyPatch>(patch);

                label patchFacei = patchFace(patchi, facei);

                facei = cpp.transformGlobalFace(facei);

                celli_ = mesh.faceOwner()[facei];

                if (!cpp.parallel())
                {
                    const tensor& T = cpp.transformT(patchFacei);

                    transformPosition(T);
                    transformProperties(T);
                }
                else if (cpp.separated())
                {
                    position_ += cpp.separation(patchFacei);
                }
            }
            else if (typeid(patch) == typeid(processorPolyPatch))
            {
                // particle has hit a processor boundary
                celli_ = facei;
            }
        }        
    }

    fraction = fracOut;
    return facei;
}


template<class particleType>
scalar particle<particleType>::wallImpactDistance(const vector&) const
{
    return 0.0;
}


template<class particleType>
void particle<particleType>::transformPosition(const tensor& T)
{
    position_ = transform(T, position_);
}


template<class particleType>
void particle<particleType>::transformProperties(const tensor&)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "particleIO.C"

// ************************************************************************* //
