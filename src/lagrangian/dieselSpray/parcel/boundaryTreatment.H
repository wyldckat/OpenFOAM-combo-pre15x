if (isType<wallPolyPatch>(pbMesh[patch(facei)]))
{
    keepParcel = sDB.wall().wallTreatment(*this, facei);

    if (sDB.twoD())
    {
        scalar vs = U() & sDB.axisOfSymmetry();
        vector Us = U() - (vs * sDB.axisOfSymmetry());
        scalar magUs = mag(Us);
        scalar v1 = Us & n();
        scalar magV1 = mag(v1);
        U() = (magUs/magV1)*v1*n() + vs*sDB.axisOfSymmetry();
    }
}
else if (isType<wedgePolyPatch>(pbMesh[patch(facei)]))
{
    // check if parcel is trying to move out of the domain
    label patchi = patch(facei);
    label patchFacei = patchFace(patchi, facei);
    const polyPatch& patch = mesh.boundaryMesh()[patchi];
    vector nf = patch.faceAreas()[patchFacei];

    scalar Un = U() & nf;
    if (Un > 0)
    {
        scalar Un2 = U() & n(); 
        U() -= 2.0*Un2*n();
    }
}
else if (isType<symmetryPolyPatch>(pbMesh[patch(facei)]))
{
    // check if parcel is trying to move out of the domain
    label patchi = patch(facei);
    label patchFacei = patchFace(patchi, facei);
    const polyPatch& patch = mesh.boundaryMesh()[patchi];
    vector nf = patch.faceAreas()[patchFacei];

    scalar Un = U() & nf;
    if (Un > 0)
    {
        if (sDB.twoD())
        {
            U() -= 2.0*(U() & n())*n();
        }
        else
        {
            nf /= mag(nf);
            U() -= 2.0*(U() & nf)*nf;
        }
    }
}
else if (isType<polyPatch>(pbMesh[patch(facei)]))
{
    // Parcel has hit an inlet or outlet or some such so throw it away
    keepParcel = false;
}
