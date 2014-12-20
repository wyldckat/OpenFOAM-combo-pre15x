fvPatchScalarFieldField& phiPatches = phi.boundaryField();
const fvPatchScalarFieldField& rhoPatches = rho.boundaryField();
const fvPatchVectorFieldField& Upatches = U.boundaryField();
const fvPatchVectorFieldField& SfPatches = mesh.Sf().boundaryField();

forAll(phiPatches, patchI)
{
    if (phi.boundaryField().types()[patchI] == "calculated")
    {
        calculatedFvPatchScalarField& phiPatch =
            refCast<calculatedFvPatchScalarField>(phiPatches[patchI]);

        phiPatch == ((rhoPatches[patchI]*Upatches[patchI]) & SfPatches[patchI]);
    }
}
