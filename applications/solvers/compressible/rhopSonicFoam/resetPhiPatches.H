fvPatchScalarFieldField& phiPatches = phi.boundaryField();
const fvPatchVectorFieldField& rhoUpatches = rhoU.boundaryField();
const fvPatchVectorFieldField& SfPatches = mesh.Sf().boundaryField();

forAll(phiPatches, patchI)
{
    if (phi.boundaryField().types()[patchI] == "calculated")
    {
        calculatedFvPatchScalarField& phiPatch =
            refCast<calculatedFvPatchScalarField>(phiPatches[patchI]);

        phiPatch == (rhoUpatches[patchI] & SfPatches[patchI]);
    }
}
