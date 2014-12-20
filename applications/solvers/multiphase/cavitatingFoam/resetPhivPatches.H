fvPatchScalarFieldField& phivPatches = phiv.boundaryField();
const fvPatchVectorFieldField& Upatches = U.boundaryField();
const fvPatchVectorFieldField& SfPatches = mesh.Sf().boundaryField();

forAll(phivPatches, patchI)
{
    if (phiv.boundaryField().types()[patchI] == "calculated")
    {
        calculatedFvPatchScalarField& phivPatch =
            refCast<calculatedFvPatchScalarField>(phivPatches[patchI]);

        phivPatch == (Upatches[patchI] & SfPatches[patchI]);
    }
}
