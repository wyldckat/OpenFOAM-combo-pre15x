if(turbulence)
{
    if (mesh.moving())
    {
        y.correct();
    }

    tmp<volTensorField> tgradUb = fvc::grad(Ub);
    volScalarField G = 2*nutb*(tgradUb() && dev(symm(tgradUb())));
    tgradUb.clear();

// Dissipation equation

#   include "wallDissipation.H"

    fvScalarMatrix epsEqn
    (
        fvm::ddt(beta, epsilon)
      + fvm::div(phib, epsilon)
      - fvm::laplacian(alphaEps*nuEffb, epsilon)
      ==
         C1*beta*G*epsilon/k + beta*boundarySource
       - fvm::Sp(C2*beta*epsilon/k + beta*boundaryCentral, epsilon)
    );
    epsEqn.relax();
    epsEqn.solve();

    epsilon.max(dimensionedScalar("zero", epsilon.dimensions(), 1.0e-15));


// Turbulent kinetic energy equation

    fvScalarMatrix kEqn
    (
        fvm::ddt(beta, k)
      + fvm::div(phib, k)
      - fvm::laplacian(alphak*nuEffb, k)
      ==
        beta*G
      - fvm::Sp(beta*epsilon/k, k)
    );
    kEqn.relax();
    kEqn.solve();

    k.max(dimensionedScalar("zero", k.dimensions(), 1.0e-8));

    //- Re-calculate turbulence viscosity
    nutb = Cmu*sqr(k)/epsilon;

#   include "wallViscosity.H"

}

nuEffa = sqr(Ct)*nutb + nua;
nuEffb = nutb + nub;
