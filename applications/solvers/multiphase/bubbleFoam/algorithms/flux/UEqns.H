{
    volTensorField Rca = -nuEffa*(fvc::grad(Ua)().T());
    Rca += (2.0/3.0)*sqr(Ct)*I*k - (2.0/3.0)*I*tr(Rca);

    surfaceScalarField phiRa =
      - fvc::interpolate(nuEffa)
       *mesh.magSf()*fvc::snGrad(alpha)/fvc::interpolate(alpha + 0.001);

    UaEqn =
    (
        (1.0 + Cvm*rhob*beta/rhoa)*
        (
            fvm::ddt(Ua)
          + fvm::div(phia, Ua, "div(phia,Ua)")
          - fvm::Sp(fvc::div(phia), Ua)
        )

      - fvm::laplacian(nuEffa, Ua)
      + fvc::div(Rca)

      + fvm::div(phiRa, Ua, "div(phia,Ua)")
      - fvm::Sp(fvc::div(phiRa), Ua)
      + (fvc::grad(alpha)/(fvc::average(alpha) + 0.001) & Rca)
     ==
    //  g                        // Buoyancy term transfered to p-equation
      - fvm::Sp(beta/rhoa*dragCoef, Ua)
    //+ beta/rhoa*dragCoef*Ub    // Explicit drag transfered to p-equation
      - beta/rhoa*(liftCoeff - Cvm*rhob*DDtUb)
    );

    UaEqn.relax();


    volTensorField Rcb = -nuEffb*fvc::grad(Ub)().T();
    Rcb += (2.0/3.0)*I*k - (2.0/3.0)*I*tr(Rcb);

    surfaceScalarField phiRb =
       - fvc::interpolate(nuEffb)
        *mesh.magSf()*fvc::snGrad(beta)/fvc::interpolate(beta + 0.001);

    UbEqn = 
    (
        (1.0 + Cvm*rhob*alpha/rhob)*
        (
            fvm::ddt(Ub)
          + fvm::div(phib, Ub, "div(phib,Ub)")
          - fvm::Sp(fvc::div(phib), Ub)
        )

      - fvm::laplacian(nuEffb, Ub)
      + fvc::div(Rcb)

      + fvm::div(phiRb, Ub, "div(phib,Ub)")
      - fvm::Sp(fvc::div(phiRb), Ub)

      + (fvc::grad(beta)/(fvc::average(beta) + 0.001) & Rcb)
     ==
    //  g                        // Buoyancy term transfered to p-equation
      - fvm::Sp(alpha/rhob*dragCoef, Ub)
    //+ alpha/rhob*dragCoef*Ua   // Explicit drag transfered to p-equation
      + alpha/rhob*(liftCoeff + Cvm*rhob*DDtUa)
    );

    UbEqn.relax();
}
