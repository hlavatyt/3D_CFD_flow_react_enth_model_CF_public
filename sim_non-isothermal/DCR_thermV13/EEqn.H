{
    // -- load temperature from thermo package
    T = thermo.T();
    // -- calculate effective thermal diffusivity
    volScalarField kappaEff
    (
      (1.0-max(isOnWallCoat,isInWallCoat))*kappa
      + max(isOnWallCoat,isInWallCoat) * kappa
        //   (1.0-max(isOnWallCoat,isInWallCoat))*kappa
        // + isOnWallCoat*(coatingVolFracOW*((1.0-porEpsC.value())*kappaTC + porEpsC.value()*kappa) + (1.0-coatingVolFracOW)*kappa)
        // + isInWallCoat*(coatingVolFracIW*((1.0-porEpsC.value())*kappaTC + porEpsC.value()*kappa) + (1.0-coatingVolFracIW)*((1.0-porEpsS.value())*kappaTS + porEpsS.value()*kappa))
          //~ (1.0-max(isOnWallCoat,isInWallCoat))*kappa+
          //~ max(isOnWallCoat,isInWallCoat)*kappa
        //~ + isOnWallCoat*(coatingVolFracOW*((1.0-porEpsC.value())*kappaTC + porEpsC.value()*kappa) + (1.0-coatingVolFracOW)*kappa)
        //~ + isInWallCoat*(coatingVolFracIW*((1.0-porEpsC.value())*kappaTC + porEpsC.value()*kappa) + (1.0-coatingVolFracIW)*((1.0-porEpsS.value())*kappaTS + porEpsS.value()*kappa))
    );
    
    // -- calculate heat capacity of the gas
    volScalarField Cp
    (
        (ALst[ALst.size()-1] + BLst[BLst.size()-1] * T + CLst[CLst.size()-1] * Foam::pow(T,2) + DLst[DLst.size()-1] * Foam::pow(T,3)) / molMList[molMList.size()-1]

    );
      	
    // NOTE: discuss how to treat CpEff and rhoEff
    // -- calculate efective heat capacity
    volScalarField CpEff
    (
          Cp
        //    (1.0-max(isOnWallCoat,isInWallCoat))*Cp   
        //  + isOnWallCoat*(coatingVolFracOW*((1.0-porEpsC.value())*CpC + porEpsC.value()*Cp) + (1.0-coatingVolFracOW)*Cp)
        // //  + isOnWallCoat*(CpC)
        //  + isInWallCoat*(coatingVolFracIW*((1.0-porEpsC.value())*CpC + porEpsC.value()*Cp) + (1.0-coatingVolFracIW)*((1.0-porEpsS.value())*CpS + porEpsS.value()*Cp))
        //~ + isInWallCoat*(CpS)
    );

    //mass density -- using reference gas
    //~ volScalarField rho(molMList[molMList.size()-1]*pAtm/univR/T);
    //NOTE: this is inside UEqn.H
    
    //efective rho
    volScalarField rhoEff
    (
          thermo.rho()
          //  (1.0-max(isOnWallCoat,isInWallCoat))*thermo.rho()
        // (1.0-max(isOnWallCoat,isInWallCoat))*rhoCons+
          //~ max(isOnWallCoat,isInWallCoat)*rhoCons
        //  + isOnWallCoat*(coatingVolFracOW*((1.0-porEpsC.value())*rhoC + porEpsC.value()*thermo.rho()) + (1.0-coatingVolFracOW)*thermo.rho())
        //  + isInWallCoat*(coatingVolFracIW*((1.0-porEpsC.value())*rhoC + porEpsC.value()*thermo.rho()) + (1.0-coatingVolFracIW)*((1.0-porEpsS.value())*rhoS + porEpsS.value()*thermo.rho()))
        //~ + isOnWallCoat*(rhoC)
        //~ + isInWallCoat*(rhoS)
    );    
    
    
    // dissipation rate 
    volScalarField K (0.5*magSqr(U));
    
    // flux field at surface face without heat capacity 
    // surfaceScalarField phiRhoMy
    // (
    //     IOobject
    //     (
    //         "phiRho",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     linearInterpolate(U*rhoEff) & mesh.Sf()
    //     //linearInterpolate(U*thermo.rho()) & mesh.Sf()
    // );
    // -- know part of heat flux at the surface of the cell
    surfaceScalarField phiRhoCpMy
    (
        IOobject
        (
            "phiRhoCp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(U*rhoEff*CpEff) & mesh.Sf()
        //~ phi * linearInterpolate(CpEff) & mesh.Sf()
        //linearInterpolate(U*thermo.rho()*Cp) & mesh.Sf()
    );

    // write fields for corection
    Info << "kappa " << min(kappaEff).value()<<" "<<max(kappaEff).value() << endl;
    Info << "CpEff " << min(CpEff).value()<<" "<<max(CpEff).value() << endl;
    Info << "Cp " << min(Cp).value()<<" "<<max(Cp).value() << endl;
    Info << "rho " << min(thermo.rho()).value()<<" " <<max(thermo.rho()).value() << endl;
    Info << "rhoEff " << min(rhoEff).value()<<" " <<max(rhoEff).value() << endl;
    Info << "T "    << min(T) <<" "<< max(T) << endl;
    Info << "T thermo"    << min(thermo.T()) <<" "<< max(thermo.T()) << endl;
    
    // enthalpy balance
    fvScalarMatrix EEqn
    (
        fvm::div(phiRhoCpMy, T)
      //  + fvc::div(phiRhoMy, K)
      //~ //= fvm::laplacian(kappaEff, T)
        ==
        fvm::laplacian(kappaEff, T)
        + fvOptions(dummyEEqn,T)
        // fvm::div(phiRhoCpMy, thermo.T())
      // + fvc::div(phiRhoMy, K)
      //= fvm::laplacian(kappaEff, T)
      
      // - fvm::laplacian(kappaEff, thermo.T())
      // ==
        
      //  fvOptions(dummy2,thermo.T())
      //    fvm::div(phiRhoMy, thermo.he())
      //  + fvc::div(phiRhoMy, K)
      //~ //= fvm::laplacian(kappaEff, T)
      
      // - fvm::laplacian(kappaEff/CpEff, thermo.he())
      //  ==
        
      //  fvOptions(dummy2/CpEff,thermo.he())
    );
    
    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    // fvOptions.correct(thermo.T());
    fvOptions.correct(T);
    Info << "T "    << min(T) <<" "<< max(T) << endl;



    // thermo.T() = T;
    
    
    
    // T = thermo.he()/CpEff;
    // thermo.he() = thermo.he() * CpEff;
    // thermo.rho() = thermo.p();

    Info << "thermo.T "    << min(thermo.T()) <<" "<< max(thermo.T()) << endl;

    // fvScalarMatrix EEqnDva
    // (
    //     thermo.he()
    //    ==
    //     T * CpEff
    // );

    // EEqnDva.solve();


    thermo.he() = T * Cp;
    
    thermo.correct();
    
    Info << "thermo.T "    << min(thermo.T()) <<" "<< max(thermo.T()) << endl;

    
    //~ radiation->correct();
    // T = thermo.T();
    
}
