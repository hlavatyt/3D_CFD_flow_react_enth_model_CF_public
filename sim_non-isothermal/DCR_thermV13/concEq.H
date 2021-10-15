// -- solve the concentration equation for documentation see Hlavaty 2021
forAll(specieNames,lstI)
{
    // -- get current name of the specie that is solved for
    word currSpName = specieNames[lstI];                            	

    Info<<"\nWorking on the specie: "<<currSpName<<endl;

    volScalarField& currSp	= MFrList[lstI];                	//current concentration field

    // -- DIFUSION COEFICIENT CALCULATION 
    // -- difussivity for mixture -- according to XPF see XPF doc. or 20_hlavaty.pdf

    // -- scalar constant --> add temperature dependence as volScalarField
    PtrList<dimensionedScalar> DABLst(specieNames.size()+1);
    forAll(DABLst,lstII)
    {
	    // not whole difusion coefficient -- see units, not dependent on temperature
	    DABLst.set
	    (
            lstII,
            new dimensionedScalar(
                "DAB",
                dimLength*dimLength/dimTime/Foam::pow(dimTemperature,1.75),
                ((0.014)/(pAtm*Foam::sqrt(2.0/(1e-3/molMList[lstII] + 1e-3/molMList[lstI]))*Foam::pow(Foam::pow(volDiffLst[lstI],0.3333333)+Foam::pow(volDiffLst[lstII],0.3333333),2))).value()
            )
	    );
    }
    
    //~ forAll(DABLst,lstII)
    //~ {
	//~ if(lstII < specieNames.size()){
	    //~ Info<<"D of " <<specieNames[lstI] << " with " <<specieNames[lstII]<< " is " << DABLst[lstII]<<endl;
	//~ }
	//~ else{
	    //~ Info<<"DRef is " << DABLst[lstII]<<endl;
	//~ }
    //~ }
    
    //Note: not a very nice thing inicializing denominator field to 0
    // volScalarField denominator(MFrList[lstI]/(DABLst[lstI]*Foam::pow(T,1.75)));
    // denominator -= MFrList[lstI]/(DABLst[lstI]*Foam::pow(T,1.75));

    // forAll(specieNames,lstII)
    // {
    //     if (lstII != lstI && min(MFrList[lstII]).value() > 0.0)
    //     {	
    //         denominator += MFrList[lstII] / (DABLst[lstII]*Foam::pow(T,1.75));
    //     }
    // }
    
    //Diffusion of the mixture -- default reference gas diffusion
    // adding temperature dependence
    volScalarField DifMix ((1-max(isOnWallCoat,isInWallCoat))*DABLst[DABLst.size()-1]*Foam::pow(T,1.75) + max(isOnWallCoat,isInWallCoat)*DABLst[DABLst.size()-1]*Foam::pow(T,1.75));
    // if (min(denominator).value() > 0){
//	volScalarField DifMix ((1.0 - MFrList[lstI])/denominator);
	//Info<<min(denominator).value()<<" "<<max(denominator).value()<<endl;
    // }

    // -- prepare Knudsen diffusivity based on XPF
    volScalarField DKn((1.0-max(isOnWallCoat,isInWallCoat))*(0.33333333333*2.0*dpC*sqrt((T*8*univR)/(3.141592653*molMList[lstI]))));
    DKn = (1.0-max(isOnWallCoat,isInWallCoat))*DifMix + (0.33333333333*sqrt((T*8*univR)/(3.141592653*molMList[lstI])))*isOnWallCoat*2.0*dpC + (0.33333333333*sqrt((T*8*univR)/(3.141592653*molMList[lstI])))*isInWallCoat*2.0*dpS;

    // -- prepare the effective diffusivity - Knudsen XPF
    volScalarField DEff(DifMix*(1.0-max(isOnWallCoat,isInWallCoat)));
    DEff  = DifMix*(1.0-max(isOnWallCoat,isInWallCoat));
    DEff += isOnWallCoat*(1.0/(tortC/porEpsC*(1.0/DifMix + 1.0/DKn)));
    DEff +=  isInWallCoat*(1.0/(tortS/porEpsS*(1.0/DifMix + 1.0/DKn)));
    DEff.correctBoundaryConditions();

    // --  calculation of total concentration
    // volScalarField cTot(p/univR/T);
    volScalarField cTot(thermo.rho()/molMList[molMList.size()-1]);
    // Info<<"max rhoTu1: "<<max(rhoTu).value()<< "min rhoTu: "<<min(rhoTu).value()<<endl;
    // Info<<"max rho: "<<max(thermo.rho()).value()<< "min rho: "<<min(thermo.rho()).value()<<endl;

    // -- inform about calculated properties
    Info<<"max T: "<<max(T).value()<<" min T: "<<min(T).value()<<endl;
    Info<<"max cTot: "<<max(cTot).value()<< "min cTot: "<<min(cTot).value()<<endl;
    Info<<"DEff "<<specieNames[lstI]<<" max: "<<max(DEff).value()<<" DEff "<<specieNames[lstI]<<" min: "<<min(DEff).value()<<endl;
    Info<<"ymax "<<specieNames[lstI]<<": "<<max(currSp).value()<<" ymin "<<specieNames[lstI]<<": "<<min(currSp).value()<<endl;
    // prepare the effective diffusivity - previous
    //~ volScalarField DEff(currD*(1.0-max(isOnWallCoat,isInWallCoat)));
    //~ DEff  = currD*(1.0-max(isOnWallCoat,isInWallCoat));
    //~ DEff += currD*isOnWallCoat*(coatingVolFracOW*porEpsC/constrC + (1.0-coatingVolFracOW));
    //~ DEff += currD*isInWallCoat*(coatingVolFracIW*porEpsC/constrC + (1.0-coatingVolFracIW)*(1.0-porEpsS-coatingVolFracIW)/constrS);
    //~ Info<<"cTemp: "<<cTemp.value()<<"DEff "<<specieNames[lstI]<<" max: "<<max(DEff).value()<<"DEff "<<specieNames[lstI]<<" min: "<<min(DEff).value()<<endl;

    // -- calculate molar flux at the face -- concEqn.H
    surfaceScalarField phiUCTot
    (
        IOobject
        (
            "phiU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(U*cTot) & mesh.Sf()
    );

    // -- balance of the specie equation
    fvScalarMatrix concEq
    (
        //fvm::ddt(currSp)                                               // fvm ... implicit
        fvm::div(phiUCTot, currSp)                                       // fvm ... implicit, convection term, phiUCTot initialized in createFields
        ==
        fvm::laplacian(DEff*cTot, currSp)                                // fvm ... implicit, first parameter = known, second parameter = unknown
        +fvc::laplacian(DEff*currSp, cTot)                               // fvm ... implicit, first parameter = known, second parameter = unknown
        + fvOptions(dummyConcEqn,currSp)                                       // source term defined in fvOptions, linearized -- part implicit, part explicit
    );
	    
    fvOptions.constrain(concEq);

    concEq.relax();

    concEq.solve();

    fvOptions.correct(currSp);
	    
    forAll(currSp,cellI)
    { 
    //	 currSp[cellI] = min(max(currSp[cellI],0),1);                    //clip unreal molar fractions
        currSp[cellI] = max(currSp[cellI],0); 
    }

    currSp.correctBoundaryConditions();
}

// forAll(specieNames,lstI)
// {//this should keep the system consistent even for the case of relaxation
//     volScalarField&     currSp      = MFrList[lstI];                //current concentration field
//     currSp.relax();                                                 //relax AFTER the first cycle
// }

// Note: update all -> relax
