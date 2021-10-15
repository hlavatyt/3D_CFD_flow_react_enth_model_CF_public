    forAll(specieNames,lstI)
    {
        word currSpName = specieNames[lstI];                            //get current specie name
        
        Info<<"Working on the species: "<<currSpName<<endl;

        //get the current variables
        //~dimensionedScalar   currD       = DList[lstI];                  //current diffusivity (free)
        
        volScalarField&     currSp      = MFrList[lstI];                //current concentration field
        
        //~volScalarField   currD        = A0List[lstI]*Foam::pow(T/dimensionedScalar("aux",dimTemperature,1.0),1.75)/pAtm*dimensionedScalar("aux",dimPressure,1.0);                  //current diffusivity (free)
         dimensionedScalar   currD        = A0List[lstI]*Foam::pow(cTemp/dimensionedScalar("aux",dimTemperature,1.0),1.75)/pAtm*dimensionedScalar("aux",dimPressure,1.0);                  //current diffusivity (free)
        
        // Note: not a very nice thing. A0 carries the diffusivity dimensions,
        //       so I make the temperature and pressure dimensionless
        
        //~Info<< "Dvol max/min : " << max(currD).value() << " " << min(currD).value()
        //~<< endl;
        
        //difussivity for mixture -- XPF
        //~volScalarField denominator(MFrList[lstI]/(A0List[lstI]*Foam::pow(cTemp/dimensionedScalar("aux",dimTemperature,1.0),1.75)/pAtm*dimensionedScalar("aux",dimPressure,1.0)));
        //~denominator -= MFrList[lstI]/(A0List[lstI]*Foam::pow(cTemp/dimensionedScalar("aux",dimTemperature,1.0),1.75)/pAtm*dimensionedScalar("aux",dimPressure,1.0));
		 //Info<<"denominatormin: "<<min(denominator).value()<<"denominatormax: "<<max(denominator).value()<<endl;
        
        //~volScalarField soucet(MFrList[lstI]);
        //~soucet -= MFrList[lstI];
        //~forAll(specieNames,lstI2)
        //~{
			//~if (min(MFrList[lstI2]).value() > 1.0e-11)
			//~{
				//~denominator += MFrList[lstI2] /  (A0List[lstI2]*Foam::pow(cTemp/dimensionedScalar("aux",dimTemperature,1.0),1.75)/pAtm*dimensionedScalar("aux",dimPressure,1.0));
                //~soucet += MFrList[lstI2];
            //~}
		//~}
        //~denominator += soucet /  (A0List[0]*Foam::pow(cTemp/dimensionedScalar("aux",dimTemperature,1.0),1.75)/pAtm*dimensionedScalar("aux",dimPressure,1.0));
		
		//~Info<<"denominatormin: "<<min(denominator).value()<<endl;
		//~volScalarField   currD ((1.0-max(isOnWallCoat,isInWallCoat))*A0List[lstI]*Foam::pow(cTemp/dimensionedScalar("aux",dimTemperature,1.0),1.75)/pAtm*dimensionedScalar("aux",dimPressure,1.0)+max(isOnWallCoat,isInWallCoat)*A0List[lstI]*Foam::pow(cTemp/dimensionedScalar("aux",dimTemperature,1.0),1.75)/pAtm*dimensionedScalar("aux",dimPressure,1.0));
		//~Info<<"currDMIND: "<<min(currD).value()<<endl;
		//~if (min(denominator).value() != 0.0)
		//~{
				//~currD = (1.0 - MFrList[lstI]) / denominator;
				//~Info<<"In if: currDMIND: "<<min(currD).value()<<endl;
        //~}
        //~Info<<"afterIFcurrDMIND: "<<min(currD).value()<<endl;
        //prepare Knudsen diffusivity based on XPF
        volScalarField DKn((1.0-max(isOnWallCoat,isInWallCoat))*(0.33333333333*2.0*dpC*sqrt((cTemp*8*univR)/(3.141592653*molMList[lstI]))));
		DKn = (1.0-max(isOnWallCoat,isInWallCoat))*currD + (0.33333333333*sqrt((cTemp*8*univR)/(3.141592653*molMList[lstI])))*isOnWallCoat*2.0*dpC +  (0.33333333333*sqrt((cTemp*8*univR)/(3.141592653*molMList[lstI])))*isInWallCoat*2.0*dpS;
        
        Info<<"cTemp: "<<cTemp.value()<<"M "<<specieNames[lstI]<<" = "<<molMList[lstI].value()<<"DKn "<<specieNames[lstI]<<" max: "<<max(DKn).value()<<"DEff "<<specieNames[lstI]<<" min: "<<min(DKn).value()<<endl;
        //prepare the effective diffusivity -Knudsen XPF
        volScalarField DEff(currD*(1.0-max(isOnWallCoat,isInWallCoat)));
        DEff  = currD*(1.0-max(isOnWallCoat,isInWallCoat));
        DEff += isOnWallCoat*(1.0/(tortC/porEpsC*(1.0/currD + 1.0/DKn)));
        DEff +=  isInWallCoat*(1.0/(tortS/porEpsS*(1.0/currD + 1.0/DKn)));
        
        // prepare the effective diffusivity - previous
        //~ volScalarField DEff(currD*(1.0-max(isOnWallCoat,isInWallCoat)));
        //~ DEff  = currD*(1.0-max(isOnWallCoat,isInWallCoat));
        //~ DEff += currD*isOnWallCoat*(coatingVolFracOW*porEpsC/constrC + (1.0-coatingVolFracOW));
        //~ DEff += currD*isInWallCoat*(coatingVolFracIW*porEpsC/constrC + (1.0-coatingVolFracIW)*(1.0-porEpsS-coatingVolFracIW)/constrS);
        DEff.correctBoundaryConditions();
        
        Info<<"cTemp: "<<cTemp.value()<<"DEff "<<specieNames[lstI]<<" max: "<<max(DEff).value()<<"DEff "<<specieNames[lstI]<<" min: "<<min(DEff).value()<<endl;
        
        fvScalarMatrix concEq
        (
            //fvm::ddt(currSp)                                            // fvm ... implicit
           fvm::div(phi, currSp)                                       // fvm ... implicit, convection term, phi has to be initialized in createFields.H
         ==
            //~ fvm::laplacian(effD, currSp)                                // fvm ... implicit, first parameter = known, second parameter = unknown
            fvm::laplacian(DEff, currSp)                                // fvm ... implicit, first parameter = known, second parameter = unknown
          + fvOptions(currSp)
        );
                
        fvOptions.constrain(concEq);
        
        concEq.relax();
        
        concEq.solve();

        fvOptions.correct(currSp);
                
        //~ forAll(currSp,cellI)
        //~ {
            //~ currSp[cellI] = min(max(currSp[cellI],0),1);                //clip unreal molar fractions
        //~ }

        currSp.correctBoundaryConditions();
    }
    
    forAll(specieNames,lstI)
    {//this should keep the system consistent even for the case of relaxation
        volScalarField&     currSp      = MFrList[lstI];                //current concentration field
        currSp.relax();                                                 //relax AFTER the first cycle
    }
    
    // Note: update all -> relax
