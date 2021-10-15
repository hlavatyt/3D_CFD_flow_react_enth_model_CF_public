    forAll(specieNames,lstI)
    {
        word currSpName = specieNames[lstI];                            //get current specie name
        
        Info<<"Working on the species: "<<currSpName<<endl;
        //get the current variables
        volScalarField&     currSp      = MFrList[lstI];                //current concentration field
        dimensionedScalar   currD       = DList[lstI];                  //current diffusivity (free)
        
        // prepare the effective diffusivity
        //~ volScalarField      effD(currD*isCoat*porEps*constr/tort);      //effective diffusivity in coating
        //~ volScalarField      effD(currD*porEps*constr/tort);      //effective diffusivity in coating
        //~ forAll(effD,cellI)
        //~ {
            //~ if(effD[cellI] == 0)                                        //0 -> original diffusivity
            //~ {
                //~ effD[cellI] = currD.value();
            //~ }
        //~ }
        //~ effD.correctBoundaryConditions();
        
        fvScalarMatrix concEq
        (
            fvm::ddt(currSp)                                            // fvm ... implicit
          + fvm::div(phi, currSp)                                       // fvm ... implicit, convection term, phi has to be initialized in createFields.H
         ==
            //~ fvm::laplacian(effD, currSp)                                // fvm ... implicit, first parameter = known, second parameter = unknown
            fvm::laplacian(currD, currSp)                                // fvm ... implicit, first parameter = known, second parameter = unknown
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
