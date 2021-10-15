// read the corresponding row of stoechiometric matrix
List<scalar> nuLst(transportProperties.subDict(currSpName).lookup("nuVec"));

// load the current specie
const volScalarField&  currSp = mesh_.lookupObject<volScalarField>(currSpName);

// load data for numerics form transportProperties
scalar  minCForGradND(readScalar(transportProperties.subDict(currSpName).lookup("minForGrad")));
dimensionedScalar minCForGrad("minCForGrad",currSp.dimensions(),minCForGradND);

volScalarField cForGrad(max(currSp,minCForGrad));


forAll(nuLst,rnI)                                                       //for all reactions
{    
    if (mag(nuLst[rnI]) > SMALL and activeReacts[rnI])                  //tango only in active reactions
    {
        Pout << currSpName << " is active in reaction " << rnI;
        Pout << " => loading kinetics parameters (k0 and EA)" << endl;
        
        char rName [50];
        label nChars = sprintf(rName,"reaction%02d",rnI);
        
        dimensionedScalar k0  = transportProperties.subDict("reactProps").subDict(rName).lookup("k0");
        dimensionedScalar EA  = transportProperties.subDict("reactProps").subDict(rName).lookup("EA");
        
        //counting exponent for total concetration
        //NOTE: I start with -1.0 because of units of kinetic constants 
        //      in transport properties (maybe there can be problem with
        //      non-integer alpha)
        float cTotExp = -1.0;
        forAll(specieNames,countI){
            //this can be probably done nicer, but I used the code above  
            List<scalar> alphaLstCurSp(transportProperties.subDict(specieNames[countI]).lookup("alphaVec"));    
            cTotExp += alphaLstCurSp[rnI];
        }
        
        // calculating real kinetic constant of reaction with arrhenius kinetics and scaling with concentration for right dimension
        dimensionedScalar kReact = k0*Foam::exp(-EA/(univR*cTemp))*Foam::pow(cTot,cTotExp);
        
        volScalarField auxR(0.0*R);
        volScalarField auxGradR(0.0*gradR);
    
        auxR        += nuLst[rnI]*kReact;
        auxGradR    += min(0.0,nuLst[rnI])*kReact;
                
        Pout << "reactants in reaction " << rnI << ": ";
        forAll(specieNames,countI)                                      //forAll other species
        {
            word testedSpName = specieNames[countI];
            List<scalar> nuLstTested(transportProperties.subDict(testedSpName).lookup("nuVec"));
            
            if (nuLstTested[rnI] < -SMALL)                              // if testedSpecie is reactant
            {
                Pout << testedSpName << " ";
                
                List<scalar> alphaLst(transportProperties.subDict(testedSpName).lookup("alphaVec"));
                const volScalarField& testedSp = mesh_.lookupObject<volScalarField>(testedSpName); 
            
                auxR *= Foam::pow(testedSp,alphaLst[rnI]);
                
                if (testedSpName != currSpName)
                {
                    auxGradR *= Foam::pow(testedSp,alphaLst[rnI]);
                }
                else
                {
                    auxGradR *= alphaLst[rnI]*Foam::pow(testedSp,alphaLst[rnI]-1.0);
                }
            }
        }
        Pout << endl;

        R     += auxR;
        gradR += auxGradR;
    }
    else
    {
        if (activeReacts[rnI])
        {
            Pout << currSpName << " is not active in reaction " << rnI << endl;
        }
        else
        {
            Pout << "reaction " << rnI << " switched-off" << endl;
        }
    }
}

// bounding the gradR field by maximal and minimal values            
dimensionedScalar minGradR("minGrad",gradR.dimensions(), clipGradMin);
dimensionedScalar maxGradR("maxGrad",gradR.dimensions(), clipGradMax);
gradR = min(max(gradR,minGradR),maxGradR);
            
eqn += isCoat*coatingVolFrac*fvm::Sp(gradR,currSp);
eqn += isCoat*coatingVolFrac*R - isCoat*coatingVolFrac*gradR*currSp;
// Note: Taylor expansion of kinetics
