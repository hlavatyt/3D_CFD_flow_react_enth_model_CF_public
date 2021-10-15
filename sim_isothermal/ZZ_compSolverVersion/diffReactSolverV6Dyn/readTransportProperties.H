    Info<< "Reading transport properties\n" << endl;
    //Note: only for 2 materials (not at all nice)   
    
    //get the list of materials
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            mesh,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    );
    
    //prepare the lists for materials
    wordList specieNames
    (
        transportProperties.lookup("species")
    );
    
    //prepare the output variables - not very nice but hopefully working
    PtrList<dimensionedScalar> A0List(specieNames.size());
     PtrList<dimensionedScalar> molMList(specieNames.size());
    
    //prepare the output variables - not very nice but hopefully working
    //~PtrList<dimensionedScalar> DList(specieNames.size());
        
    //read the material properties for each of them
    forAll(specieNames,countI)
    {
        word currentSp = specieNames[countI];
        Info<< "Processing the specie: " << currentSp << endl;
        
        Info<< "    Preparing constant A0 for diffusivity computations\n" << endl;
        scalar sigmaV(readScalar(transportProperties.subDict(currentSp).lookup("sigmaV")));
        scalar sigmaVR(readScalar(transportProperties.subDict(currentSp).lookup("sigmaVR")));
        scalar molM(readScalar(transportProperties.subDict(currentSp).lookup("molM")));
        scalar molMR(readScalar(transportProperties.subDict(currentSp).lookup("molMR")));
        scalar sqrtMAB(Foam::sqrt(2.0/(1e-3/molM + 1e-3/molMR)));
        scalar sigmaA0(Foam::sqr(Foam::pow(sigmaV,0.33333333333333)+Foam::pow(sigmaVR,0.333333333)));
        
        
        A0List.set
        (
            countI,
            new dimensionedScalar("A0",dimLength*dimLength/dimTime,0.01449/(sqrtMAB*sigmaA0))
        );
        
        molMList.set
        (
            countI,
            new dimensionedScalar("molM", dimMass /dimMoles,molM)
           );
        
        // Note: this is rather ugly. I use empirical correlations, so my
        //       dimensions do not quite match. Consequently, I just prepare
        //       the constant A0 such that it carries the dimensions of diffusivity
    }
    
    // load the coating locations inwall
    const dictionary& isOnWallCoatDict(transportProperties.subDict("isOnWallCoat"));
    word isOnWallCoatType(isOnWallCoatDict.lookup("type"));

    autoPtr<volScalarField> isOnWallCoatPtr;

    IOobject isOnWallCoatIO
    (
        "isOnWallCoat",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (isOnWallCoatType == "uniform")
    {
        scalar isOnWallCoatValue(readScalar(isOnWallCoatDict.lookup("value")));

        isOnWallCoatPtr.reset
        (
            new volScalarField
            (
                isOnWallCoatIO,
                mesh,
                dimensionedScalar
                (
                    "isOnWallCoat",
                    dimless/dimTime,
                    isOnWallCoatValue
                ),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );
    }
    else if (isOnWallCoatType == "field")
    {
        isOnWallCoatIO.readOpt() = IOobject::MUST_READ;

        isOnWallCoatPtr.reset
        (
            new volScalarField
            (
                isOnWallCoatIO,
                mesh
            )
        );
    }
    else
    {
        FatalErrorIn
        (
            "readTransportProperties.H"
        )   << "Valid type entries are uniform or field for isOnWallCoat"
            << abort(FatalError);
    }

    volScalarField& isOnWallCoat = isOnWallCoatPtr();
    
    // load the coating locations
    const dictionary& isInWallCoatDict(transportProperties.subDict("isInWallCoat"));
    word isInWallCoatType(isInWallCoatDict.lookup("type"));

    autoPtr<volScalarField> isInWallCoatPtr;

    IOobject isInWallCoatIO
    (
        "isInWallCoat",
        runTime.timeName(0),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (isInWallCoatType == "uniform")
    {
        scalar isInWallCoatValue(readScalar(isInWallCoatDict.lookup("value")));

        isInWallCoatPtr.reset
        (
            new volScalarField
            (
                isInWallCoatIO,
                mesh,
                dimensionedScalar
                (
                    "isInWallCoat",
                    dimless/dimTime,
                    isInWallCoatValue
                ),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );
    }
    else if (isInWallCoatType == "field")
    {
        isInWallCoatIO.readOpt() = IOobject::MUST_READ;

        isInWallCoatPtr.reset
        (
            new volScalarField
            (
                isInWallCoatIO,
                mesh
            )
        );
    }
    else
    {
        FatalErrorIn
        (
            "readTransportProperties.H"
        )   << "Valid type entries are uniform or field for isInWallCoat"
            << abort(FatalError);
    }

    volScalarField& isInWallCoat = isInWallCoatPtr();
    
// load the coating properties
dimensionedScalar porEpsC(transportProperties.subDict("coatingProps").lookup("porEps"));
dimensionedScalar constrC(transportProperties.subDict("coatingProps").lookup("constr"));
dimensionedScalar tortC(transportProperties.subDict("coatingProps").lookup("tort"));
dimensionedScalar kappaTC(transportProperties.subDict("coatingProps").lookup("kappaT"));
dimensionedScalar CpC(transportProperties.subDict("coatingProps").lookup("Cp"));
dimensionedScalar rhoC(transportProperties.subDict("coatingProps").lookup("rho"));
dimensionedScalar dpC(transportProperties.subDict("coatingProps").lookup("dP"));

// load the substrate properties
dimensionedScalar porEpsS(transportProperties.subDict("substrateProps").lookup("porEps"));
dimensionedScalar constrS(transportProperties.subDict("substrateProps").lookup("constr"));
dimensionedScalar tortS(transportProperties.subDict("substrateProps").lookup("tort"));
dimensionedScalar kappaTS(transportProperties.subDict("substrateProps").lookup("kappaT"));
dimensionedScalar CpS(transportProperties.subDict("substrateProps").lookup("Cp"));
dimensionedScalar rhoS(transportProperties.subDict("substrateProps").lookup("rho"));
dimensionedScalar dpS(transportProperties.subDict("substrateProps").lookup("dP"));

// load the on-wall and in-wall properties
scalar  coatingVolFracOW(readScalar(transportProperties.subDict("onWallProperties").lookup("coatingVolFrac")));
scalar  coatingVolFracIW(readScalar(transportProperties.subDict("inWallProperties").lookup("coatingVolFrac")));
    
dimensionedScalar cTemp(transportProperties.subDict("genProps").lookup("cTemp"));//current temperature
dimensionedScalar pAtm(transportProperties.subDict("genProps").lookup("pAtm"));//pressure
dimensionedScalar univR(transportProperties.subDict("genProps").lookup("univR"));//R constant
    
    
    
    
