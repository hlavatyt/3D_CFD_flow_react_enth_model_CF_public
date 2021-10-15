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
    PtrList<dimensionedScalar> DList(specieNames.size());
        
    //read the material properties for each of them
    forAll(specieNames,countI)
    {
        word currentSp = specieNames[countI];
        Info<< "Processing the specie: " << currentSp << endl;
        
        Info<< "    Reading diffusivity\n" << endl;
        // result of "write quick, write dirty"
        DList.set
        (
            countI,
            new dimensionedScalar(transportProperties.subDict(currentSp).lookup("D"))
        );
    }
    
    // load the coating locations inwall (wall of the channel)
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
    dimensionedScalar porEps(transportProperties.subDict("coatingProps").lookup("porEps"));
    dimensionedScalar constr(transportProperties.subDict("coatingProps").lookup("constr"));
    dimensionedScalar tort(transportProperties.subDict("coatingProps").lookup("tort"));
    
    dimensionedScalar cTemp(transportProperties.subDict("genProps").lookup("cTemp"));//current temperature
    cTemp.value() = -1.0;
    // Note: cTem is, at this moment, loaded just as a placeholder
    //       (I need this for ROM)
    
    
    
    
    
