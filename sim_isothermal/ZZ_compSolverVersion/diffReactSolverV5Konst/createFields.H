    Info<< "Reading the concentrations \n" << endl;

    // prepare list to store the concentrations in
    PtrList<volScalarField> MFrList(specieNames.size());

    forAll(specieNames,countI)
    {
        word currentSp = specieNames[countI];
        Info<< "Processing the specie: " << currentSp << endl;

        MFrList.set
        (
            countI,
            new volScalarField
            (
                IOobject
                (
                    currentSp,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }

    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"
