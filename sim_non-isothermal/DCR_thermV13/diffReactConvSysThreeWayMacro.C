/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    diffReactConvSysThreeWayMacro

Description
    Diffusion-reaction-convection solver for macro-scale model of a 
    monolithic catalytic filter. Three way catalytic coating is assumed.
    
    Solver can process in-wall and on-wall coating separately with
    different reaction constants (based on the volumetric fraction of 
    catalytic coating in on-wall and in-wall layer)
    
    In-wall/on-wall regions are stored in different cellZones

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "simpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "IOporosityModelList.H"
#include "Switch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //~ #define CREATE_FIELDS_2 createPorousZones.H

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "readTransportPropertiesMy.H"
    #include "createFields.H"
    //~ #include "createPorousZones.H"
    #include "createZones.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "postProcess.H"
    
    
    turbulence->validate();
    //~ simpleControl simple(mesh);
                                       
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    
    while (simple.loop(runTime))
    // while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl; 
        
        
    //    while (simple.correctNonOrthogonal())
        {
            Info << "\nSolving NS equations"<<endl;
            
            #include "UEqn.H"

            
            // solve the equations for the concentrations
            Info << "\nSolving concentration equation "<<endl;
            
            #include "concEq.H"
            
            for (label auxI=0;auxI < 3;auxI++)
            {   
                Info << "\nCorrection of temperature"<<endl;
                #include "EEqn.H"
            }
            
            Info << "\nCorrection of pressure"<<endl;
            #include "pEqn.H"   
            //~ for (label auxI=0;auxI < 3;auxI++)
            //~ {   
                
        }              
                
            //~ }
            turbulence->correct();

        
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
        
        
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
