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
#include "fvOptions.H"
#include "Switch.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readTransportProperties.H"
    
    simpleControl simple(mesh);
    
    #include "createFields.H"
    #include "createFvOptions.H"                                        
    
    #include "CourantNo.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    
    // while (simple.loop())
    while  (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl; 
        
        Info<< "\nCalculating the concentrations\n" <<endl;
        
        while (simple.correctNonOrthogonal())
        {
            // solve the equations for the concentrations
            #include "concEq.H"
        }
        
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
        
        
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
