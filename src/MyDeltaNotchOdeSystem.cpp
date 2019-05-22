/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "MyDeltaNotchOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"
#include "Debug.hpp"

MyDeltaNotchOdeSystem::MyDeltaNotchOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(6)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<MyDeltaNotchOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Notch concentration for this cell
     * 1 - Delta concentration for this cell
     *
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    SetDefaultInitialCondition(0, 1.0); // soon overwritten
    SetDefaultInitialCondition(1, 1.0); // soon overwritten
    SetDefaultInitialCondition(2, 1.0); // soon overwritten
    SetDefaultInitialCondition(4, 1.0); // soon overwritten
    SetDefaultInitialCondition(3, 1.0); // soon overwritten
    SetDefaultInitialCondition(5, 1.0); // soon
    this->mParameters.push_back(0.5);

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

MyDeltaNotchOdeSystem::~MyDeltaNotchOdeSystem()
{
}

void MyDeltaNotchOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    // first define each dynamic component of the system
    double cell_surface_notch = rY[0];
    double sudx_dependent_notch = rY[1];
    double dx_dependent_early_endosome_notch = rY[2];
    double dx_dependent_late_endosome_notch = rY[3];
    double notch_intracellular_domain = rY[4];
    double delta = rY[5];
    double mean_delta = this->mParameters[0]; // Shorthand for "this->mParameter("mean delta");"
    double x_distance = this->mParameters[1];

    // define the components of the fluxes
    double k_1 = 14.0;
    double k_2 = 10.0;
    double k_3 = 240.0;
    double k_4 = 420.0;
    double k_5 = 100.0;
    double k_6 = 500.0;
    double k_7 = 15.0;
    double k_8 = 1.2;
    double k_9 = 108.0;
    double k_10 = 250.0;
    double k_11 = 1.0;
    double k_12 = 70.0;
    double k_13 = 0.06;
    double c_3 = 14.0;
    double c_4 = 14.0;
    double c_8a = 14.0;
    double c_8b = 14.0;
    double c_9 = 14.0;
    double c_10 = 14.0;
    double beta_N = 10.0;
    double f = 5.0;
    double k_c = 0.001;
    double fb_D = 10.0;
    double fb_N = 10.0;
    double fb_5 = 10.0;
    double fb_10 = 10.0;
    // double f_bs = 10.0;
    double gamma = 0.25;
    double dx = 10.0;
    double sudx = 10.0;

    double beta_D = beta_N * x_distance * (1 - f/12) * (fb_D / (fb_D + notch_intracellular_domain));

    // define the fluxes using the above components
    double r_1 = k_1 * (2 - fb_N/(fb_N + notch_intracellular_domain));
    double r_2 = k_2 * cell_surface_notch;
    double r_3 = (k_3 * sudx + c_3) * cell_surface_notch;
    double r_4 = (k_4 * dx + c_4) * cell_surface_notch;
    double r_5 = k_5 * sudx * (1 - fb_5/(fb_5 + delta)) * dx_dependent_early_endosome_notch;
    double r_6 = k_6 * mean_delta * cell_surface_notch; // need to incorporate neighbour delta
    double r_7 = k_7 * sudx_dependent_notch;
    double r_8 = k_8 * dx_dependent_early_endosome_notch + (c_8a * dx_dependent_early_endosome_notch) / (c_8b + dx_dependent_early_endosome_notch);
    double r_9 = (k_9 * sudx + c_9) * dx_dependent_late_endosome_notch;
    double r_10 = (k_10 * sudx + c_10) * (1 - fb_10/(fb_10 + delta)) * sudx_dependent_notch;
    double r_11 = k_11 * dx_dependent_early_endosome_notch;
    double r_12 = k_12 * dx_dependent_late_endosome_notch;
    double r_13 = k_13 * notch_intracellular_domain;
    double r_c = cell_surface_notch * delta/k_c;

    // The next 6 lines define the ODE system by Shimizu et al. (2014)
    rDY[0] = r_1 - r_2 - r_3 - r_4 - r_6;  // d[Notch_1]/dt
    rDY[1] = r_3 + r_5 + r_7 - r_10;  // d[Notch_2]/dt
    rDY[2] = r_4 - r_5 - r_8 - r_11;  // d[Notch_3]/dt
    rDY[3] = r_8 - r_9 - r_12;  // d[Notch_4]/dt
    rDY[4] = r_6 + r_7 + r_9 - r_13;  // d[NICD]/dt
    rDY[5] = beta_D - gamma*delta - r_6 - r_c;  // d[Delta]/dt
}

template<>
void CellwiseOdeSystemInformation<MyDeltaNotchOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("cell surface notch");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("sudx dependent notch");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("dx dependent early endosome notch");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("dx dependent late endosome notch");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("notch intracellular domain");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("delta");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    // If this is ever not the first parameter change the line
    // double mean_delta = this->mParameters[0]; in EvaluateYDerivatives().
    this->mParameterNames.push_back("mean delta");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("x distance");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MyDeltaNotchOdeSystem)
