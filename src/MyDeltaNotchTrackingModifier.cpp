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

#include "MyDeltaNotchTrackingModifier.hpp"
#include "MyDeltaNotchSrnModel.hpp"
#include "Debug.hpp"

template<unsigned DIM>
MyDeltaNotchTrackingModifier<DIM>::MyDeltaNotchTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
MyDeltaNotchTrackingModifier<DIM>::~MyDeltaNotchTrackingModifier()
{
}

template<unsigned DIM>
void MyDeltaNotchTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void MyDeltaNotchTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void MyDeltaNotchTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();
    c_vector<double,2> population_centroid = rCellPopulation.GetCentroidOfCellPopulation();
    // First recover each cell's Notch and Delta concentrations from the ODEs and store in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        MyDeltaNotchSrnModel* p_model = static_cast<MyDeltaNotchSrnModel*>(cell_iter->GetSrnModel());
        double this_delta                             = p_model->GetDelta();
        double this_cell_surface_notch                = p_model->GetCellSurfaceNotch();
        double this_sudx_dependent_notch              = p_model->GetSudxDependentNotch();
        double this_dx_dependent_early_endosome_notch = p_model->GetDxDependentEarlyEndosomeNotch();
        double this_dx_dependent_late_endosome_notch  = p_model->GetDxDependentLateEndosomeNotch();
        double this_notch_intracellular_domain        = p_model->GetNotchIntracellularDomain();

        c_vector<double,2> this_centroid = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        c_vector<double,2> this_distance_to_tissue_centre;
        this_distance_to_tissue_centre(0) = 0.0;
        this_distance_to_tissue_centre(1) = 0.0;
        this_distance_to_tissue_centre = this_centroid - population_centroid;
        double this_x_distance = fabs(this_distance_to_tissue_centre[0]);

        double total_notch = this_cell_surface_notch + this_sudx_dependent_notch +
                             this_dx_dependent_early_endosome_notch + this_dx_dependent_late_endosome_notch +
                             this_notch_intracellular_domain;

        // Note that the state variables must be in the same order as listed in DeltaNotchOdeSystem
        cell_iter->GetCellData()->SetItem("cell surface notch", this_cell_surface_notch);
        cell_iter->GetCellData()->SetItem("sudx dependent notch", this_sudx_dependent_notch);
        cell_iter->GetCellData()->SetItem("dx dependent early endosome notch", this_dx_dependent_early_endosome_notch);
        cell_iter->GetCellData()->SetItem("dx dependent late endosome notch", this_dx_dependent_late_endosome_notch);
        cell_iter->GetCellData()->SetItem("notch intracellular domain", this_notch_intracellular_domain);
        cell_iter->GetCellData()->SetItem("total notch", total_notch);
        cell_iter->GetCellData()->SetItem("delta", this_delta);
        cell_iter->GetCellData()->SetItem("x distance", this_x_distance);
    }

    // Next iterate over the population to compute and store each cell's neighbouring Delta concentration in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the set of neighbouring location indices
        std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);

        // Compute this cell's average neighbouring Delta concentration and store in CellData
        if (!neighbour_indices.empty())
        {
            double mean_delta = 0.0;
            for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                 iter != neighbour_indices.end();
                 ++iter)
            {
                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                double this_delta = p_cell->GetCellData()->GetItem("delta");
                mean_delta += this_delta/neighbour_indices.size();
            }
            cell_iter->GetCellData()->SetItem("mean delta", mean_delta);
        }
        else
        {
            // If this cell has no neighbours, such as an isolated cell in a CaBasedCellPopulation, store 0.0 for the cell data
            cell_iter->GetCellData()->SetItem("mean delta", 0.0);
        }
    }
}

template<unsigned DIM>
void MyDeltaNotchTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class MyDeltaNotchTrackingModifier<1>;
template class MyDeltaNotchTrackingModifier<2>;
template class MyDeltaNotchTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MyDeltaNotchTrackingModifier)
