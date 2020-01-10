/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

// declare cell definitions here 
Cell_Definition CSC; 
Cell_Definition DCC; 

Cycle_Model csc_cycle_model;

// rf. my_mutation_function in User Guide
void csc_differentiation_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// recall
	// void standard_live_phase_entry_function( Cell* pCell, Phenotype& phenotype, double dt )
        // the cell wants to double its volume
        // phenotype.volume.target_solid_nuclear *= 2.0;
        // phenotype.volume.target_solid_cytoplasmic *= 2.0;

	std::cout << "-------csc_differentiation_function ----------\n";
    phenotype.volume.target_solid_nuclear *= 2.0;
    phenotype.volume.target_solid_cytoplasmic *= 2.0;

	if (UniformRandom() < 0.5)
	{
		pCell->type = 1;
	}
	return;
}
void create_csc_cycle_model( void )
{
	std::cout << "-------create_csc_cycle_model() \n";
	csc_cycle_model.code = PhysiCell_constants::custom_cycle_model;
	csc_cycle_model.name = "CSC cycle";

	csc_cycle_model.data.time_units = "min";

	// csc_cycle_model.add_phase( PhysiCell_constants::live , "Live" );
	csc_cycle_model.add_phase( PhysiCell_constants::custom_cycle_model , "CSC cycle" );

	csc_cycle_model.phases[0].division_at_phase_exit = true;

	csc_cycle_model.add_phase_link( 0 , 0 , NULL );

	csc_cycle_model.transition_rate(0,0) = 0.0432 / 60.0; // MCF10A have ~0.04 1/hr net birth rate
	csc_cycle_model.transition_rate(0,0) *= 2.0;

	// csc_cycle_model.phases[0].entry_function = standard_live_phase_entry_function;  
	csc_cycle_model.phases[0].entry_function = csc_differentiation_function;  

	return;
}

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	// cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// set default cell cycle model 
	// cell_defaults.functions.cycle_model = live; 
	// vs.  ??
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	// set default_cell_functions; 
	
	// cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based; 
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	// _cells_cycle_model
	int live_index = live.find_phase_index( PhysiCell_constants::live );
	// cell_defaults.phenotype.cycle.data.transition_rate(live_index,live_index) *= parameters.doubles( "CSC_relative_cycle_entry_rate" );

	std::cout << "--------------- \n"; 
	std::cout << "parameters.doubles( 'CSC_relative_cycle_entry_rate' ) = " << parameters.doubles( "CSC_relative_cycle_entry_rate" ) << std::endl; 
	std::cout << "cell_defaults.phenotype.cycle.data.transition_rate(live_index,live_index) = " <<  cell_defaults.phenotype.cycle.data.transition_rate(live_index,live_index)  << std::endl;
	std::cout << "--------------- \n\n"; 

	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "CSC_apoptosis_rate" );

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38; 
	
	// add custom data here, if any 
	

	// -------- define cell types -----------
	// It's best to just copy the default and modify it. 

	live.display(std::cout);
	
	CSC = cell_defaults; 
	CSC.type = 0; 
	CSC.name = "CSC"; 

	create_csc_cycle_model();
	CSC.phenotype.cycle.sync_to_cycle_model( csc_cycle_model );
	
	CSC.functions.cycle_model = csc_cycle_model; 
	
	// csc_cycle_model = live;
	// csc_cycle_model.code = PhysiCell_constants::custom_cycle_model;
	// csc_cycle_model.name = "CSC cycle model";
	// // csc_cycle_model.phases[0].entry_function = csc_differentiation_function;
	// csc_cycle_model.phase_link(0,0).exit_function = csc_differentiation_function;
	// CSC.functions.cycle_model = csc_cycle_model; 
	// CSC.phenotype.cycle.data.transition_rate(0,0) = 0.00042; 	
	
	csc_cycle_model.display(std::cout);
	
	CSC.phenotype.sync_to_functions( CSC.functions ); 

	// rgb_model.phase_link(blue_index,red_index).exit_function = my_mutation_function;
	// live.phases[0].entry_function = standard_live_phase_entry_function;
	// CSC.phenotype.cycle.data.transition_rate(live_index,live_index) 
	// CSC.phenotype.cycle.model.ph 
	// int live_index = live.find_phase_index( PhysiCell_constants::live );
	// live.phases[0].entry_function = standard_live_phase_entry_function;    // in create_live_model()
	// CSC.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0; 	
	
	// make this cell type randomly motile, less adhesive, greater survival, 
	// and less proliferative 
	
	DCC = cell_defaults; 
	DCC.type = 1; 
	DCC.name = "DCC"; 
	
	// make sure the new cell type has its own reference phenotype
	
	// DCC.parameters.pReference_live_phenotype = &( DCC.phenotype ); 

	// DCC.phenotype.cycle.data.transition_rate(live_index,live_index) *= parameters.doubles( "DCC_relative_cycle_entry_rate" );
	
	// Set apoptosis to 
	// DCC.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "DCC_apoptosis_rate" ); // 0.0; 
	
	// Set proliferation to 10% of other cells. 
	// Alter the transition rate from G0G1 state to S state
	DCC.phenotype.cycle.data.transition_rate(live_index,live_index) *= parameters.doubles( "DCC_relative_cycle_entry_rate" );

	return; 
}

void setup_microenvironment( void )
{
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	initialize_microenvironment(); 	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pC;

	pC = create_cell(); 
	pC->assign_position( 0.0, 0.0, 0.0 );

	pC = create_cell(); 
	pC->assign_position( -300, 0, 0.0   );
	
	pC = create_cell(); 
	pC->assign_position( 0, 120, 0.0    );
	
	pC = create_cell(); 
	pC->assign_position( 0, 240, 0.0    );

	pC = create_cell(); 
	pC->assign_position( 0, 360, 0.0    );

	pC = create_cell(); 
	pC->assign_position( -120, 0, 0.0    );

	pC = create_cell(); 
	pC->assign_position( -240, 0, 0.0    );

	pC = create_cell(); 
	pC->assign_position( -300, 100, 0.0   );

	pC = create_cell(); 
	pC->assign_position( -60, 200, 0.0 );

	pC = create_cell(); 
	pC->assign_position( -400, 300, 0.0 );
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
		
	if(pCell->type == 0 )
	{
		 output[0] = "red";
		 output[2] = "red"; 
	} else if(pCell->type == 1 ){
		 output[0] = "white";
		 output[2] = "white"; 
	}
	return output; 
}
