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

double prob_p = 0.5;
double one_minus_prob_p = 0.5;

void check_cell_cycle_and_differentiate( Cell* pCell, Phenotype& phenotype, double dt )
{
	double elapsed_time = phenotype.cycle.data.elapsed_time_in_phase;
	// if ((elapsed_time < 0.6) && (pCell->type == 0))
	// phenotype.execute_cycle_phase_entry_function = true; 
	// std::cout << __FUNCTION__ << ": pCell->ID = "<< pCell->ID << ", elapsed_time_in_phase = " << elapsed_time << ", current_time = " << PhysiCell_globals.current_time<< std::endl;
	// if (phenotype.execute_cycle_phase_entry_function && elapsed_time < 0.6) 
	if (elapsed_time < 0.1) 
	{
		// std::cout << __FUNCTION__ << "pCell->ID = "<< pCell->ID << ", elapsed_time_in_phase = " << elapsed_time << ", current_time = " << PhysiCell_globals.current_time<< std::endl;
		std::cout << "pCell->ID = "<< pCell->ID << ", elapsed_time_in_phase = " << elapsed_time;
		std::cout << ", current_time = " << PhysiCell_globals.current_time<< std::endl;
		std::cout << ", flagged_for_division = " << phenotype.flagged_for_division << std::endl;

		if (UniformRandom() <= one_minus_prob_p)
		{
			// pCell->type = 1;
			// if (pCell->type == 0)
			std::cout << "*******  CSC --> DCC\n";
			pCell->convert_to_cell_definition(DCC);
		}

		phenotype.execute_cycle_phase_entry_function = false;
	}
}

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 

	prob_p = parameters.doubles("prob_p"); 
	one_minus_prob_p = 1.0 - prob_p; // rwh: probability that CSC becomes a DCC
	std::cout << "--------- prob_p = " << prob_p << std::endl;
	std::cout << "--------- 1 - prob_p = " << one_minus_prob_p << std::endl;

	initialize_default_cell_definition();
	
	// set default cell cycle model 
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live );   // rwh: can I override this later??
	// set default_cell_functions; 
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 

	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// turn off death
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0;


	int live_index = live.find_phase_index( PhysiCell_constants::live );

	std::cout << "--------------- \n"; 
	std::cout << "parameters.doubles( 'CSC_relative_cycle_entry_rate' ) = " << parameters.doubles( "CSC_relative_cycle_entry_rate" ) << std::endl; 
	std::cout << "cell_defaults.phenotype.cycle.data.transition_rate(live_index,live_index) = " <<  cell_defaults.phenotype.cycle.data.transition_rate(live_index,live_index)  << std::endl;
	std::cout << "--------------- \n\n"; 


	// set oxygen uptake / secretion parameters for the default cell type 
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
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
	CSC.functions.update_phenotype = check_cell_cycle_and_differentiate; 
	

	DCC = cell_defaults; 
	DCC.type = 1; 
	DCC.name = "DCC"; 
	
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

	pC = create_cell(CSC); 
	pC->assign_position( 0, -200, 0.0    );
	pC = create_cell(CSC); 
	pC->assign_position( -200, 0, 0.0    );
	pC = create_cell(CSC); 
	pC->assign_position( 0, 200, 0.0    );
	pC = create_cell(CSC); 
	pC->assign_position( 200, 0, 0.0    );
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	// /modules/*_pathology.cpp
	// std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
	std::vector< std::string > output( 4 , "rgb(0,0,0)" );
		
	if(pCell->type == 0 )
	{
		 output[0] = "red";
		 output[2] = "red"; 
	} else if(pCell->type == 1 ){
		 output[0] = "cyan";
		 output[2] = "cyan"; 
	}
	return output; 
}
