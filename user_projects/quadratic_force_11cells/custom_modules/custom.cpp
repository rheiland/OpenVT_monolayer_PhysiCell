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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

double time_90pct;  // time at which the outer cells reach 90% max width
double xpos_90pct;

void create_cell_types( void )
{
	// set the random seed 
	// SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	// cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
    cell_defaults.functions.cell_division_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	// cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_cell_rule; 
    cell_defaults.functions.update_velocity = custom_update_cell_velocity;

    xpos_90pct = (0.9 * cell_defaults.phenotype.geometry.radius * 2 * 10.0) / 2.0;
    std::cout << "-------- xpos_90pct= " << xpos_90pct << std::endl;
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

// do every mechanics dt
void custom_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{
    pCell->custom_data["cell_ID"] = pCell->ID;

    // Chaste RepulsionForce / GeneralisedLinearSpringForce cell-cell repulsion.
    // Parameters matching Relax11.cpp:
    //   mu    = 5.0  (SetMeinekeSpringStiffness)
    //   alpha = 5.0  (hard-coded in GeneralisedLinearSpringForce.cpp, unused in pure repulsion)
    // Note: Chaste uses damping = 1.0 by default, so velocity = force / 1.0 = force.
    //       Set PhysiCell's damping coefficient to 1.0 in the XML config to match.
    // static double mu    = 5.0;
    // static double alpha = 5.0;  // exponential decay constant (attraction branch only)

	static double m_stiffness = parameters.doubles("m_stiffness");   // 5, 10, ?

    double r_a = phenotype.geometry.radius;

    // Reset velocity; motility not included here (matches Chaste's RelaxationForce-only setup)
    pCell->velocity = {0.0, 0.0, 0.0};

    for( auto pNeighbor : *PhysiCell::all_cells )   // rwh: optimize
    {
        if( pNeighbor == pCell ) continue;

        double r_b = pNeighbor->phenotype.geometry.radius;
        double rest_length = r_a + r_b;

        if(std::abs(pCell->position[0]-pNeighbor->position[0]) > rest_length) continue;
        else if(std::abs(pCell->position[1]-pNeighbor->position[1]) > rest_length) continue;

        // if (pCell->ID <= 1)
        // {
        //     std::cout << "---   nbr ID= " << pNeighbor->ID << ", rest_length= " << rest_length << std::endl;
        // }

        // Vector from A (pCell) toward B (pNeighbor)
        double dx = pNeighbor->position[0] - pCell->position[0];
        double dy = pNeighbor->position[1] - pCell->position[1];
        double dz = 0.; // pNeighbor->position[2] - pCell->position[2];
        //jdouble distance = std::sqrt( dx*dx + dy*dy + dz*dz );
        double distance = std::sqrt( dx*dx + dy*dy);

        if( distance < 1e-12 ) continue;  // avoid division by zero

        double overlap = distance - rest_length;  // negative when cells overlap

        if( overlap >= 0.0 )
        {
            continue;

            // Attraction branch (linear * exponential decay); only reached if
            // distance >= rest_length -- RepulsionForce skips this, but kept
            // for completeness to match GeneralisedLinearSpringForce exactly.

            // double magnitude = mu * overlap * std::exp( -alpha * overlap / rest_length );
            // pCell->velocity[0] += magnitude * dx / distance;
            // pCell->velocity[1] += magnitude * dy / distance;
            // pCell->velocity[2] += magnitude * dz / distance;
        }
        else
        {
            // Repulsion branch: overlap/rest_length in (-1, 0)
            // log argument in (0, 1) -> magnitude < 0 -> force points away from B

            // logarithmic (matches Chaste's "stable" form)
            // double magnitude = mu * std::log( 1.0 + overlap / rest_length );

            // quadratic magnitude, preserving sign for direction
            double magnitude = m_stiffness * overlap * std::abs(overlap) / rest_length;

            pCell->velocity[0] += magnitude * dx / distance;
            pCell->velocity[1] += magnitude * dy / distance;
            // pCell->velocity[2] += magnitude * dz / distance;
            // pCell->velocity[2] = 0.0;
        }

    }
}

// called every dt_mech; pC->functions.custom_cell_rule( pC,pC->phenotype,time_since_last_mechanics );
void custom_cell_rule( Cell* pCell, Phenotype& phenotype , double dt )
{ 
    static bool reached_90 = false;
    std::stringstream ss;

    // if (pCell->ID == 10 && pCell->position[0] >= 90.0)  // 90%, 90pct
    if (pCell->ID == 10)
    {
        // std::cout << "------- " << __FUNCTION__ << ":  ID= " << pCell->ID <<":  x= " << pCell->position[0] << std::endl;
        // full width: -5*d to 5*d, e.g. d=10 -> 50-(-50)=100; 90% is 90; divide by 2 = 45 (x-position)
        // full width: -5*d to 5*d, e.g. d=10 -> 50-(-50)=100; 90% is 90; divide by 2 = 45 (x-position)
        // if (pCell->position[0] >= 45.0)  // for symmetric test (hard-coded for diam=10)
        if (pCell->position[0] >= xpos_90pct)  // for symmetric test (hard-coded for diam=10)
        {
            if (!reached_90)
            {
                std::cout <<"---- "<< __FUNCTION__ << ": cell radius = " << phenotype.geometry.radius << std::endl;
                time_90pct = PhysiCell_globals.current_time;
                std::cout <<"---- "<< __FUNCTION__ << ": Width reached 90% , t= " << time_90pct << std::endl;

                std::ofstream outFile;
                // sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );
                ss << PhysiCell_settings.folder.c_str()  << "/time_90pct.txt" ;
                std::string out_file = ss.str();
                std::cout <<"---- "<< __FUNCTION__ << ": time_pct out_file=" << out_file << std::endl;
                outFile.open(out_file);   // read by  ../analysis/plot_11cells_sweep.py  for plotting results
                outFile << time_90pct;
                outFile.close();

                reached_90 = true;
                // std::exit();
            }
            else
            {
                if (PhysiCell_globals.current_time / time_90pct > 10.0)
                {
                    std::cout <<"---- "<< __FUNCTION__ << "  calibrated T > 10.  Exit simulation!" << std::endl;
                    std::exit(-1);
                }

            }
        }
    }
} 