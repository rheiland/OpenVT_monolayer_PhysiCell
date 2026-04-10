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

#include <algorithm>  // to "find" element in vector
#include "./custom.h"

// int param_matrix_elm; // 0,1,2 (top row); 3,4,5 (middle); 6,7,8 (bottom)
double beta_threshold, gamma_threshold;

double target_volume;
double target_area;

void create_cell_types( void )
{
	// set the random seed 
	// SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	// param_matrix_elm = parameters.ints("param_matrix_elm");  
	beta_threshold = parameters.doubles("beta_threshold");  
	gamma_threshold = parameters.doubles("gamma_threshold");  
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// cell_defaults.functions.volume_update_function = standard_volume_update_function;
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

	// cell_defaults.functions.custom_cell_rule = custom_function;   // do every mechanics dt
	cell_defaults.functions.custom_cell_rule = custom_cell_rule;   // do every mechanics dt

	// cell_defaults.functions.contact_function = contact_function; 

	cell_defaults.functions.update_velocity = custom_update_cell_velocity;

    // rwh: simple test - hard-coded division direction
    // cell_defaults.functions.cell_division_direction_function = custom_division_dir_function; 

    cell_defaults.functions.cell_division_function = custom_division_function; 
    cell_defaults.functions.volume_update_function = custom_volume_function; 

    // target_volume = cell_defaults.phenotype.volume.total;  // 523.6
    target_area = 3.1415927 * cell_defaults.phenotype.geometry.radius * cell_defaults.phenotype.geometry.radius;
    // pCell->custom_data["cell_area"] =  3.1415927 * pCell->custom_data["cell_radius"] * pCell->custom_data["cell_radius"];
    // std::cout << "target_volume  = " << target_volume <<std::endl;
    std::cout << "target_area  = " << target_area <<std::endl;   // 78.5399
    // rwh: uh, but this will change, based on the norm_randd_norm

	
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


// cell_defaults.functions.update_velocity = standard_update_cell_velocity;

// in container:  (kinda strange, but don't even need this for the monolayer? Could update pos here!)
		// // update velocities 
		// #pragma omp parallel for 
		// for( int i=0; i < (*all_cells).size(); i++ )
		// {
		// 	Cell* pC = (*all_cells)[i]; 
		// 	if( pC->functions.update_velocity && pC->is_out_of_domain == false && pC->is_movable )
		// 	{ pC->functions.update_velocity( pC,pC->phenotype,time_since_last_mechanics ); }
		// }

		// // new March 2023: 
		// // dynamic spring attachments, followed by built-in springs
		// if( PhysiCell_settings.disable_automated_spring_adhesions == false )
		// {
		// 	#pragma omp parallel for 
		// 	for( int i=0; i < (*all_cells).size(); i++ )
		// 	{
		// 		Cell* pC = (*all_cells)[i]; 
		// 		dynamic_spring_attachments(pC,pC->phenotype,time_since_last_mechanics); 
		// 	}		
		// 	#pragma omp parallel for 
		// 	for( int i=0; i < (*all_cells).size(); i++ )
		// 	{
		// 		Cell* pC = (*all_cells)[i]; 
		// 		if( pC->is_movable )
		// 		{
		// 			for( int j=0; j < pC->state.spring_attachments.size(); j++ )
		// 			{
		// 				Cell* pC1 = pC->state.spring_attachments[j]; 
		// 				standard_elastic_contact_function(pC,pC->phenotype,pC1,pC1->phenotype,time_since_last_mechanics);  
		// 			}
		// 		}
		// 	}	
		// }

		// // new March 2022: 
		// // run standard interactions (phagocytosis, attack, fusion) here 
		// #pragma omp parallel for 
		// for( int i=0; i < (*all_cells).size(); i++ )
		// {
		// 	Cell* pC = (*all_cells)[i]; 
		// 	standard_cell_cell_interactions(pC,pC->phenotype,time_since_last_mechanics); 
		// }
		// // super-critical to performance! clear the "dummy" cells from phagocytosis / fusion
		// // otherwise, comptuational cost increases at polynomial rate VERY fast, as O(10,000) 
		// // dummy cells of size zero are left ot interact mechanically, etc. 
		// if( cells_ready_to_die.size() > 0 )
		// {
		// 	for( int i=0; i < cells_ready_to_die.size(); i++ )
		// 	{ cells_ready_to_die[i]->die(); }
		// 	cells_ready_to_die.clear();
		// }

		// // update positions 

//----------------------------------------
void generate_neighbor_list( Cell* pCell)
{
    // static double tmin= 1330;
    // static double tmin= 1.0e99;
    static double tmin= 4.0;
    static double tmax= 1331;
    // if (PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
    //     std::cout << "-------- " <<__FUNCTION__ << ":  t= "<<PhysiCell_globals.current_time << ": pCell->ID= " << pCell->ID<< std::endl;

	pCell->state.neighbors.clear();
    std::vector<PhysiCell::Cell *> true_nbrs;

    // pCell->custom_data["num_nbrs"] = 0;
	
	// ---------- First check the neighbors in my current voxel --------------
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		// add_potentials_monolayer(pCell, *neighbor);

        Cell* pN = *neighbor;
        if (pCell->ID == pN->ID) continue;

        bool nbr_exists = false;
        for (const auto& nbr : pCell->state.neighbors)
        {
            if (nbr->ID == pN->ID)
            {
                nbr_exists = true;
                break;
            }
        }
        if (!nbr_exists)
            pCell->state.neighbors.push_back(*neighbor);

        // if (PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
        //     std::cout << "   (current voxel) ID= " << pCell->ID<< ",  neighbor ID=" << pN->ID << std::endl;
        // std::cout << __FUNCTION__ << "  t= "<<PhysiCell_globals.current_time << ": (current voxel) ID= " << pCell->ID<< ",  neighbor ID=" << pN->ID << ", pCell num_nbrs= "<< pCell->custom_data["num_nbrs"] << std::endl;
                // if (fabs( - 120.0) < 1.e-4)
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
            Cell* pN = *neighbor;
            if (pCell->ID == pN->ID) continue;

            bool nbr_exists = false;
            for (const auto& nbr : pCell->state.neighbors)
            {
                if (nbr->ID == pN->ID)
                {
                    nbr_exists = true;
                    // if (pCell->ID == 1 && PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
                    // {
                    //     std::cout << "       -- for ID=1, nbr_exists ID=" << nbr->ID << std::endl;
                    // }
                    break;
                }
            }
            if (!nbr_exists)
            {
                pCell->state.neighbors.push_back(*neighbor);
                // if (pCell->ID == 1 && PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
                // {
                //     std::cout << "       -- for ID=1, push_back ID=" << (*neighbor)->ID << std::endl;
                // }
            }
            // else
            // {
            //     if (pCell->ID == 1 && PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
            //     {
            //         std::cout << "   ------ for ID=1, nbr_exists ID=" << (*neighbor)->ID << ", so not adding " <<std::endl;
            //     }
            // }

            // rwh_standard_update_cell_velocity  t= 0: (current voxel) ID= 0,  neighbor ID=0, pCell num_nbrs= 1
            // if (PhysiCell_globals.current_time < 0.03)
            // if (PhysiCell_globals.current_time > 1030.)
            // if (PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
            //     std::cout << "   (nbr voxel) ID= " << pCell->ID<< ",  neighbor ID=" << pN->ID << std::endl;
            //     true_nbrs.push_back(*neighbor);
                // std::cout << __FUNCTION__ << "  t= "<<PhysiCell_globals.current_time << ": (outer voxel) ID= " << pCell->ID<< ",  neighbor ID=" << pN->ID << std::endl;

			// add_potentials_monolayer(pCell, *neighbor);
            // if (PhysiCell_globals.current_time > 1030.)
            //     std::cout << "      -- after add_potentials_monolayer(): (outer voxel) ID= " << pCell->ID<< ",  neighbor ID=" << pN->ID << ", pCell num_nbrs= "<< pCell->custom_data["num_nbrs"] << std::endl;
		}
	}

    // if (PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
    // {
    //     std::cout << " -- nbrs, minus duplicates:" << std::endl;
    //     for (const auto& nbr : pCell->state.neighbors)
    //         std::cout <<  nbr->ID << std::endl;
    // }

    // if (pCell->ID == 1 && PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
    // {
    //     std::cout << " -------- ID 1 has prelim nbrs: " << std::endl;
    //     for (const auto& nbr : pCell->state.neighbors)
    //         std::cout <<  nbr->ID << std::endl;
    // }


    // -----------  OK, now let's cull the nbrs to match our distance constraint  --------------
    for (const auto& nbr : pCell->state.neighbors)
    {
        double touching_dist = pCell->phenotype.geometry.radius + nbr->phenotype.geometry.radius;

        double dx = nbr->position[0] - pCell->position[0];
        double dy = nbr->position[1] - pCell->position[1];
        // double dz = 0.; // pNeighbor->position[2] - pCell->position[2];
        //jdouble distance = std::sqrt( dx*dx + dy*dy + dz*dz );
        double distance = std::sqrt( dx*dx + dy*dy);
        if (distance > touching_dist)
        {
            // if (pCell->ID == 1 && PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
            // {
            //     std::cout << " -- for ID=1, removing ID " << nbr->ID << std::endl;
            // }

            auto new_end = std::remove(pCell->state.neighbors.begin(), pCell->state.neighbors.end(), nbr);
            // Use vector::erase to remove the elements from the new logical end to the actual end.
            pCell->state.neighbors.erase(new_end, pCell->state.neighbors.end());
        }
    }

    // if (pCell->ID == 1 && PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
    // {
    //     std::cout << " ----- ID 1 has true nbrs, overlapping :" << std::endl;
    //     for (const auto& nbr : pCell->state.neighbors)
    //         std::cout <<  nbr->ID << std::endl;
    // }

}

// do every mechanics dt
void custom_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{
    static double effective_repulsion = phenotype.mechanics.cell_cell_repulsion_strength;  // e.g., 10
    static double rest_length_factor = parameters.doubles("rest_length_factor");

    pCell->custom_data["cell_ID"] = pCell->ID;

    generate_neighbor_list(pCell);   // <---------------- NEW (rwh)

    // Chaste RepulsionForce / GeneralisedLinearSpringForce cell-cell repulsion.
    // Parameters matching Relax11.cpp:
    //   mu    = 5.0  (SetMeinekeSpringStiffness)
    //   alpha = 5.0  (hard-coded in GeneralisedLinearSpringForce.cpp, unused in pure repulsion)
    // Note: Chaste uses damping = 1.0 by default, so velocity = force / 1.0 = force.
    //       Set PhysiCell's damping coefficient to 1.0 in the XML config to match.
    // static double mu    = 5.0;
    // static double alpha = 5.0;  // exponential decay constant (attraction branch only)

    // static double mStiffness = 10;  // 5.0;
    // static double m_stiffness = parameters.doubles("m_stiffness");   // 5, 10, ?

    double r_a = phenotype.geometry.radius;

    // Reset velocity; motility not included here (matches Chaste's RelaxationForce-only setup)
    pCell->velocity = {0.0, 0.0, 0.0};

    pCell->custom_data["num_nbrs"] = 0;
    // int nbr_ID = -1;
    // double t_saved;
    // for( auto pNeighbor : *PhysiCell::all_cells )   // rwh: optimize
    for (const auto& pNeighbor : pCell->state.neighbors)
    {
        // if( pNeighbor == pCell ) continue;

        double r_b = pNeighbor->phenotype.geometry.radius;
        // double rest_length = r_a + r_b;
        double rest_length = (r_a + r_b) * rest_length_factor;

        // if(std::abs(pCell->position[0]-pNeighbor->position[0]) > rest_length) continue;
        // else if(std::abs(pCell->position[1]-pNeighbor->position[1]) > rest_length) continue;

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
            pCell->custom_data["num_nbrs"] += 1;   // now performed in custom_cell_rule[_slow], after update_pos

            // if ((pCell->custom_data["num_nbrs"] == 1) && pCell->ID == 1)
            // {
            //     // std::cout << "---   cell ID= " << pCell->ID << " has just 1 nbr: ID = " << pNeighbor->ID << ", t="<<PhysiCell_globals.current_time << std::endl;
            //     nbr_ID = pNeighbor->ID;
            //     t_saved = PhysiCell_globals.current_time;
            //     // true_nbrs.push_back();

            //     // std::cout << "---   cell ID= " << pCell->ID << " had nbrs = " << pCell->state.neighbors <<", latest has ID = " << pNeighbor->ID << std::endl;
            //     // if (pNeighbor->ID != 0)
            //         // std::cout << "---   cell ID= " << pCell->ID << " has nbr ID = " << pNeighbor->ID << std::endl;
            // }
            // Repulsion branch: overlap/rest_length in (-1, 0)
            // log argument in (0, 1) -> magnitude < 0 -> force points away from B

            // logarithmic (matches Chaste's "stable" form)
            // double magnitude = mu * std::log( 1.0 + overlap / rest_length );

            // quadratic magnitude, preserving sign for direction
            // double magnitude = m_stiffness * overlap * std::abs(overlap) / rest_length;

            double temp = 1.0 - (distance / rest_length);   // normalized overlap fraction
            // double magnitude = stiffness * temp * temp;
            double magnitude = effective_repulsion * temp * temp;

            // --- cf. PhysiCell /core: void Cell::add_potentials(Cell* other_agent)
            //	velocity[i] += displacement[i] * temp_r; 
    	    // axpy( &velocity , temp_r , displacement );    # core/PhysiCell_cell.cpp

            // pCell->velocity[0] += magnitude * dx / distance;
            // pCell->velocity[1] += magnitude * dy / distance;

            pCell->velocity[0] -= magnitude * dx / distance;  // yipeee: negate for repulsion
            pCell->velocity[1] -= magnitude * dy / distance;  // yipeee: negate for repulsion

            // pCell->velocity[2] += magnitude * dz / distance;
            // pCell->velocity[2] = 0.0;
        }
    }
    // if ((pCell->custom_data["num_nbrs"] == 1) && pCell->ID == 1)
    // {
    //     std::cout << "--- FINAL:   cell ID= " << pCell->ID << " has 1 nbr: ID = " << nbr_ID << ", t_saved="<< t_saved << std::endl;
    //     // std::cout << "---   cell ID= " << pCell->ID << " had nbrs = " << pCell->state.neighbors <<", latest has ID = " << pNeighbor->ID << std::endl;
    //     // if (pNeighbor->ID != 0)
    //         // std::cout << "---   cell ID= " << pCell->ID << " has nbr ID = " << pNeighbor->ID << std::endl;
    // }
}

// do every mechanics dt
// void custom_update_cell_velocity_old( Cell* pCell, Phenotype& phenotype, double dt )
// {
//     pCell->custom_data["cell_ID"] = pCell->ID;

//     // Chaste RepulsionForce / GeneralisedLinearSpringForce cell-cell repulsion.
//     // Parameters matching Relax11.cpp:
//     //   mu    = 5.0  (SetMeinekeSpringStiffness)
//     //   alpha = 5.0  (hard-coded in GeneralisedLinearSpringForce.cpp, unused in pure repulsion)
//     // Note: Chaste uses damping = 1.0 by default, so velocity = force / 1.0 = force.
//     //       Set PhysiCell's damping coefficient to 1.0 in the XML config to match.
//     // static double mu    = 5.0;
//     // static double alpha = 5.0;  // exponential decay constant (attraction branch only)

//     // static double mStiffness = 10;  // 5.0;
//     static double m_stiffness = parameters.doubles("m_stiffness");   // 5, 10, ?

//     double r_a = phenotype.geometry.radius;

//     // Reset velocity; motility not included here (matches Chaste's RelaxationForce-only setup)
//     pCell->velocity = {0.0, 0.0, 0.0};

//     // Iterate over neighboring cells via the mechanics grid (Moore neighborhood)
//     // Cell_Container* container = (Cell_Container*) pCell->get_container();
//     // int mech_voxel = pCell->current_mechanics_voxel_index;
//     // int mech_voxel = pCell->current_mechanics_voxel_index;
//     // int mech_voxel = get_container()->underlying_mesh.nearest_voxel_index( pCell->position );
//     // int mech_voxel = container->underlying_mesh.nearest_voxel_index( pCell->position );

//     // if (pCell->ID <= 1)
//     // {
//     //     std::cout << __FUNCTION__ << "-------   ID= " << pCell->ID << std::endl;
//     // }

//     // for( int n_vox : container->underlying_mesh.moore_connected_voxel_indices[mech_voxel] )
//     // {
//     //     for( Cell* pNeighbor : container->agent_grid[n_vox] )

//     // for( auto cell : *PhysiCell::all_cells )
//         double gamma = 0.0;
//         double beta = 0.0;
//     {
//         static double pi = 3.1415926535897932384626433832795; 
//         static double one_div_pi = 0.31830989; 

//         for( auto pNeighbor : *PhysiCell::all_cells )   // rwh: optimize
//         {
//             if( pNeighbor == pCell ) continue;


//             double r_b = pNeighbor->phenotype.geometry.radius;
//             double rest_length = r_a + r_b;

//             if(std::abs(pCell->position[0]-pNeighbor->position[0]) > rest_length) continue;
//             else if(std::abs(pCell->position[1]-pNeighbor->position[1]) > rest_length) continue;

//             // if (pCell->ID <= 1)
//             // {
//             //     std::cout << "---   nbr ID= " << pNeighbor->ID << ", rest_length= " << rest_length << std::endl;
//             // }

//             // Vector from A (pCell) toward B (pNeighbor)
//             double dx = pNeighbor->position[0] - pCell->position[0];
//             double dy = pNeighbor->position[1] - pCell->position[1];
//             double dz = 0.; // pNeighbor->position[2] - pCell->position[2];
//             //jdouble distance = std::sqrt( dx*dx + dy*dy + dz*dz );
//             double distance = std::sqrt( dx*dx + dy*dy);

//             if( distance < 1e-12 ) continue;  // avoid division by zero

//             double overlap = distance - rest_length;  // negative when cells overlap

//             if( overlap >= 0.0 )
//             {
//                 continue;

//                 // Attraction branch (linear * exponential decay); only reached if
//                 // distance >= rest_length -- RepulsionForce skips this, but kept
//                 // for completeness to match GeneralisedLinearSpringForce exactly.

//                 // double magnitude = mu * overlap * std::exp( -alpha * overlap / rest_length );
//                 // pCell->velocity[0] += magnitude * dx / distance;
//                 // pCell->velocity[1] += magnitude * dy / distance;
//                 // pCell->velocity[2] += magnitude * dz / distance;
//             }
//             else
//             {
//                 // Repulsion branch: overlap/rest_length in (-1, 0)
//                 // log argument in (0, 1) -> magnitude < 0 -> force points away from B

//                 // logarithmic (matches Chaste's "stable" form)
//                 // double magnitude = mu * std::log( 1.0 + overlap / rest_length );

//                 // quadratic magnitude, preserving sign for direction
//                 double magnitude = m_stiffness * overlap * std::abs(overlap) / rest_length;

//                 pCell->velocity[0] += magnitude * dx / distance;
//                 pCell->velocity[1] += magnitude * dy / distance;
//                 // pCell->velocity[2] += magnitude * dz / distance;
//                 // pCell->velocity[2] = 0.0;
//             }

//             if( false )
//             {

//             // ----------  compute f_i and a_i
//             static double pi = 3.1415926535897932384626433832795; 
//             static double one_div_pi = 0.31830989; 
//             // static double two_pi = 6.283185307179586476925286766559; 

//             if ((*all_cells).size() < 2)
//                 return;

//             // pCell->custom_data["cell_ID"] =  pCell->ID;

//             // double r1 = phenotype.geometry.radius;
//             double ra_2 = r_a*r_a;
//             // double x1 = (*pCell).position[0];
//             // double y1 = (*pCell).position[1];
//             // double gamma = 0.0;
//             // double beta = 0.0;

//             // pCell->custom_data["num_nbrs"] = pCell->state.neighbors.size();   // new: Nov 4 '25
//             pCell->custom_data["num_nbrs"] = 0;

//             // if (d < r1+r2)
//             if (distance < rest_length)
//             {
//                 pCell->custom_data["num_nbrs"] += 1;
//                 // phi_ij = ( d_ij*d_ij - r_j*r_j + r_i*r_i ) / (2 * d_ij * r_i)
//                 double phi = (distance*distance - r_b*r_b + ra_2 ) / (2 * distance * r_a);
//                 // f_i = 1.0 - 1.0/np.pi * np.sqrt( 1.0 - phi_ij*phi_ij )   # for all nbrs j:  1 - 1/pi * SUM_j (sqrt(1 - phi_ij)
//                 gamma += sqrt( 1.0 - phi*phi);
//                 // if (fabs(PhysiCell_globals.current_time - 120.0) < 1.e-4)
//                 //     std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << ", gamma= " << gamma << std::endl;

//                 beta += acos(phi) - phi * sqrt(1 - phi*phi);
//             }

//             double gamma_inv = one_div_pi * gamma;     // 1.0/pi * gamma;
//             gamma = 1.0 - gamma_inv;   // free surface fraction
//             if (gamma < 0.0)  gamma = 0.0;

//             pCell->custom_data["f_i"] = gamma;
//             // pCell->custom_data["gamma_inv"] = gamma_inv;

//             beta = 1.0 - 1.0/pi * beta;
//             // pCell->custom_data["beta"] = 1.0 - 1.0/pi * beta;
//             if (beta < 0.0)  beta = 0.0;
//             pCell->custom_data["a_i"] = beta;
//             }
//         }
//     }
// }

// do every mechanics dt
void custom_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
    // std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << std::endl;

    static double pi = 3.1415926535897932384626433832795; 
    static double one_div_pi = 0.31830989; 
	// static double two_pi = 6.283185307179586476925286766559; 

    if ((*all_cells).size() < 2)
        return;

    pCell->custom_data["cell_ID"] =  pCell->ID;

    double r_a = phenotype.geometry.radius;

    double r_a_2 = r_a * r_a;
    double x1 = (*pCell).position[0];
    double y1 = (*pCell).position[1];
    double gamma = 0.0;
    double beta = 0.0;

    // pCell->custom_data["num_nbrs"] = pCell->state.neighbors.size();   // new: Nov 4 '25
    pCell->custom_data["num_nbrs"] = 0;

    // std::vector<int> my_nbrs;
	// for( int idx=0; idx<pCell->state.neighbors.size(); idx++ )  // for all j nbrs
	// {
    // for( auto pNeighbor : *PhysiCell::all_cells )
    // for( auto pNeighbor : *PhysiCell::all_cells )  // rwh: optimize!
    for (const auto& pNeighbor : pCell->state.neighbors)
    {
        if( pNeighbor == pCell ) continue;

        double r_b = pNeighbor->phenotype.geometry.radius;
        double rest_length = r_a + r_b;

        // if(std::abs(pCell->position[0]-pNeighbor->position[0]) > rest_length) continue;
        // else if(std::abs(pCell->position[1]-pNeighbor->position[1]) > rest_length) continue;

        double dx = pNeighbor->position[0] - pCell->position[0];
        double dy = pNeighbor->position[1] - pCell->position[1];
        double dz = 0.; // pNeighbor->position[2] - pCell->position[2];
        //jdouble distance = std::sqrt( dx*dx + dy*dy + dz*dz );
        double distance = std::sqrt( dx*dx + dy*dy);

		// Cell* pC = pCell->state.neighbors[idx]; 

        // // check for buggy duplicate nbr and skip if present
        // std::vector<int>::iterator it = std::find(my_nbrs.begin(), my_nbrs.end(), pC->ID);
        // if (it != my_nbrs.end())
        //     continue;
        // my_nbrs.push_back(pC->ID);

        // compute chord of intersection (if any)
        // radii of cells
        // double r2 = pC->phenotype.geometry.radius;
        // centers of cells
        // double x2 = (*pC).position[0];
        // double y2 = (*pC).position[1];
        // double xdiff = x1-x2;
        // double ydiff = y1-y2;
        // double d = sqrt(xdiff*xdiff + ydiff*ydiff);
        // if (d < r1+r2)
        if (distance < rest_length)
        {
            pCell->custom_data["num_nbrs"] += 1;
            // phi_ij = ( d_ij*d_ij - r_j*r_j + r_i*r_i ) / (2 * d_ij * r_i)
            // double phi = (d*d - r2*r2 + r1_2 ) / (2 * d * r1);
            double phi = (distance*distance - r_b*r_b + r_a_2 ) / (2 * distance * r_a);
            // f_i = 1.0 - 1.0/np.pi * np.sqrt( 1.0 - phi_ij*phi_ij )   # for all nbrs j:  1 - 1/pi * SUM_j (sqrt(1 - phi_ij)
            gamma += sqrt( 1.0 - phi*phi);
            // if (fabs(PhysiCell_globals.current_time - 120.0) < 1.e-4)
            //     std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << ", gamma= " << gamma << std::endl;

            beta += acos(phi) - phi * sqrt(1 - phi*phi);
        }
    }

    double gamma_inv = one_div_pi * gamma;     // 1.0/pi * gamma;
    gamma = 1.0 - gamma_inv;   // free surface fraction
    if (gamma < 0.0)  gamma = 0.0;

    pCell->custom_data["f_i"] = gamma;
    // pCell->custom_data["gamma_inv"] = gamma_inv;

    beta = 1.0 - 1.0/pi * beta;
    // pCell->custom_data["beta"] = 1.0 - 1.0/pi * beta;
    if (beta < 0.0)  beta = 0.0;
    pCell->custom_data["a_i"] = beta;

    {
        // ------rwh: don't need to do this for the 1000 cell, no contact inhibition model!
        pCell->custom_data["arrest_cycle"] =  0.0;  // rwh: improve?
        pCell->custom_data["beta_or_gamma"] = 0;
        if (beta < beta_threshold)
        {
            pCell->custom_data["arrest_cycle"] =  1.0;  // rwh: improve?
            pCell->custom_data["beta_or_gamma"] += 1;
        }
        if (gamma < gamma_threshold)
        {
            pCell->custom_data["arrest_cycle"] =  1.0;  // rwh: improve?
            pCell->custom_data["beta_or_gamma"] += 2;
        }
        return;
    }

    return; 
}

// do every mechanics dt
// void custom_function_v0( Cell* pCell, Phenotype& phenotype, double dt )
// {
//     // std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << std::endl;

//     static double pi = 3.1415926535897932384626433832795; 
//     static double one_div_pi = 0.31830989; 
// 	// static double two_pi = 6.283185307179586476925286766559; 

//     if ((*all_cells).size() < 2)
//         return;

//     pCell->custom_data["cell_ID"] =  pCell->ID;

//     // x_j = x0[idx]
//     // c_i = Circle((x_i, y_i), r_i, edgecolor='r', facecolor='none')
//     // c_j = Circle((x_j, y_j), r_j, edgecolor='r', facecolor='none')
//     // ax[idx].add_patch(c_i)
//     // ax[idx].add_patch(c_j)
//     // ax[idx].set_xlim(-1.5, 2.5)
//     // ax[idx].set_ylim(-1.5, 1.5)
//     // ax[idx].set_aspect('equal', adjustable='box')
//     // #plt.show()
    
//     // xd = x_j-x_i
//     // yd = y_j-y_i
//     // d_ij = np.sqrt(xd*xd + yd*yd)
//     // # print("d_ij = ",d_ij)
//     // #  *d_ij^2 - r_j^2 + r_i^2) / (2 * d_ij * r_i)
//     // phi_ij = ( d_ij*d_ij - r_j*r_j + r_i*r_i ) / (2 * d_ij * r_i)
//     // # print("phi_ij = ",phi_ij)
//     // # free surface fraction of cell_i  (gamma)
//     // f_i = 1.0 - 1.0/np.pi * np.sqrt( 1.0 - phi_ij*phi_ij )   # for all nbrs j:  1 - 1/pi * SUM_j (sqrt(1 - phi_ij)
//     // # print("f_i = ",f_i)
//     // ax[idx].text(-0.9, 1.2, f'f_i={f_i:.3f}', fontsize = 10)
    
//     // # A_i/A_i0 = 1 - 1/pi SUM_j ( arccos(phi_ij) - phi_ij * sqrt(1 - phi_ij^2) )    # beta
//     // relative_compressed_size = 1.0 - 1.0/np.pi * (np.arccos(phi_ij) - phi_ij * np.sqrt(1.0 - phi_ij*phi_ij))

//     double r1 = phenotype.geometry.radius;
//     double r1_2 = r1*r1;
//     double x1 = (*pCell).position[0];
//     double y1 = (*pCell).position[1];
//     double gamma = 0.0;
//     double beta = 0.0;

// 	// if( pCell->state.neighbors.size() == 0)  // for all j nbrs
//     // {
//     //     pCell->custom_data["gamma"] = 0.0;
//     //     pCell->custom_data["beta"] = 0.0;
//     //     return;
//     // }

//     // pCell->custom_data["num_nbrs"] = pCell->state.neighbors.size();   // new: Nov 4 '25
//     pCell->custom_data["num_nbrs"] = 42;  // rwh - fix

//     // std::vector<int> my_nbrs;
// 	for( int idx=0; idx<pCell->state.neighbors.size(); idx++ )  // for all j nbrs
// 	{
// 		Cell* pC = pCell->state.neighbors[idx]; 

//         // check for buggy duplicate nbr and skip if present
//         // std::vector<int>::iterator it = std::find(my_nbrs.begin(), my_nbrs.end(), pC->ID);
//         // if (it != my_nbrs.end())
//         //     continue;
//         // my_nbrs.push_back(pC->ID);

//         // if (fabs(PhysiCell_globals.current_time - 120.0) < 1.e-4)
//         // {
//         //     // std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << ", gamma= " << gamma << std::endl;
//         //     std::cout << __FUNCTION__ << " -------- t= " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << " has nbr ID= " << pC->ID << std::endl;

//         // }

//         // compute chord of intersection (if any)
//         // radii of cells
//         double r2 = pC->phenotype.geometry.radius;
//         // centers of cells
//         double x2 = (*pC).position[0];
//         double y2 = (*pC).position[1];
//         double xdiff = x1-x2;
//         double ydiff = y1-y2;
//         double d = sqrt(xdiff*xdiff + ydiff*ydiff);
//         if (d < r1+r2)
//         {
//             // std::cout << "cell " << pCell->ID << " intersects cell " << pC->ID << std::endl;
//             // std::cout << "x1,y1 " << x1 << ", " << y1 << std::endl;
//             // std::cout << "x2,y2 " << x2 << ", " << y2 << std::endl;
//             // std::cout << "  r1,r2 " << r1 << ", " << r2 << std::endl;
//             // std::cout << "cell " << pCell->ID << " intersects cell " << pC->ID << std::endl;

//     // phi_ij = ( d_ij*d_ij - r_j*r_j + r_i*r_i ) / (2 * d_ij * r_i)
//     // # print("phi_ij = ",phi_ij)
//     // # free surface fraction of cell_i  (gamma)
//     // f_i = 1.0 - 1.0/np.pi * np.sqrt( 1.0 - phi_ij*phi_ij )   # for all nbrs j:  1 - 1/pi * SUM_j (sqrt(1 - phi_ij)
//     // # print("f_i = ",f_i)
//     // ax[idx].text(-0.9, 1.2, f'f_i={f_i:.3f}', fontsize = 10)
    
//     // # A_i/A_i0 = 1 - 1/pi * SUM_j ( arccos(phi_ij) - phi_ij * sqrt(1 - phi_ij^2) )    # beta
//     // relative_compressed_size = 1.0 - 1.0/np.pi * (np.arccos(phi_ij) - phi_ij * np.sqrt(1.0 - phi_ij*phi_ij))

//             // phi_ij = ( d_ij*d_ij - r_j*r_j + r_i*r_i ) / (2 * d_ij * r_i)
//             double phi = (d*d - r2*r2 + r1_2 ) / (2 * d * r1);
//             // f_i = 1.0 - 1.0/np.pi * np.sqrt( 1.0 - phi_ij*phi_ij )   # for all nbrs j:  1 - 1/pi * SUM_j (sqrt(1 - phi_ij)
//             gamma += sqrt( 1.0 - phi*phi);
//             // if (fabs(PhysiCell_globals.current_time - 120.0) < 1.e-4)
//             //     std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << ", gamma= " << gamma << std::endl;

//             beta += acos(phi) - phi * sqrt(1 - phi*phi);
//         }
//     }
//     // if (fabs(PhysiCell_globals.current_time - 120.0) < 1.e-4)
//     // {
//     //     std::cout << __FUNCTION__ << " ------- t= " << PhysiCell_globals.current_time << ": gamma= " << gamma << std::endl;
//     // }
//     double gamma_inv = one_div_pi * gamma;     // 1.0/pi * gamma;   
//     gamma = 1.0 - gamma_inv;   // free surface fraction
//     if (gamma < 0.0)  gamma = 0.0;

//     pCell->custom_data["f_i"] = gamma;
//     // pCell->custom_data["gamma_inv"] = gamma_inv;

//     beta = 1.0 - 1.0/pi * beta;
//     // pCell->custom_data["beta"] = 1.0 - 1.0/pi * beta;
//     if (beta < 0.0)  beta = 0.0;
//     pCell->custom_data["a_i"] = beta;

//     // beta ~[0,1];  gamma~[0,1] (?)
//     // if ((gamma < 0.025) && (beta > 0.925) )
//     // {
//     //     set_single_behavior( pCell , "cycle entry" , 0.0); 
//     // }

//     // std::cout << "-------- gamma= " << gamma << ",  beta= " << beta << std::endl;

//     // Decide which value of the 3x3 matrix of param ranges we're running/saving
//     // Inhibition is OR; Growth is AND
// 	// if (param_matrix_elm == 3)
//     // {
//     //     set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
//     //     if ((beta > 0.9) && (gamma > 0.5))  // allow cell cycle/growth/prolif
//     //     {
//     //         // set cycle rate = 0.002257  (duration=443)
//     //         // set_single_behavior( pCell , "cycle entry" , 0.002257); 
//     //         set_single_behavior( pCell , "cycle entry" , 0.01124);  // 1/89
//     //     }
//     //     return;
//     // }
// 	// else if (param_matrix_elm == 0)
//     {
//         // set_single_behavior( pCell , "cycle entry" , 0.01124);  // 1/89
//         // set_single_behavior( pCell , "growth_rate" , 5.883);  // 1/89
//         // pCell->custom_data["growth_rate"] = 5.883;

//         // pCell->custom_data["growth_rate"] = parameters.doubles("growth_rate");
//         pCell->custom_data["arrest_cycle"] =  0.0;  // rwh: improve?

//         pCell->custom_data["beta_or_gamma"] = 0;
//         // if ((beta < 0.9) || (gamma < 0.9))
//         // if ((beta < 0.5) || (gamma < 0.9))
//         // if ((beta < beta_threshold) || (gamma < gamma_threshold))
//         // {
//         //     set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
//         // }
//         if (beta < beta_threshold)
//         {
//             // set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
//             // set_single_behavior( pCell , "growth_rate" , 0.0);
//             // pCell->custom_data["growth_rate"] = 0.0;
//             pCell->custom_data["arrest_cycle"] =  1.0;  // rwh: improve?
//             pCell->custom_data["beta_or_gamma"] += 1;
//         }
//         if (gamma < gamma_threshold)
//         {
//             // set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
//             // set_single_behavior( pCell , "growth_rate" , 0.0);
//             // pCell->custom_data["growth_rate"] = 0.0;
//             pCell->custom_data["arrest_cycle"] =  1.0;  // rwh: improve?
//             pCell->custom_data["beta_or_gamma"] += 2;
//         }
//         return;
//     }

//     //-------------
// 	// else if (param_matrix_elm == 0)
//     // {
//     //     // set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
//     //     set_single_behavior( pCell , "cycle entry" , 0.01124);  // 1/89
//     //     // set_single_behavior( pCell , "cycle entry" , 0.01124);  // allow growth, by default (=1/89)
//     //     // if ((beta < 0.8) && (gamma < 0.8))  // arrest growth
//     //     // if ((beta < 0.8) && (gamma < 0.8))  // arrest growth
//     //     // if ((beta > 0.6) && (gamma > 0.6))
//     //     if ((beta < 0.9) || (gamma < 0.9))
//     //     {
//     //         // set_single_behavior( pCell , "cycle entry" , 0.0);  
//     //         // set_single_behavior( pCell , "cycle entry" , 0.01124);  // 1/89
//     //         set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
//     //     }
//     //     return;
//     // }
//     return; 
// }

// do every phenotype dt
// No longer critical assumption that dt_phenotype = 1
void custom_volume_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	static double growth_rate = parameters.doubles("cell_area_0") / parameters.doubles("cycle_duration");   // e.g., 78.54/(88.3 * 5=441.5) = 0.1779


    if (pCell->custom_data["arrest_cycle"] >  0.0)  // rwh: improve?
    {
        return;
    }

    // if (PhysiCell_globals.current_time < 5)
    //     std::cout << " --- "<<__FUNCTION__ << ":  t= " << PhysiCell_globals.current_time << ", growth_rate= " << growth_rate << std::endl;


    // option 1: Absolute (time-anchored)
	// static double cycle_duration = parameters.doubles("cycle_duration");  
    // pCell->custom_data["time_in_cycle"] +=  dt;
    // pCell->custom_data["cell_area"] =  pCell->custom_data["cell_area_0"] * (1.0 + pCell->custom_data["time_in_cycle"] / cycle_duration );  

    // option 2:  Incremental (stateless delta
    // pCell->custom_data["cell_area"] += pCell->custom_data["growth_rate"] * dt;  
    // pCell->custom_data["cell_area"] += growth_rate * dt;  

    // option 3: 
    pCell->custom_data["time_in_cycle"] +=  dt;
    // pCell->custom_data["cell_area"] =  pCell->custom_data["cell_area_0"] +  (0.169 * pCell->custom_data["time_in_cycle"] );
    pCell->custom_data["cell_area"] =  pCell->custom_data["cell_area_0"] +  (growth_rate * pCell->custom_data["time_in_cycle"] );

    double radius = std::sqrt(pCell->custom_data["cell_area"] / 3.14159);  // A=pi * r^2 ; r= sqrt(A/pi)

    static double four_thirds_pi =  4.188790204786391;
	double volume = four_thirds_pi * radius * radius * radius;   // Volume(sphere) = (4/3) * pi * r^3
    pCell->set_total_volume( volume );

    double area_norm_rand = pCell->custom_data["norm_rand"] * target_area; // uh, only changes when norm_rand changes, so at division

    pCell->custom_data["cell_radius"] = radius;   // r of cell_area

    if (pCell->custom_data["cell_area"] >= area_norm_rand)  
    {
        #pragma omp critical
        {
            pCell->divide();
        }
    }
}

// just for debugging daughter cells separation
std::vector<double> custom_division_dir_function( Cell* pParent)
{
    // const std::vector<double> output = {cos(theta),sin(theta),0};
    const std::vector<double> output = {1,0,0};
    // return normalize(output);
    return(output);
}
// Assign N(2,0.4^2) to each daughter cell. Exit the sim at 10K cells (or "max_cells")
void custom_division_function( Cell* pCell1, Cell* pCell2 )
{ 
    static bool normal_rand_flag = parameters.bools("normal_random_flag");   // 5, 10, ?

    // if (pCell1->ID == 0)
    // {
    //     std::cout << " --- "<<__FUNCTION__ << ":  t= " << PhysiCell_globals.current_time << ", cell1 ID= " << pCell1->ID << std::endl;
    // }
    // static int monolayer_max_cells = 9999;
    static int monolayer_max_cells = 10000;
    static double draw_mean = 2.0;
    static double draw_stddev = 0.4;

    pCell1->custom_data["time_in_cycle"] =  0.0;
    pCell2->custom_data["time_in_cycle"] =  0.0;


    // std::cout << " ---- "<<__FUNCTION__ << ":  t= " << PhysiCell_globals.current_time << std::endl;
    // std::cout << " ---- "<<__FUNCTION__ << ":  pCell1 cell_radius= " << pCell1->custom_data["cell_radius"] <<  ", cell_area= " << pCell1->custom_data["cell_area"] << std::endl;
    // std::cout << "     pre: pCell1 cell_area= "<<pCell1->custom_data["cell_area"] << std::endl;
    // std::cout << "     pre: pCell1 cell_radius * radius * pi = "<<pCell1->custom_data["cell_area"]*pCell1->custom_data["cell_area"] * 3.1416 << std::endl;

    // pCell1->custom_data["cell_radius"] = pCell1->phenotype.geometry.radius;
    // pCell2->custom_data["cell_radius"] = pCell2->phenotype.geometry.radius;
    // std::cout << "     post: \n";
    pCell1->custom_data["cell_area"] /= 2.0;
    pCell2->custom_data["cell_area"] /= 2.0;

    // pCell1->custom_data["growth_rate"] = pCell1->custom_data["cell_area"] / 441.5;  // 88.3 * 5;
    // pCell2->custom_data["growth_rate"] = pCell2->custom_data["cell_area"] / 441.5;  // 88.3 * 5;

    pCell1->custom_data["cell_area_0"] =  pCell1->custom_data["cell_area"];
    pCell2->custom_data["cell_area_0"] =  pCell2->custom_data["cell_area"];

    double radius = sqrt(pCell1->custom_data["cell_area"] / 3.14159);  // A=pi * r^2 ; r= sqrt(A/pi)
    pCell1->custom_data["cell_radius"] = radius;
    pCell2->custom_data["cell_radius"] = radius;

    // std::cout << "     post: pCell1 cell_radius= "<<pCell1->custom_data["cell_radius"] << std::endl;
    // std::cout <<       ":  pCell1 cell_radius= " << pCell1->custom_data["cell_radius"] <<  ", cell_area= " << pCell1->custom_data["cell_area"] << std::endl;

    // update volumes to match desired areas
    static double four_thirds_pi =  4.188790204786391;
	double volume = four_thirds_pi * radius * radius * radius;   // Volume(sphere) = (4/3) * pi * r^3
    pCell1->set_total_volume( volume );
    pCell2->set_total_volume( volume );
    // std::cout <<       ":  pCell1 volume= " << volume <<  std::endl;

    bool custom_daughters_pos = false;
    custom_daughters_pos = true;
    if (custom_daughters_pos)
    {

    // Need to modify the center of the new daughter cell to match new areas/volumes
	double x0 = pCell1->position[0];
	double y0 = pCell1->position[1];
    // std::cout <<       ":  pCell1 x,y= " << x0 << ", "<<y0 <<  std::endl;
	double x1 = pCell2->position[0];
	double y1 = pCell2->position[1];
    // std::cout <<       ":  pCell2 x,y= " << x1 << ", "<<y1 <<  std::endl;
	double xnew = x1 - x0;
	double ynew = y1 - y0;
    std::vector<double> vec = {xnew,ynew};
    normalize(&vec);
    // std::cout <<       ":  vec= " << vec[0] << ", "<<vec[1] <<  std::endl;
    const double overlap_factor = 1.7;
	// xnew = x0 + 2*radius * vec[0];
	// ynew = y0 + 2*radius * vec[1];
	xnew = x0 + overlap_factor*radius * vec[0];
	ynew = y0 + overlap_factor*radius * vec[1];
	pCell2->assign_position(xnew,ynew,0.0);
    }


    // live.phases[0].division_at_phase_exit = true; 
    // std::cout << "       pCell1...division_at_phase_exit  = " << pCell1->phenotype.cycle.pCycle_Model->phases[0].division_at_phase_exit  << std::endl;
    // std::cout << "       pCell2...division_at_phase_exit  = " << pCell2->phenotype.cycle.pCycle_Model->phases[0].division_at_phase_exit  << std::endl;

    // This is crucial! Otherwise, "division_at_phase_exit" is true when a cell divides.
    pCell1->phenotype.cycle.pCycle_Model->phases[0].division_at_phase_exit = false;
    pCell2->phenotype.cycle.pCycle_Model->phases[0].division_at_phase_exit = false;

    if (parameters.ints.find_index("max_cells") != -1)
	{
        monolayer_max_cells = parameters.ints("max_cells");
	}

    // static int idx_default = find_cell_definition_index("default");
    // static int idx_ctype1 = find_cell_definition_index("ctype1");
    // std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell IDs= " << pCell1->ID << ", " << pCell2->ID << std::endl;

    // // Asymmetric division
    // if (UniformRandom() < 0.5)
    // {
    //     pCell2->convert_to_cell_definition( *cell_definitions_by_index[idx_default] ); 
    // }
    // else
    // {
    //     pCell2->convert_to_cell_definition( *cell_definitions_by_index[idx_ctype1] ); 
    // }

	char filename[1024];
    static std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;
    static std::string (*substrate_coloring_function)(double, double, double) = paint_by_density_percentage;

    int ncells = (*all_cells).size();
    if ( ncells >= monolayer_max_cells )
    {
        std::cout << "-------- # cells: " << ncells << std::endl;
        sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
        // sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
        save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
        
        sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index );
        // sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() );
        SVG_plot(filename, microenvironment, 0.0, PhysiCell_globals.current_time, cell_coloring_function, substrate_coloring_function);

        // timer 
        std::cout << std::endl << "Total simulation runtime: " << std::endl; 
        BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

        std::cout << std::endl; 
        exit(-1);
    }

    if (normal_rand_flag)
    {

    // set each daughter cell's "X" value to determine cell division: A(t) = X * A_0(0)
    double draw = NormalRandom(draw_mean, draw_stddev);  //     // where: NormalRandom(mean, std_dev)
    while (draw < 0.0)
    {
        draw = NormalRandom(draw_mean, draw_stddev); 
    }
    // pCell1->custom_data["growth_rate"]
    pCell1->custom_data["norm_rand"] = draw;
    // pCell1->custom_data["norm_rand"] = 3.0;  //rwh - temporary


    draw = NormalRandom(draw_mean, draw_stddev);  //     // where: NormalRandom(mean, std_dev)
    while (draw < 0.0)
    {
        draw = NormalRandom(draw_mean, draw_stddev); 
    }
    pCell2->custom_data["norm_rand"] = draw;
    // pCell2->custom_data["norm_rand"] = 1.5;  //rwh - temporary
    }
    else
    {
        pCell2->custom_data["norm_rand"] = 2.0;
    }

    return; 
}