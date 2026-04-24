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

constexpr double pi = 3.1415926535897932384626433832795; 
constexpr double one_div_pi = 0.31830989; 

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

	// cell_defaults.functions.custom_cell_rule = custom_function;   // do every mechanics dt
	// cell_defaults.functions.custom_cell_rule = custom_cell_rule;   // do every mechanics dt

	// cell_defaults.functions.contact_function = contact_function; 

	cell_defaults.functions.update_velocity = custom_update_cell_velocity;    // dt_mech

    // rwh: simple test - hard-coded division direction
    // cell_defaults.functions.cell_division_direction_function = custom_division_dir_function; 

    cell_defaults.functions.cell_division_function = custom_division_function; 

    // cell_defaults.functions.volume_update_function = double_volume_update_function; 
    cell_defaults.functions.volume_update_function = custom_volume_function;   // dt_phenotype

    // target_volume = cell_defaults.phenotype.volume.total;  // 523.6
    target_area = pi * cell_defaults.phenotype.geometry.radius * cell_defaults.phenotype.geometry.radius;

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


// do every mechanics dt
void custom_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{
    static double m_stiffness = parameters.doubles("m_stiffness");   // 5, 10, ?

    pCell->velocity = {0.0, 0.0, 0.0};
    // pCell->velocity[2] = 0.0;    // shouldn't be necessary for this model

    double r1 = phenotype.geometry.radius;
    for( auto pC2 : *PhysiCell::all_cells )   // rwh: optimize
    {
        if( pC2 == pCell ) continue;

        double r2 = pC2->phenotype.geometry.radius;
        double touching_dist = r1 + r2;

        if(std::abs(pCell->position[0]-pC2->position[0]) > touching_dist) continue;
        else if(std::abs(pCell->position[1]-pC2->position[1]) > touching_dist) continue;

        // Vector from A (pCell) toward B (pC2)
        double dx = pC2->position[0] - pCell->position[0];
        double dy = pC2->position[1] - pCell->position[1];
        double distance = std::sqrt( dx*dx + dy*dy);

        if( distance < 1e-12 ) continue;  // avoid division by zero

        double overlap = distance - touching_dist;  // negative when cells overlap
        if( overlap >= 0.0 )
        {
            continue;   // we have no cell-cell attraction in this model
        }
        else
        {
            // quadratic magnitude, preserving sign for direction
            double magnitude = m_stiffness * overlap * std::abs(overlap) / touching_dist;
            pCell->velocity[0] += magnitude * dx / distance;
            pCell->velocity[1] += magnitude * dy / distance;
        }
    }
}

// do every phenotype dt
// Seem to get better gamma, beta (no inhibition) 1000-cell results with dt_phenotype = 1 (or lower)
void custom_volume_function( Cell* pCell, Phenotype& phenotype, double dt )
{
    constexpr double four_thirds_pi =  4.188790204786391;

    if (pCell->custom_data["arrest_cycle"] >  0.0)  // rwh: improve?
    {
        return;
    }

	static double cycle_duration = parameters.doubles("cycle_duration");  

    pCell->custom_data["time_in_cycle"] +=  dt;

    pCell->custom_data["cell_area"] =  pCell->custom_data["cell_area_0"] * (1.0 + pCell->custom_data["time_in_cycle"] / cycle_duration );  

    double radius = sqrt(pCell->custom_data["cell_area"] / pi);  // A=pi * r^2 ; r= sqrt(A/pi)

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
    else // Should we compute this cell's gamma,beta values only if this cell didn't divide, or regardless??
    {
        double r1 = pCell->phenotype.geometry.radius;
        double r1_squared = radius*radius;
        double beta = 0.0;
        double gamma = 0.0;
        pCell->custom_data["num_nbrs"] = 0;
        for( auto pC2 : *PhysiCell::all_cells )  
        {
            if( pCell == pC2 ) continue;

            double r2 = pC2->phenotype.geometry.radius;
            double touching_dist = radius + r2;

            if(std::abs(pCell->position[0]-pC2->position[0]) > touching_dist) continue;
            else if(std::abs(pCell->position[1]-pC2->position[1]) > touching_dist) continue;

            double dx = pC2->position[0] - pCell->position[0];
            double dy = pC2->position[1] - pCell->position[1];
            // double dz = 0.; // pNeighbor->position[2] - pCell->position[2];
            //jdouble distance = std::sqrt( dx*dx + dy*dy + dz*dz );
            double distance = std::sqrt( dx*dx + dy*dy);
            if (distance < touching_dist)
            {
                pCell->custom_data["num_nbrs"] += 1;
                // phi_ij = ( d_ij*d_ij - r_j*r_j + r_i*r_i ) / (2 * d_ij * r_i)
                // double phi = (d*d - r2*r2 + r1_2 ) / (2 * d * r1);
                double phi = (distance*distance - r2*r2 + r1_squared ) / (2 * distance * radius);
                // f_i = 1.0 - 1.0/np.pi * np.sqrt( 1.0 - phi_ij*phi_ij )   # for all nbrs j:  1 - 1/pi * SUM_j (sqrt(1 - phi_ij)
                gamma += std::sqrt( 1.0 - phi*phi);
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


        // ------Note: don't need to do this for the 1000 cell, no contact inhibition model!
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

    // Let's match cell division behavior of most other frameworks in OpenVT
    // bool custom_daughters_pos = false;
    // custom_daughters_pos = true;
    // if (custom_daughters_pos)
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

    return; 
}