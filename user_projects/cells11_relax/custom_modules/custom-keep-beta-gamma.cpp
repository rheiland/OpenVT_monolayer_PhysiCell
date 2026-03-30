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

// int param_matrix_elm; // 0,1,2 (top row); 3,4,5 (middle); 6,7,8 (bottom)
double beta_threshold, gamma_threshold;

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
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

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
	cell_defaults.functions.custom_cell_rule = custom_function; 
	// cell_defaults.functions.contact_function = contact_function; 
    cell_defaults.functions.cell_division_function = custom_division_function; 
    // cell_defaults.functions.volume_update_function = double_volume_update_function; 
    cell_defaults.functions.volume_update_function = custom_volume_function; 
	
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
void custom_function( Cell* pCell, Phenotype& phenotype, double dt )
{ 
    // std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << std::endl;

    static double pi = 3.1415926535897932384626433832795; 
	// static double two_pi = 6.283185307179586476925286766559; 

    if ((*all_cells).size() < 2)
        return;

    // x_j = x0[idx]
    // c_i = Circle((x_i, y_i), r_i, edgecolor='r', facecolor='none')
    // c_j = Circle((x_j, y_j), r_j, edgecolor='r', facecolor='none')
    // ax[idx].add_patch(c_i)
    // ax[idx].add_patch(c_j)
    // ax[idx].set_xlim(-1.5, 2.5)
    // ax[idx].set_ylim(-1.5, 1.5)
    // ax[idx].set_aspect('equal', adjustable='box')
    // #plt.show()
    
    // xd = x_j-x_i
    // yd = y_j-y_i
    // d_ij = np.sqrt(xd*xd + yd*yd)
    // # print("d_ij = ",d_ij)
    // #  *d_ij^2 - r_j^2 + r_i^2) / (2 * d_ij * r_i)
    // phi_ij = ( d_ij*d_ij - r_j*r_j + r_i*r_i ) / (2 * d_ij * r_i)
    // # print("phi_ij = ",phi_ij)
    // # free surface fraction of cell_i  (gamma)
    // f_i = 1.0 - 1.0/np.pi * np.sqrt( 1.0 - phi_ij*phi_ij )   # for all nbrs j:  1 - 1/pi * SUM_j (sqrt(1 - phi_ij)
    // # print("f_i = ",f_i)
    // ax[idx].text(-0.9, 1.2, f'f_i={f_i:.3f}', fontsize = 10)
    
    // # A_i/A_i0 = 1 - 1/pi SUM_j ( arccos(phi_ij) - phi_ij * sqrt(1 - phi_ij^2) )    # beta
    // relative_compressed_size = 1.0 - 1.0/np.pi * (np.arccos(phi_ij) - phi_ij * np.sqrt(1.0 - phi_ij*phi_ij))

    double r1 = phenotype.geometry.radius;
    double r1_2 = r1*r1;
    double x1 = (*pCell).position[0];
    double y1 = (*pCell).position[1];
    double gamma = 0.0;
    double beta = 0.0;

	// if( pCell->state.neighbors.size() == 0)  // for all j nbrs
    // {
    //     pCell->custom_data["gamma"] = 0.0;
    //     pCell->custom_data["beta"] = 0.0;
    //     return;
    // }

    pCell->custom_data["num_nbrs"] = pCell->state.neighbors.size();   // new: Nov 4 '25
	for( int idx=0; idx<pCell->state.neighbors.size(); idx++ )  // for all j nbrs
	{
		Cell* pC = pCell->state.neighbors[idx]; 

        // compute chord of intersection (if any)
        // radii of cells
        double r2 = pC->phenotype.geometry.radius;
        // centers of cells
        double x2 = (*pC).position[0];
        double y2 = (*pC).position[1];
        double xdiff = x1-x2;
        double ydiff = y1-y2;
        double d = sqrt(xdiff*xdiff + ydiff*ydiff);
        if (d < r1+r2)
        {
            // std::cout << "cell " << pCell->ID << " intersects cell " << pC->ID << std::endl;
            // std::cout << "x1,y1 " << x1 << ", " << y1 << std::endl;
            // std::cout << "x2,y2 " << x2 << ", " << y2 << std::endl;
            // std::cout << "  r1,r2 " << r1 << ", " << r2 << std::endl;
            // std::cout << "cell " << pCell->ID << " intersects cell " << pC->ID << std::endl;

    // phi_ij = ( d_ij*d_ij - r_j*r_j + r_i*r_i ) / (2 * d_ij * r_i)
    // # print("phi_ij = ",phi_ij)
    // # free surface fraction of cell_i  (gamma)
    // f_i = 1.0 - 1.0/np.pi * np.sqrt( 1.0 - phi_ij*phi_ij )   # for all nbrs j:  1 - 1/pi * SUM_j (sqrt(1 - phi_ij)
    // # print("f_i = ",f_i)
    // ax[idx].text(-0.9, 1.2, f'f_i={f_i:.3f}', fontsize = 10)
    
    // # A_i/A_i0 = 1 - 1/pi * SUM_j ( arccos(phi_ij) - phi_ij * sqrt(1 - phi_ij^2) )    # beta
    // relative_compressed_size = 1.0 - 1.0/np.pi * (np.arccos(phi_ij) - phi_ij * np.sqrt(1.0 - phi_ij*phi_ij))

            // phi_ij = ( d_ij*d_ij - r_j*r_j + r_i*r_i ) / (2 * d_ij * r_i)
            double phi = (d*d - r2*r2 + r1_2 ) / (2 * d * r1);
            // f_i = 1.0 - 1.0/np.pi * np.sqrt( 1.0 - phi_ij*phi_ij )   # for all nbrs j:  1 - 1/pi * SUM_j (sqrt(1 - phi_ij)
            gamma += sqrt( 1.0 - phi*phi);
            // std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << ", gamma= " << gamma << std::endl;
            beta += acos(phi) - phi * sqrt(1 - phi*phi);

        }
    }
    double gamma_inv = 1.0/pi * gamma;   
    gamma = 1.0 - gamma_inv;   // free surface fraction
    // pCell->custom_data["gamma"] = 1.0 - 1.0/pi * gamma;
    pCell->custom_data["f_i"] = gamma;
    // pCell->custom_data["gamma_inv"] = gamma_inv;

    beta = 1.0 - 1.0/pi * beta;
    // pCell->custom_data["beta"] = 1.0 - 1.0/pi * beta;
    pCell->custom_data["a_i"] = beta;

    // beta ~[0,1];  gamma~[0,1] (?)
    // if ((gamma < 0.025) && (beta > 0.925) )
    // {
    //     set_single_behavior( pCell , "cycle entry" , 0.0); 
    // }

    // std::cout << "-------- gamma= " << gamma << ",  beta= " << beta << std::endl;

    // Decide which value of the 3x3 matrix of param ranges we're running/saving
    // Inhibition is OR; Growth is AND
	// if (param_matrix_elm == 3)
    // {
    //     set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
    //     if ((beta > 0.9) && (gamma > 0.5))  // allow cell cycle/growth/prolif
    //     {
    //         // set cycle rate = 0.002257  (duration=443)
    //         // set_single_behavior( pCell , "cycle entry" , 0.002257); 
    //         set_single_behavior( pCell , "cycle entry" , 0.01124);  // 1/89
    //     }
    //     return;
    // }
	// else if (param_matrix_elm == 0)
    {
        // set_single_behavior( pCell , "cycle entry" , 0.01124);  // 1/89
        // set_single_behavior( pCell , "growth_rate" , 5.883);  // 1/89
        pCell->custom_data["growth_rate"] = 5.883;
        pCell->custom_data["beta_or_gamma"] = 0;
        // if ((beta < 0.9) || (gamma < 0.9))
        // if ((beta < 0.5) || (gamma < 0.9))
        // if ((beta < beta_threshold) || (gamma < gamma_threshold))
        // {
        //     set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
        // }
        if (beta < beta_threshold)
        {
            // set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
            // set_single_behavior( pCell , "growth_rate" , 0.0);
            pCell->custom_data["growth_rate"] = 0.0;
            pCell->custom_data["beta_or_gamma"] += 1;
        }
        if (gamma < gamma_threshold)
        {
            // set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
            // set_single_behavior( pCell , "growth_rate" , 0.0);
            pCell->custom_data["growth_rate"] = 0.0;
            pCell->custom_data["beta_or_gamma"] += 2;
        }
        return;
    }

    //-------------
	// else if (param_matrix_elm == 0)
    // {
    //     // set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
    //     set_single_behavior( pCell , "cycle entry" , 0.01124);  // 1/89
    //     // set_single_behavior( pCell , "cycle entry" , 0.01124);  // allow growth, by default (=1/89)
    //     // if ((beta < 0.8) && (gamma < 0.8))  // arrest growth
    //     // if ((beta < 0.8) && (gamma < 0.8))  // arrest growth
    //     // if ((beta > 0.6) && (gamma > 0.6))
    //     if ((beta < 0.9) || (gamma < 0.9))
    //     {
    //         // set_single_behavior( pCell , "cycle entry" , 0.0);  
    //         // set_single_behavior( pCell , "cycle entry" , 0.01124);  // 1/89
    //         set_single_behavior( pCell , "cycle entry" , 0.0);  // arrest the cell cycle, by default
    //     }
    //     return;
    // }
    return; 
}

// void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
// { 
//     std::cout << __FUNCTION__ << ": " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << std::endl;
//     return; 
// } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

// void double_volume_update_function( Cell* pCell, Phenotype& phenotype , double dt )
// {
//     // Warning! Do not use get_total_volume!
//     // Use (some_cell).phenotype.volume.total instead!
//     // pCell->set_total_volume( pCell->phenotype.volume.total + pCell->phenotype.volume.total * 0.1);
//     // if ( get_single_signal( pCell, "pressure" ) < 3.0 )
//     // {
//     //     pCell->set_total_volume( pCell->phenotype.volume.total + pCell->phenotype.volume.total * pCell->custom_data["growth_rate"]);
//     // }
//     pCell->set_total_volume( pCell->phenotype.volume.total + pCell->phenotype.volume.total * pCell->custom_data["growth_rate"]);

//     // if (pCell->phenotype.volume.total > 1047)    //rwh: hard-coded; fix!
    // double NormalRandom( double mean, double standard_deviation )
    // if (pCell->phenotype.volume.total > NormalRandom(2.0, 0.4) * 523.6)    //rwh: hard-coded; fix!
//     // if (pCell->phenotype.volume.total > NormalRandom(2.0, 1.0) * 523.6)    //rwh: hard-coded; fix!
//     {
//         // std::cout << "------- " << __FUNCTION__ << ":  ID= " << pCell->ID <<":  volume.total= " << pCell->phenotype.volume.total << std::endl;
//         pCell->flag_for_division();
//     }
// }

// do every phenotype dt
void custom_volume_function( Cell* pCell, Phenotype& phenotype, double dt )
{ 
    // std::cout << " -------- "<<__FUNCTION__ << ":  t= " << PhysiCell_globals.current_time << ": cell ID= " << pCell->ID << std::endl;

    // pCell->set_total_volume( pCell->phenotype.volume.total + pCell->phenotype.volume.total * pCell->custom_data["growth_rate"]);
    pCell->set_total_volume( pCell->phenotype.volume.total + phenotype_dt * pCell->custom_data["growth_rate"]);

    // if (pCell->phenotype.volume.total > NormalRandom(2.0, 0.4) * 523.6)    //rwh: hard-coded; fix!
    double vol_norm_rand = pCell->custom_data["norm_rand"] * 523.6;
    if (pCell->phenotype.volume.total >= vol_norm_rand)  // hard-code
    {
        // std::cout << "-- flag for division because vol(total)= " <<  pCell->phenotype.volume.total<< " and vol_norm_rand= " <<  vol_norm_rand << std::endl;
        // pCell->flag_for_division();
        #pragma omp critical
        {
            // std::cout << "-- pCell->divide() because vol(total)= " <<  pCell->phenotype.volume.total<< " and vol_norm_rand= " <<  vol_norm_rand << std::endl;
            pCell->divide();
        }
    }
}

// Assign N(2,0.4^2) to each daughter cell. Exit the sim at 10K cells.
void custom_division_function( Cell* pCell1, Cell* pCell2 )
{ 
    // static int monolayer_max_cells = 9999;
    static int monolayer_max_cells = 10000;
    static double draw_mean = 2.0;
    static double draw_stddev = 0.4;

    // std::cout << " ---- "<<__FUNCTION__ << ":  t= " << PhysiCell_globals.current_time << std::endl;

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


    draw = NormalRandom(draw_mean, draw_stddev);  //     // where: NormalRandom(mean, std_dev)
    while (draw < 0.0)
    {
        draw = NormalRandom(draw_mean, draw_stddev); 
    }
    pCell2->custom_data["norm_rand"] = draw;

    return; 
}