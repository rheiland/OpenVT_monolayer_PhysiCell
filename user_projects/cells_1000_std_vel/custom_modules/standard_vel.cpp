in core/P*_standard_models.cpp:

void standard_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	pCell->state.neighbors.clear(); // new 1.8.0
	
	//First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		pCell->add_potentials(*neighbor);
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
			pCell->add_potentials(*neighbor);
		}
	}

	pCell->update_motility_vector(dt); 
	pCell->velocity += phenotype.motility.motility_vector; 
	
	return; 
}

---------------------------------------------
void Cell::add_potentials(Cell* other_agent)
{
	// if( this->ID == other_agent->ID )
	if( this == other_agent )
	{ return; }
/*
	// new April 2022: don't interact with cells with 0 volume 
	// does not seem to really help 
	if( other_agent->phenotype.volume.total < 1e-15 )
	{ std::cout << "zero size cell in mechanics!" << std::endl; return; }
*/
	// 12 uniform neighbors at a close packing distance, after dividing out all constants
	static double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2 
	// 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2

	double distance = 0; 
	for( int i = 0 ; i < 3 ; i++ ) 
	{ 
		displacement[i] = position[i] - (*other_agent).position[i]; 
		distance += displacement[i] * displacement[i]; 
	}
	// Make sure that the distance is not zero
	
	distance = std::max(sqrt(distance), 0.00001); 
	
	//Repulsive
	double R = phenotype.geometry.radius+ (*other_agent).phenotype.geometry.radius; 
	
	// double RN = phenotype.geometry.nuclear_radius + (*other_agent).phenotype.geometry.nuclear_radius;	
	double temp_r, c;
	if( distance > R ) 
	{
		temp_r=0;
	}
	// else if( distance < RN ) 
	// {
		// double M = 1.0; 
		// c = 1.0 - RN/R; 
		// c *= c; 
		// c -= M; 
		// temp_r = ( c*distance/RN  + M  ); 
	// }
	else
	{
		// temp_r = 1 - distance/R;
		temp_r = -distance; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2 
		
		// add the relative pressure contribution 
		state.simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017 
	}
	
	// August 2017 - back to the original if both have same coefficient 

	double effective_repulsion = sqrt( phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength ); 
	temp_r *= effective_repulsion; 
	
	// temp_r *= phenotype.mechanics.cell_cell_repulsion_strength; // original 
	//////////////////////////////////////////////////////////////////
	
	// Adhesive
	//double max_interactive_distance = parameters.max_interaction_distance_factor * phenotype.geometry.radius + 
	//	(*other_agent).parameters.max_interaction_distance_factor * (*other_agent).phenotype.geometry.radius;
		
	double max_interactive_distance = phenotype.mechanics.relative_maximum_adhesion_distance * phenotype.geometry.radius + 
		(*other_agent).phenotype.mechanics.relative_maximum_adhesion_distance * (*other_agent).phenotype.geometry.radius;
		
	if(distance < max_interactive_distance ) 
	{	
		// double temp_a = 1 - distance/max_interactive_distance; 
		double temp_a = -distance; // -d
		temp_a /= max_interactive_distance; // -d/S
		temp_a += 1.0; // 1 - d/S 
		temp_a *= temp_a; // (1-d/S)^2 
		// temp_a *= phenotype.mechanics.cell_cell_adhesion_strength; // original 
		
		// August 2017 - back to the original if both have same coefficient 
		// May 2022 - back to oriinal if both affinities are 1
		int ii = find_cell_definition_index( this->type ); 
		int jj = find_cell_definition_index( other_agent->type ); 

		double adhesion_ii = phenotype.mechanics.cell_cell_adhesion_strength * phenotype.mechanics.cell_adhesion_affinities[jj]; 
		double adhesion_jj = other_agent->phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_adhesion_affinities[ii]; 

		// double effective_adhesion = sqrt( phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_cell_adhesion_strength ); 
		double effective_adhesion = sqrt( adhesion_ii*adhesion_jj ); 
		temp_a *= effective_adhesion; 
		
		temp_r -= temp_a;

		state.neighbors.push_back(other_agent); // move here in 1.10.2 so non-adhesive cells also added. 
	}
	/////////////////////////////////////////////////////////////////
	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distance;
	// for( int i = 0 ; i < 3 ; i++ ) 
	// {
	//	velocity[i] += displacement[i] * temp_r; 
	// }
	axpy( &velocity , temp_r , displacement ); 
	
	// state.neighbors.push_back(other_agent); // new 1.8.0
	return;
}
