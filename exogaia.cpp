#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <bitset>
#include <string>
#include <ios>

using namespace std;


//********* MICROBES ****************

class microbe {  // class for our microbe species
    
    public:
    int genome;
    int population;
    double nutrient; 
    int biomass;
    int waste;
    
};

//********** SORTING ****************

struct myclass {  // rank species by population, largest first
    bool operator() ( microbe i, microbe j) { return (i.population > j.population);}
} sorting_rule;


//********* CHOOSE INDIVIDUAL *********

// code to chose a random individual from the system
int chooseAgent(vector< microbe > &species, int p) {
    
    double r = 0;
    
    double p_num = drand48();
    
    for (int j = 0; j < species.size(); j++){
        r+= species[j].population;
        if (p*p_num <= r){
            return j;
        }
        
    }
    
    return -1; // happens if life goes extinct
    
}

//********* GREATER COMMON DENOMENATOR ***************

int GCD(int a, int b)
{
    while( 1 )
    {
        a = a % b;
		if( a == 0 )
			return b;
		b = b % a;

        if( b == 0 )
			return a;
    }
}

//******** NUTRIENT GENOME INTERACTIONS **************

vector < vector <int> > nutrient_genome_interactions(int genome_length, int num_nutrients, default_random_engine &generator){
    // this maps which geomones will code for which metabolisms
  
    vector <int> cons_ex_vector (2*num_nutrients, 0);  // consumption excretion vector
    vector <int> temp(num_nutrients, 0);
    
    vector < vector <int> > all_metabolisms;


    for (int i = 0; i < pow(2,genome_length); i++){

      int cons = 0;
      int exc = 0;

      for (int j = 0; j < 2*num_nutrients; j++) { cons_ex_vector[j] = 0; }

      cons = floor(drand48()*num_nutrients);
      exc  = floor(drand48()*num_nutrients);

      while (cons == exc) { exc = floor( drand48()*num_nutrients ); } // canot eat and excrete the same chemical species / nutrient

      cons_ex_vector[cons]  = 1;
      cons_ex_vector[exc+num_nutrients] = 1;
      
      all_metabolisms.push_back(cons_ex_vector);

    }
    
    return all_metabolisms;
    
}



//************ CREATE GEOLOGICAL LINKS  **************

vector < vector <double> >  environment_setup_links(double num_nutrients, double link_probability) {

  vector< vector <double> > geological_links( num_nutrients , vector<double>(num_nutrients, 0));

  for (int k = 0; k < num_nutrients; k++){

    for (int l = k+1; l < num_nutrients; l++){ // prevents repeats ie [1][2] = [2][1]

      if (drand48() <= link_probability && k != l) { // are the two linked at all
	double wt = drand48();   // 50 50 chance of nutrients flowing one way or the other
	double link_strength = drand48();

	if (wt < 0.5) {
	  geological_links[k][l] = link_strength; // flow from k to l
	  geological_links[l][k] = 0;
	} else {
	  geological_links[l][k] = link_strength; // flow from l to k
	  geological_links[k][l] = 0; // doesn't flow both ways
	} 
	
      }
      
    }
  }

  for (int k = 0; k < num_nutrients; k++){
    double total_out = 0.0;
    for (int l = 0; l < num_nutrients; l++){
      if (k != l) { total_out += geological_links[k][l]; }
      
    }
    geological_links[k][k] = total_out;
    if (total_out > 1) {
      for (int l = 0; l < num_nutrients; l++) {

	geological_links[k][l] = geological_links[k][l] / double(total_out);
	
      }
    }
  }

  return geological_links;

}

//****************************** ABIOTIC TRICKLE *********************

double update_abiotic(vector <double> &environment, int num_nutrients, vector <double> node_abiotic, double abiotic_env, double abiotic_T){

      double reflected = 0.0;
      double insulated = 0.0;
      
      for (int j = 0; j < num_nutrients; j++){
	
	if (node_abiotic[j] < 0) { reflected += -1.0*(node_abiotic[j]*100) * tanh(environment[j]/75000.0); } 
	else if (node_abiotic[j] > 0) { insulated += (node_abiotic[j]*100) * tanh(environment[j]/75000.0); } 

      }

      if (reflected > 100) { reflected = 100; }
      if (insulated > 100) { insulated = 100; }
      
      double current_env = abiotic_T;

      double incoming = abiotic_env * (100.0 - reflected)/100.0;
      double residual = current_env * insulated / 100.0;

      return incoming + residual - abiotic_T;
      
}


//****************** UPDATE ENVIRONMENT *****************
vector <double> update_environment(vector <double> &environment, int num_nutrients, double percentage_outflow, vector < vector <double> > &geological_links, vector <double> &influx_nodes){

  vector <double> nutrient_trickle (num_nutrients,0);


  for (int j = 0; j < num_nutrients; j++) { // outflow inflow
    
    nutrient_trickle[j] = environment[j]*(1.0-percentage_outflow) + influx_nodes[j];

  }
  
   for (int j = 0; j < num_nutrients; j++) { // geological processes
    double cur = nutrient_trickle[j];
 
    for (int l = 0; l < num_nutrients; l++){
      
      nutrient_trickle[l] += cur*geological_links[j][l];
      nutrient_trickle[j] -= cur*geological_links[j][l];
      
    }
	      	      
  }

  for (int j = 0; j < num_nutrients; j++){
    nutrient_trickle[j] -= environment[j]; // difference with the current environment
  }

  return nutrient_trickle;
  
}


// ************* SET UP INFLUX NODES / CHEMICAL SPECIES **************************
vector < double > create_nodes(int num_nutrients, int num_source, int max_nutrient_inflow){
  
  vector < double > influx_nodes (num_nutrients, 0);
 
  
  for (int j = 0; j < num_source; j++){  // create the source nodes randomly

    int new_nut = 0;
    int loc_n = floor(drand48()*num_nutrients);
    while (new_nut == 0) {
      if (influx_nodes[loc_n] == 0) { influx_nodes[loc_n] = max_nutrient_inflow; new_nut = 1; }
      else { loc_n = floor(drand48()*num_nutrients); }
	
    }
      
  }
  return influx_nodes;

}


//***** SET UP INSULATING / REFLECTING PROPERTIES OF EACH CHEMICAL SPECIES ******
vector < double > ab_nodes(int num_nutrients, double abiotic_prob){

  vector < double > node_abiotic (num_nutrients, 0);
 
  for (int j = 0; j < num_nutrients; j++){

    if (drand48() < abiotic_prob) { // chance chemical species will have non-zero effect on temp
      node_abiotic[j] = (2.0*drand48()-1.0); // +ve means insulating, -ve means reflecting 
    } 

  }

  return node_abiotic;

}


//************ MAIN CODE ************************

int main(int argc, char **argv) {
    
    
    //***** MICROBE PARAMETERS*****
    int i;                                        // this is a marker for chosing and individual
    const int initial_population = 100;           // initial population
    const int genome_length = 8;                  // length of genes in genome, each gene can be either '0' or '1'
    const int reproduction_thresh = 120;          // biomass threashold to reproduce
    const int starve_thresh = 50;                 // if somethings biomass drops below this, it dies
    const int max_consumption = 10;               // maximum number of nutrients a microbe can eat at once
    const double nutrient_conversion_eff = 0.6;   // efficiency of microbe metabolism
    const int maintainence_cost = 1;              // how much biomass it costs per timestep to live
    const double p_mut = 0.01;                    // probability per gene of mutation in reproduction event
    const double p_kill = 0.002;                  // probability of dead due to causes other than starvation
    const double prefered_abiotic = 1000.0;       // abstract temperature
    const bool reseed_on = false;                 // reseed with life after extinction?
    
    
    //***** ENVIRONMENT PARAMETERS *******
    const int num_nutrients = 8;                  // number of different types of nutrients/chemicals in system
    const int num_source = 2;                     // number of nutrients with influx. Must be <= num_nutrients!
    const double percentage_outflow = 0.0001;     // to calculate outflow
    const int max_nutrient_inflow = 75;           // per source node
    const double abiotic_start = 500.0;           // incoming heat from 'sun'. If start and end are different, world wil gradually warm/cool
    const double abiotic_end = 500.0;
    double abiotic_env = abiotic_start;           // starting environmental temperature
    double abiotic_trickle = 0;                   // used to update environment between timestep iterations
    double abiotic_prob = 1.0;                    // probability of heating / cooling for each node (chemical species)

    //***** RANDOM NUMBER GENERATORS ****
    int t1 = atoi(argv[1]);                       // number to initialise the chemical set (heating / cooling properties)
    srand48 (t1);
    mt19937 rng(t1);
    default_random_engine generator;
    generator.seed(t1);                           // provide seed for randomness in different runs
    
    vector < double > node_abiotic = ab_nodes(num_nutrients, abiotic_prob); // abiotic affect of the nodes
    vector < double > influx_nodes = create_nodes(num_nutrients, num_source, max_nutrient_inflow);
    double abiotic_T = abiotic_start;                 // set initial temperature to the start temperature
    vector<double> environment(num_nutrients, 0);     // initally no chemicals in 'atmosphere'
    vector<double> nutrient_trickle(num_nutrients,0); // for updating environment in between iterations

    int t2 = atoi(argv[2]);                           // number to initialise geochemistry (links)
    srand48 (t2);
    generator.seed(t2); 
    
    const double link_probability = 0.4;              // how likely for two nodes to be connected
    vector< vector<double> > geological_links = environment_setup_links(num_nutrients, link_probability);    


    //***** RANDOM NUMBER GENERATORS ****
    int t = atoi(argv[3]);  // number to initialise microbe metabolisms
    srand48 (t);
    generator.seed(t);      // provide seed for randomness in different runs

    
    
    // METABOLISM SET UP
    vector < vector< int > >  n_g_interacts = nutrient_genome_interactions(genome_length, num_nutrients, generator);
    
    //***** SPECIES VARIABLES ********
    vector < microbe > species; 
    int total_population = initial_population;
    double average_biomass;
    int i_biomass;
    double species_nutrient_avg;
    double species_biomass_avg;
    int nutrient_available;
    const double abiotic_scaling = 0.015; // microbe sensitivity to temperature
    double satisfaction;                  // a measure of how 'fit' the microbes are in their current environment
    double factor_i;
    microbe temp_mutant;
    int did_we_mutate = 0;


    int genome_new = floor(drand48()*pow(2,genome_length)); // randomly generate microbes
    microbe new_microbe;
    new_microbe.population = initial_population;
    new_microbe.genome = genome_new;
    new_microbe.nutrient = 0;
    new_microbe.biomass = 80*initial_population;
    new_microbe.waste = 0;
    species.push_back(new_microbe);

    
    // VARIABLES FOR KEEPING TRACK OF TIME
    int timestep_length = total_population; // number of timestep iterations determined by total population at start of timestep
    int number_gens = 0;                    // number of timesteps that have passed
    int timestep_counter = 0;               // for counting iterations within a timestep
    int max_timesteps = 50*pow(10,4);       // max length of experiment
    int init_period = 5*pow(10,4);          // after this time has passed if habitable conditions haven't been reached, seed anyway
    int init_counter = 0;                   // for tracking time before life is seeded
    int non_ideal = 1;                      // switches to 0 when environment is ideal and life can be seeded
    int death_iteration = 0; 		    // if a microbe dies at the start of iteration, skip to next iteration (a dead microbe can't eat etc)
    
    // DATA FILES
    int file_num = atoi(argv[4]);                                                     // data file number
    ofstream macro_data ("exogaia_macro_data_"+to_string(file_num)+".txt");           // macro properties - total pop, temp, etc
    ofstream pop_data ("exogaia_pop_data_"+to_string(file_num)+".txt");               // population of each species alive at timestep
    ofstream nutrient_data ("exogaia_nutrient_data_"+to_string(file_num)+".txt");     // chemical species levels over time
    ofstream genome_data ("exogaia_genome_data_"+to_string(file_num)+".txt");         // which genomes exist over time
    ofstream nutrient_genome ("exogaia_nutrient_genome_"+to_string(file_num)+".txt"); // chemicals being consumed at each timestep
    ofstream waste_genome ("exogaia_waste_data_"+to_string(file_num)+".txt");         // chemicals being excreted over time
    ofstream geological_net ("exogaia_geological_network_dat.txt");                   // the geochemical network

    for (int j = 0; j < num_nutrients; j++){
      geological_net << influx_nodes[j] << " ";

    }
    geological_net << endl;
    for (int j = 0; j < num_nutrients; j++){
      geological_net << node_abiotic[j] << " ";

    }
    geological_net << endl;
    
    for (int k = 0; k < num_nutrients; k++){
      for (int l = 0; l < num_nutrients; l++){
	geological_net << geological_links[k][l] << " ";

      }
      geological_net << endl;
    }

    geological_net.close();

    while (init_counter < init_period && non_ideal == 1) {   // INITIALISE OUR ENVIRONMENT

      /* ********************************************************************************
                                    NUTRIENT FLOW
      ********************************************************************************/

            
      // only update the flow once every time step       
      // nutrient outflow

      abiotic_trickle = update_abiotic(environment, num_nutrients, node_abiotic, abiotic_env, abiotic_T);
      abiotic_T += abiotic_trickle; 
      nutrient_trickle = update_environment(environment, num_nutrients, percentage_outflow, geological_links, influx_nodes);
      for (int j = 0; j < num_nutrients; j++){

	environment[j] += nutrient_trickle[j];
	nutrient_trickle[j] = 0;

      }
	abiotic_trickle = 0;
	
       
	if (abiotic_T >= 1000 && abiotic_T <= 1050) {      // seeding window
		for (int k = 0; k < num_nutrients; k++){
			if (environment[k] > 1000) {
				 non_ideal = 0;            // environment is suitable for seeding with life
			}
		}
	} // seed once conditions are habitable for life 
      init_counter++;

    	  double reflected = 0.0;
    	  double insulated = 0.0;
      
      	for (int j = 0; j < num_nutrients; j++){
	
		if (node_abiotic[j] < 0) { reflected += -1.0*(node_abiotic[j]*100) * tanh(environment[j]/75000.0); }
		else if (node_abiotic[j] > 0) { insulated += (node_abiotic[j]*100) * tanh(environment[j]/75000.0); }
	
      	}

	if (reflected > 100) { reflected = 100; }
      	if (insulated > 100) { insulated = 100; }

	macro_data << init_counter << " 0 0 0 " << abiotic_T << " " << insulated << " " << reflected << endl;

    }

    //************************** SEED **************************************


	int suit_metab = 0; // find a suitable metabolism - food suitable for the species' metabolism must be available

	while (suit_metab == 0) {
		species[0].genome = floor(drand48()*pow(2,genome_length));
		for (int q = 0; q < num_nutrients; q++){
			if (n_g_interacts[species[0].genome][q] > 0 && environment[q] > 1000) {  // food source available?
				suit_metab = 1;
			}
		}
	}   

    
    while (number_gens < max_timesteps) {
          
        /* ********************************************************************************
                                RECORD DATA
         ********************************************************************************/
        if (timestep_counter >= timestep_length){
	  
            timestep_counter = 0;
            timestep_length = total_population;
            number_gens++;
 
            stable_sort (species.begin(), species.end(), sorting_rule);

    	  double reflected = 0.0;
    	  double insulated = 0.0;
      
      	for (int j = 0; j < num_nutrients; j++){
	
		if (node_abiotic[j] < 0) { reflected += -1.0*(node_abiotic[j]*100) * tanh(environment[j]/75000.0); }
		else if (node_abiotic[j] > 0) { insulated += (node_abiotic[j]*100) * tanh(environment[j]/75000.0); }
	
      	}

	if (reflected > 100) { reflected = 100; }
      	if (insulated > 100) { insulated = 100; }

	  factor_i = abiotic_scaling*sqrt(pow(abiotic_T - prefered_abiotic, 2.0));
	  satisfaction = exp (-1.0*pow(factor_i,2.0));

	  int total_count_eat = floor(max_consumption * satisfaction);

            // Record data here!!!!!
            macro_data << number_gens+init_counter << " " << total_population << " " << species[0].population << " " << species.size() << " " << abiotic_T << " " << insulated << " " << reflected <<  endl;

	    genome_data << number_gens;
	    pop_data << number_gens;
	    
	    nutrient_genome << number_gens;
	    waste_genome << number_gens;

            for (int j = 0; j < species.size(); j++){
		  pop_data << " " << species[j].population;
		  genome_data << " " << species[j].genome;
		  
		  int eat_in = 0;
		  int waste_out = 0;
		  
		  for (int l = 0; l < num_nutrients; l++) {
		    if (n_g_interacts[species[j].genome][l] > 0) { eat_in = l+1; }
		    if (n_g_interacts[species[j].genome][l+num_nutrients] > 0) { waste_out = l+1; }
		  }

		  nutrient_genome << " " << eat_in;
		  waste_genome << " " << waste_out; 
            }

            pop_data << endl;
	    genome_data << endl;
	    nutrient_genome << endl;
	    waste_genome << endl;

	    nutrient_data << number_gens;
	    
	    for (int j = 0; j < num_nutrients; j++) {

	      nutrient_data << " " << environment[j];
	      
	    }
	    
	    nutrient_data << endl; 

	    /************************************************************************************
	                       CALCULATE NUTRIENT TRICKLE
	    **************************************************************************************/

	    abiotic_env += (abiotic_end - abiotic_start) / max_timesteps;
	    abiotic_trickle = update_abiotic(environment, num_nutrients, node_abiotic, abiotic_env, abiotic_T);

	    if (timestep_length > 0) { abiotic_trickle = abiotic_trickle / (1.0*timestep_length); }

	    nutrient_trickle = update_environment(environment, num_nutrients, percentage_outflow, geological_links, influx_nodes);
	    if (timestep_length > 0) {
	      
	      for (int j = 0; j < num_nutrients; j++){
		nutrient_trickle[j] = nutrient_trickle[j]/(1.0*timestep_length);
	      }
	    }

	}
        
            /* ********************************************************************************
                                    NUTRIENT FLOW
             ********************************************************************************/

            // have a trickle every iteration adding up to the alloted count per timestep     
            // nutrient outflow
            for (int j = 0; j < num_nutrients; j++ ){

	      environment[j] += nutrient_trickle[j];
	      if (environment[j] < 0) { environment[j] = 0; }
	      
            }

	    
	    abiotic_T += abiotic_trickle; 


        /*********************************************************************************
                 RESEED IF PLANET IS EXTINCT (Only happens if reseed_on == True)
        *********************************************************************************/

	if (abiotic_T >= 1000 && abiotic_T < 1050 && total_population == 0 && reseed_on){  

	     int genome_news = floor(drand48()*pow(2,genome_length)); // randomly generate microbes
	     microbe try_microbe;
	     try_microbe.population = initial_population;
	     try_microbe.genome = genome_news;
	     try_microbe.nutrient = 0;
	     try_microbe.biomass = 80*initial_population;
	     try_microbe.waste = 0;
	     species.push_back(try_microbe);
	     total_population = initial_population;

              	int suit_metab = 0; // find a suitable metabolism for current environment
	
		while (suit_metab == 0) {
			species[0].genome = floor(drand48()*pow(2,genome_length));
			for (int q = 0; q < num_nutrients; q++){
				if (n_g_interacts[species[0].genome][q] > 0 && environment[q] > 1000) {
					suit_metab = 1;
				}
			}
		}   


        }
        
        /* ********************************************************************************
                                    KILL
         ********************************************************************************/
        
        // NEED TO REMOVE BIOMASS WHEN AN INDIVIDUAL DIES

        i = chooseAgent(species, total_population);

	death_iteration = 0; // Reset at start of each iteration. Becomes 1 if chosen microbe at start of iteration dies
			     // If a microbe dies it cannot eat / reproduce etc therefore there is one less possible eating etc event 
			     // within the current timestep
        if (i > -1) {

	  // death event starvation
	  average_biomass = species[i].biomass/(1.0*species[i].population);
	  normal_distribution<double> biomass_dist( average_biomass, average_biomass*0.01 ); // distribution of biomass in population
	  i_biomass = floor(biomass_dist(generator));

	  species_nutrient_avg = 1.0*species[i].nutrient/species[i].population;
	  normal_distribution<double> nutrient_species_dist( species_nutrient_avg, species_nutrient_avg*0.1);
	  nutrient_available = floor(nutrient_species_dist(generator)); // we'll just round down as can't use half a nutrient

	  if (i_biomass <= starve_thresh) {
            // dies if biomass lower than starvation threshold 

            species[i].population--;
            species[i].biomass -= i_biomass;           // remove biomass of dead microbe
	    species[i].nutrient -= nutrient_available; // remove the undigested food of the dead microbe
            total_population--;
            if (species[i].biomass < 1) { species[i].biomass = 0; species[i].population = 0; } // if there is no biomass, extinct
            if (species[i].population == 0){
	      species.erase(species.begin() + i);      // remove species from list if extinct

            }
	    death_iteration = 1;
	  }

	  else if (drand48() <= p_kill && species[i].population > 0) {            
            species[i].population--;
            species[i].biomass -= i_biomass;           // remove biomass of dead microbe
	    species[i].nutrient -= nutrient_available; // remove the undigested food of the dead microbe
            total_population--;
            if (species[i].population == 0){
	      species.erase(species.begin() + i);      // remove from species list if extinct
            }
	    death_iteration = 1;
	  }
	}


        /* ********************************************************************************
                                MAINTENANCE COST
         ********************************************************************************/
        i = chooseAgent(species, total_population);
	
        if (i > -1 && death_iteration == 0) {
	  species[i].biomass--;
	}

        /* ********************************************************************************
                                METABOLISM
         ********************************************************************************/

        // metabolism event
        i = chooseAgent(species, total_population);
	
	if( i > -1 && death_iteration == 0) {
	  
	  factor_i = abiotic_scaling*sqrt(pow(abiotic_T - prefered_abiotic, 2.0));
	  satisfaction = exp (-1.0*pow(factor_i,2.0));

	  int minimum_count_eat = 0; // the minumum total number of nutrients microbe can intake
	  int max_count_eat = 0;
	  int nut_num;
	  for (int k = 0; k < num_nutrients; k++) {

	    if (n_g_interacts[species[i].genome][k] > 0) { nut_num = k;  } // which nutrient / chemical species does this microbe eat?
	  }

	  double total_count_eat = max_consumption * satisfaction;
	 if (environment[nut_num] < total_count_eat) { total_count_eat = environment[nut_num]; }
	 environment[nut_num] -= total_count_eat;
	 species[i].nutrient += total_count_eat;

	 if (environment[nut_num] < 0) { cout << "NUTRIENT EATING PROBLEM" << endl; } // bug check - has never happened

	}
	
        /* ********************************************************************************
                                BIOMASS CREATION
         ********************************************************************************/

        i = chooseAgent(species, total_population);

	if (i > -1 && death_iteration == 0){

	  species_nutrient_avg = 1.0*species[i].nutrient/species[i].population;

	  normal_distribution<double> nutrient_species_dist( species_nutrient_avg, species_nutrient_avg*0.1);
	  nutrient_available = floor(nutrient_species_dist(generator)); // we'll just round down as can't use half a nutrient

	  while ( nutrient_available >= 5) { 
           
            species[i].nutrient -= 5;
            nutrient_available -= 5;
            species[i].biomass += int(5.0*nutrient_conversion_eff);
	    species[i].waste += int(5*(1.0 - nutrient_conversion_eff));
	    
	  }
	}


	/*********************************************************************************
                                  WASTE 
	 **********************************************************************************/

        i = chooseAgent(species, total_population);
	if (i > -1 && death_iteration == 0){

	  double species_waste_avg = 1.0*species[i].waste/species[i].population;
	  normal_distribution<double> waste_species_dist( species_waste_avg, species_waste_avg*0.1);
	  int waste_available = floor(waste_species_dist(generator));

	  if (waste_available > species[i].waste) { waste_available = species[i].waste; }

	  for (int k = 0; k < num_nutrients; k++) {

		if (n_g_interacts[species[i].genome][k+num_nutrients] > 0) {  // microbe excretes this chemical species as waste

			species[i].waste -= waste_available;
			environment[k] += waste_available;
		}
	  }

	}


        /* ********************************************************************************
                                REPRODUCTION
         ********************************************************************************/
        
        // reproduction event

        i = chooseAgent(species, total_population);

        if (i > -1 && death_iteration == 0){
        
	  species_biomass_avg = (1.0*species[i].biomass)/species[i].population; // average biomass per indiviual
        
        
	  normal_distribution<double> biomass_species_dist( species_biomass_avg, species_biomass_avg*0.01 );
	  i_biomass = floor(biomass_species_dist(generator));
        
	  if (i_biomass >= reproduction_thresh) {

            
            // DO WE MUTATE?
            
            bitset<genome_length> mutant_genome(species[i].genome);
            did_we_mutate = 0;
            
            for (int j = 0; j < genome_length; j++){
	      if (drand48() <= p_mut){
		did_we_mutate = 1;
		if (mutant_genome[j] == 1) {
		  mutant_genome[j] = 0;
		} else {
		  mutant_genome[j] = 1;
		}
	      }
            }
            
            if (did_we_mutate == 1) {
	      int mutant_number = int(mutant_genome.to_ulong());
	      int species_exists = 0;
	      for (int q = 0; q < species.size(); q++){ // check to see if species exists
		if (species[q].genome == mutant_number){
		  species[q].population++;
		  species[q].biomass += int(i_biomass / 2.0); // half biomass goes to new mutant
		  species[i].biomass -= int(i_biomass / 2.0);
		  species_exists = 1;
		  break;
		}
	      }
                
	      if (species_exists == 0){ // add species if it doesn't exist
		temp_mutant.genome = mutant_number;
		temp_mutant.nutrient = 0; //  no nutrient count to begin with
		temp_mutant.biomass = int(i_biomass / 2.0);  // half biomass goes to new mutant
		species[i].biomass -= int(i_biomass / 2.0);
		temp_mutant.population = 1; // initial population of 1
		temp_mutant.waste = 0;
		species.push_back(temp_mutant);
	      }
                
                
	      total_population++;
                
            } else { // no mutation takes place, we add one to the population
	      species[i].population++;
	      total_population++;
            }

        
	  }
        }

        timestep_counter++; // increment our timestep counter

        /*********** END OF WHILE LOOP ***************/
        
    }
    
    macro_data.close();
    pop_data.close(); 
    nutrient_data.close();
    genome_data.close();
    nutrient_genome.close();
    waste_genome.close();
    
    return 0;
}
