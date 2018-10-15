

#include <cstring>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/feature/FeatureNNInteractionSc.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>

#include "catchorg/clara/clara.hpp"

#include "Analyzer_EigenvaluesRouseMatrix.h"

int main(int argc, char* argv[])
{
  try{
	std::string infile  = "input.bfm";
	
    bool showHelp = false;
    
    auto parser
    = clara::Opt( infile, "input (=input.bfm)" )
        ["-i"]["--infile"]
        ("BFM-file to load and analyzer.")
        .required()
    | clara::Help( showHelp );
        
    auto result = parser.parse( clara::Args( argc, argv ) );
    if( !result ) {
    std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
    exit(1);
    }
    else if(showHelp == true)
    {
        std::cout << "Analyzer for calculating the eigenvalues of Rouse Matrix" << std::endl
                  << "providing the (ideal) <Rg²> of the topology" << std::endl
                  << "Features used: FeatureMoleculesIO, FeatureAttributes" << std::endl
		          << "Updaters used: ReadBFMFile, SimpleSimulator" << std::endl
		          << "Analyzers used: WriteBfmFile" << std::endl;
        
        parser.writeToStream(std::cout);
        exit(0);
    }
    else
    {
        std::cout << "infile:        " << infile << std::endl;
    }
       
	//seed the globally available random number generators
	RandomNumberGenerators rng;
	rng.seedAll();
	
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	typedef LOKI_TYPELIST_2(FeatureMoleculesIO, FeatureAttributes) Features;
	
	typedef ConfigureSystem<VectorInt3,Features, 6> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	TaskManager taskmanager;
	taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE));
	
	taskmanager.addAnalyzer(new Analyzer_EigenvaluesRouseMatrix<Ing>(myIngredients, "."));

	taskmanager.initialize();
	taskmanager.run(1);
	taskmanager.cleanup();
	
	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

