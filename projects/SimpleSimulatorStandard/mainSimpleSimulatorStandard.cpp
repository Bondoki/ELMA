

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

int main(int argc, char* argv[])
{
  try{
	std::string infile  = "input.bfm";
	std::string outfile = "outfile.bfm";
	uint32_t max_mcs=100;
	uint32_t save_interval=100;
	
    bool showHelp = false;
    
    auto parser
    = clara::Opt( infile, "input (=input.bfm)" )
        ["-i"]["--infile"]
        ("BFM-file to load.")
        .required()
    | clara::Opt( outfile, "output (=outfile.bfm)" )
        ["-o"]["--outfile"]
        ("BFM-file to save.")
        .required()
    | clara::Opt( [&max_mcs](int const m)
        {
         if (m <= 0)
         {
            return clara::ParserResult::runtimeError("Simulation time must be greater than 0");
         }
         else
         {
            max_mcs = m;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }
        }, "max MCS(=100)" )
        ["-m"]["--max-mcs"]
        ("(required) specifies the total Monte-Carlo steps to simulate.")
        .required()
    | clara::Opt( [&save_interval](int const s)
        {
         if (s < 0)
         {
            return clara::ParserResult::runtimeError("Save intervall must be greater than 0");
         }
         else
         {
            save_interval = s;
            return clara::ParserResult::ok(clara::ParseResultType::Matched);
         }
      }, "save MCS(=100)")
        ["-s"]["--save-mcs"]
        ("(required) Save after every <integer> Monte-Carlo steps to the output file." )
        .required()
    | clara::Help( showHelp );
        
    auto result = parser.parse( clara::Args( argc, argv ) );
    if( !result ) {
    std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
    exit(1);
    }
    else if(showHelp == true)
    {
        std::cout << "Simulator for the ScBFM with Ex.Vol and BondCheck" << std::endl
                  << "maximum number of connections per monomer is 6" << std::endl
                  << "Features used: FeatureBondset, FeatureAttributes, FeatureExcludedVolumeSc<FeatureLattice<bool> >" << std::endl
		          << "Updaters used: ReadFullBFMFile, SimpleSimulator" << std::endl
		          << "Analyzers used: WriteBfmFile" << std::endl;
        
        parser.writeToStream(std::cout);
        exit(0);
    }
    else
    {
        std::cout << "infile:        " << infile << std::endl
                  << "outfile:       " << outfile << std::endl
                  << "max_mcs:       " << max_mcs << std::endl
                  << "save_interval: " << save_interval << std::endl;
    }
       
	//seed the globally available random number generators
	RandomNumberGenerators rng;
	rng.seedAll();
	
	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes< >,FeatureExcludedVolumeSc<>) Features;
	
	typedef ConfigureSystem<VectorInt3,Features, 6> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;

	TaskManager taskmanager;
	taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
	//here you can choose to use MoveLocalBcc instead. Careful though: no real tests made yet
	//(other than for latticeOccupation, valid bonds, frozen monomers...)
	taskmanager.addUpdater(new UpdaterSimpleSimulator<Ing,MoveLocalSc>(myIngredients,save_interval));

	taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));
	
	taskmanager.initialize();
	taskmanager.run(max_mcs/save_interval);
	taskmanager.cleanup();
	
	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

