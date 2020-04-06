

#include <omp.h>

#include <cstring>
#include <sstream>      // std::stringstream, std::stringbuf


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

#include "FeatureWangLandauNextNeighbor.h"
#include "UpdaterAdaptiveWangLandauSamplingNextNeighbor.h"
#include "ReadInHGLnDOS.h"


int main(int argc, char* argv[])
{
	try{
		std::string infile  = "input.bfm";
		std::string outfile = "outfile.bfm";
		std::string HGLnDOSfile = "";

		uint64_t max_mcs=100;
		uint32_t save_interval=100;
		uint32_t bias_update_interval=100;
		double min_histogram = -100.0;
		double max_histogram = +100.0;
		uint32_t bins_histogram = 200;
		double modFactor = 1.01;

		double minWin = -100.0;
		double maxWin = +100.0;

		bool showHelp = false;

		bool filedump = false;

		double overlap = 0.75;

		double lengthIncrease = 0.25;

		bool readinBFMinWin = false;

		double flatness = 0.85;

		auto parser
		= clara::Opt( infile, "input (=input.bfm)" )
		["-i"]["--infile"]
			   ("BFM-file to load.")
			   .required()
		| clara::Opt( outfile, "output (=outfile.bfm)" )
		["-o"]["--outfile"]
			   ("BFM-file to save.")
			   .required()
		| clara::Opt(  min_histogram, "min histogram (=-100.0)" )
		["--min"]
			   ("minimal histogram boundary (=-100.0)")
			   .required()
		| clara::Opt(  max_histogram, "max histogram (=+100.0)" )
		["--max"]
			   	("maximal histogram boundary (=+100.0)")
			   	.required()
		| clara::Opt(  bins_histogram, "bins histogram (=+200.0)" )
		["--bins"]
				("bins histogram boundary (=+200.0)")
				.required()
		| clara::Opt(  modFactor, "modification factor (=1.01)" )
		["-f"]["--mod-factor"]
				("initial modification factor for update DOS (=1.01)")
				.required()
		| clara::Opt( [&max_mcs](uint64_t const m)
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
		["-r"]["--replica-exchange-time"]
			   ("(required) Time in MCS to try a replica exchange move." )
			   .required()
		| clara::Opt( [&bias_update_interval](int const b)
			   	{
			   					if (b < 0)
			   					{
			   						return clara::ParserResult::runtimeError("bias_update_interval must be greater than 0");
			   					}
			   					else
			   					{
			   						bias_update_interval = b;
			   						return clara::ParserResult::ok(clara::ParseResultType::Matched);
			   					}
			   	}, "bias_update_interval MCS(=100)")
			   		["-b"]["--bias-update-interval"]
			   			   ("(required) Update histogram interval after every <integer> Monte-Carlo steps to check the DOS." )
			   			   .required()
		| clara::Opt(  minWin, "min window (=-100.0)" )
			["--min-win"]
			 ("minimal window boundary (=-100.0)")
			 .required()
		| clara::Opt(  maxWin, "max window (=+100.0)" )
			["--max-win"]
			("maximal window boundary (=+100.0)")
			.required()
		| clara::Opt( HGLnDOSfile, "HGLnDOS" )
			["--HGLnDOS"]
			("Histogram file for (logarithmic) DOS to load.")
			.required()
		| clara::Opt( filedump, "filedump" )
			["--dump"]
			("boolean: dump file every save_interval.")
			.required()
		| clara::Opt( overlap, "overlap (=0.75)" )
			["--overlap"]
			 ("Fraction of overlap [0,1] between neighboring windows.")
			.required()
		| clara::Opt( lengthIncrease, "lengthIncreaseWindow (=0.25)" )
			["--length-increase"]
			("Fractional increase of successive window length (0,1] between neighboring windows.")
			.required()
		| clara::Opt(  readinBFMinWin, "read separate bfm file for every window (=false)" )
			["--read-in-BFM"]
			("every window is initialized with separate bfm file (=false)")
			.required()
		| clara::Opt(  flatness, "flatness criterion of histogram iteration (=0.85)" )
			["--flatness"]
			("flatness criterion of histogram iteration (=0.85)")
			.required()
		| clara::Help( showHelp );

		auto result = parser.parse( clara::Args( argc, argv ) );
		if( !result ) {
			std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
			exit(1);
		}
		else if(showHelp == true)
		{
			std::cout << "Simulator for the ScBFM with Ex.Vol and BondCheck and WangLandau in NN-shell" << std::endl
					<< "maximum number of connections per monomer is 6" << std::endl
					<< "Features used: FeatureBondset, FeatureAttributes, FeatureExcludedVolumeSc<FeatureLattice<uint8_t> >, FeatureWangLandauNextNeighbor" << std::endl
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
					<< "save_interval: " << save_interval << std::endl
					<< "bias_update_interval: " << bias_update_interval << std::endl
					<< "min_histogram: " << min_histogram << std::endl
					<< "max_histogram: " << max_histogram << std::endl
					<< "bins_histogram: " << bins_histogram << std::endl
					<< "modFactor: " << modFactor << std::endl
					<< "min_win: "	<< minWin << std::endl
					<< "max_win: " 	<< maxWin << std::endl
					<< "HGLnDOS:" << HGLnDOSfile << std::endl
					<< "filedump:" << filedump << std::endl
					;
		}
	
	int nthreads, tid;

	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	//typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureFixedMonomers,FeatureAttributes,FeatureExcludedVolumeSc<>) Features;
	//typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes, FeatureNNInteractionSc< FeatureLattice >) Features;
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes, FeatureWangLandauNextNeighbor< FeatureLattice >) Features;
	//typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureWangLandau, FeatureAttributes) Features;

	typedef ConfigureSystem<VectorInt3,Features, 6> Config;
	typedef Ingredients<Config> Ing;

	Ingredients<Config>::molecules_type mol[2];

	double energyState[2];
	double energyWinStart[2];
	double energyWinEnd[2];
	double lnDOSenergyOld[2];
	double lnDOSenergyNew[2];

	uint32_t rnd_tid;
	bool acceptExchange[2];
	//mol[1]= myIngredients.getMolecules();

	int counterCovergedIteration = 0;


#pragma omp parallel private(nthreads, tid) shared(mol, energyState, energyWinStart, energyWinEnd, lnDOSenergyOld, lnDOSenergyNew, rnd_tid, acceptExchange, counterCovergedIteration)
	{
		/* Obtain thread number */
		tid = omp_get_thread_num();

		std::cout << "Hello World from thread = " <<  tid << std::endl;

		/* Only master thread does this */
		if (tid == 0)
		{
			nthreads = omp_get_num_threads();
			std::cout << "Number of threads = " << nthreads << std::endl;
		}

		//seed the globally available random number generators
		RandomNumberGenerators rng;
		rng.seedAll();

		Ing myIngredients;


		myIngredients.modifyVisitsEnergyStates().reset(min_histogram, max_histogram, bins_histogram);//-128*6*4,1,4*6*4*4*8*2*256);
		myIngredients.modifyTotalVisitsEnergyStates().reset(min_histogram, max_histogram, bins_histogram);
		//myIngredients.modifyHGLnDOS().reset(min_histogram, max_histogram, bins_histogram);

		//read-in the histogram of the DOS
		if(HGLnDOSfile.empty())
		{
			myIngredients.modifyHGLnDOS().reset(min_histogram, max_histogram, bins_histogram);
			throw std::runtime_error("File HGLnDOS has to be provided. EXITing...\n");
		}
		else
		{
			std::cout << "ReadIn of HGLnDOS:        " << HGLnDOSfile << std::endl;
			ReadInHGLnDOS in(min_histogram, max_histogram, bins_histogram, HGLnDOSfile);
			in.readin();

			//copy HGLnDOS into ingredients
			myIngredients.modifyHGLnDOS().reset(min_histogram, max_histogram, bins_histogram);

			for(size_t n=0;n<in.getHGLnDOS().getVectorValues().size();n++){

				if(in.getHGLnDOS().getVectorValues()[n].ReturnN() != 0)
					myIngredients.modifyHGLnDOS().addValue(in.getHGLnDOS().getCenterOfBin(n), in.getHGLnDOS().getVectorValues()[n].ReturnM1());
			}
		}

		// run the simulation and gather the information

		/* if lengthIncrease == 0
		double lengthWindow = (maxWin-minWin)/(omp_get_num_threads()-(omp_get_num_threads()-1)*overlap);
		double minWinThread = lengthWindow*tid*(1.0-overlap)+minWin;
		double maxWinThread = minWinThread+lengthWindow;
		*/

		double lengthFirst=(maxWin-minWin)/(overlap+(1.0-overlap)*(std::pow((1.0+lengthIncrease), 1.0*omp_get_num_threads())-1.0)/lengthIncrease );

		double minWinThread = minWin+lengthFirst*(std::pow((1.0+lengthIncrease), tid)-1.0)/lengthIncrease  - overlap*lengthFirst*(std::pow((1.0+lengthIncrease), (tid+1.0))-(1.0+lengthIncrease))/lengthIncrease;
		double maxWinThread = minWinThread+lengthFirst*std::pow((1.0+lengthIncrease), 1.0*tid);

		// reorder in case the boundaries are used in wrong order
		if(maxWinThread < minWinThread)
		{
			double tmp = minWinThread;
			minWinThread = maxWinThread;
			maxWinThread = tmp;
		}


		if(readinBFMinWin == false)
		{
			UpdaterReadBfmFile<Ing> UR(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE);
			UR.initialize();
			UR.execute();
			UR.cleanup();
		}
		else
		{
			std::stringstream ssprefixWindowBFM;
			ssprefixWindowBFM << "_idxWin" << std::setw(2) << std::setfill('0') << tid;

			UpdaterReadBfmFile<Ing> UR(std::string(infile + ssprefixWindowBFM.str() + ".bfm"),myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE);
			UR.initialize();
			UR.execute();
			UR.cleanup();

			myIngredients.setName(infile);
		}

		UpdaterAdaptiveWangLandauSamplingNextNeighbor<Ing,MoveLocalSc> UWL(myIngredients,
						save_interval,
						bias_update_interval, flatness, modFactor, max_mcs, minWinThread, maxWinThread, tid);




		UWL.initialize();

		std::cout << "Thread " << tid << " with window [ " <<   minWinThread << " ; " << maxWinThread << " ] " << std::endl;

		// if some BFM-files are read in its highly likely that they already in window
		if(readinBFMinWin == true)
		{
			UWL.checkForFirstWindowAppearance();
		}
		// run as long to reach desired window
		do
		{
			UWL.execute();
		} while(!myIngredients.isEnergyInWindow());

#pragma omp barrier
		// here all threads are sync and in therir desried window
		// run the simulations

		do
		{
#pragma omp barrier
			// run the one iterartion until histogram converged

			do
			{
#pragma omp barrier
				// output configuration
				/*if(filedump)
				{
					std::stringstream ss;
					ss << infile << "_" << tid << ".bfm";

					AnalyzerWriteBfmFile<Ing> ABFM(ss.str(),myIngredients);
					ABFM.initialize();
					ABFM.execute();
					ABFM.cleanup();
				}*/

				//for(int count = 0; count < 1; count++)
				UWL.execute();

				#pragma omp barrier
				//exchange configurations


				#pragma omp single
					{
						rnd_tid=rng.r250_rand32()%(nthreads-1);
				#pragma  omp flush(rnd_tid)
					}

				#pragma omp barrier
					if (tid == rnd_tid)
					{
						mol[0]= myIngredients.getMolecules();
						energyState[0] =  myIngredients.getInternalEnergyCurrentConfiguration(myIngredients);
						energyWinStart[0]= myIngredients.getMinWin();
						energyWinEnd[0]= myIngredients.getMaxWin();
						lnDOSenergyOld[0] = myIngredients.getHGLnDOS().getCountAt(energyState[0]);
						//lnDOSenergyNew
						std::cout << "copy molecules tid " <<  tid << " with energy " << energyState[0] << " and lnDOS " << lnDOSenergyOld[0] << " in win [ " << energyWinStart[0] << " ; " << energyWinEnd[0] << " ] " << std::endl;
					#pragma  omp flush
					}
					if (tid == rnd_tid+1)
					{
						mol[1]= myIngredients.getMolecules();
						energyState[1] =  myIngredients.getInternalEnergyCurrentConfiguration(myIngredients);
						energyWinStart[1]= myIngredients.getMinWin();
						energyWinEnd[1]= myIngredients.getMaxWin();
						lnDOSenergyOld[1] = myIngredients.getHGLnDOS().getCountAt(energyState[1]);
						std::cout << "copy molecules tid " <<  tid << " with energy " << energyState[1] << " and lnDOS " << lnDOSenergyOld[1] << " in win [ " << energyWinStart[1] << " ; " << energyWinEnd[1] << " ] " << std::endl;
					#pragma  omp flush
					}
				#pragma omp barrier
					if (tid == rnd_tid)
					{
						lnDOSenergyNew[0] = myIngredients.getHGLnDOS().getCountAt(energyState[1]);
						acceptExchange[0] = false;
					#pragma  omp flush
					}

					if (tid == rnd_tid+1)
					{
						lnDOSenergyNew[1] = myIngredients.getHGLnDOS().getCountAt(energyState[0]);
						acceptExchange[1] = false;
					#pragma  omp flush
					}
				#pragma omp barrier
					if (tid == rnd_tid)
					{

						// check boundary windows
						// only one decision for both
						if((energyState[1] < energyWinEnd[0]) && (energyState[0] > energyWinStart[1]) )
							if((energyState[0] < energyWinEnd[1]) && (energyState[1] > energyWinStart[0]) )
						{
							double diffLnDOS = lnDOSenergyOld[0]+lnDOSenergyOld[1]-lnDOSenergyNew[0]-lnDOSenergyNew[1];

							double p = 1.0;
							if(diffLnDOS < 0.0)
								p=std::exp(diffLnDOS);

							if(rng.r250_drand() < p){
								acceptExchange[0] = true;
								acceptExchange[1] = true;
							}
							else
							{
								//reject exchange
							}
						
						#pragma  omp flush
						}
					}

					/*if (tid == rnd_tid+1)
					{
						if((energyState[1] < energyWinEnd[0]) && (energyState[0] > energyWinStart[1]) )
						{
							double diffLnDOS = lnDOSenergyOld[0]+lnDOSenergyOld[1]-lnDOSenergyNew[0]-lnDOSenergyNew[1];

							double p = 1.0;
							if(diffLnDOS < 0.0)
								p=std::exp(diffLnDOS);

							if(rng.r250_drand() < p)
								acceptExchange[1] = true;
						
						#pragma  omp flush
						}
					}
					*/
				#pragma omp barrier
					if (tid == rnd_tid)
					{
						if((acceptExchange[0] == true) && (acceptExchange[1] == true) )
						{
							std::cout  << std::endl << std::endl << std::endl << "swap molecules tid " <<  tid << std::endl << std::endl << std::endl;
							myIngredients.modifyMolecules() = mol[1];
							myIngredients.synchronize();
						#pragma  omp flush
						}
						else
						{
							myIngredients.rejectMove(myIngredients);
							#pragma  omp flush
						}
					}
					if (tid == rnd_tid+1)
					{
						if( (acceptExchange[0] == true) && (acceptExchange[1] == true) )
						{
							std::cout << std::endl << std::endl << "swap molecules tid " <<  tid << std::endl << std::endl << std::endl;
							myIngredients.modifyMolecules() = mol[0];
							myIngredients.synchronize();
						#pragma  omp flush
						}
						else
						{
							myIngredients.rejectMove(myIngredients);
							#pragma  omp flush
						}
					}

				#pragma omp barrier



				#pragma omp barrier

				//UWL.histogramConverged();
				//counter++;
				if(UWL.histogramConverged() && UWL.isFirstConverged())
				{
					UWL.outputConvergedIteration();

					UWL.unsetFirstConverged();

				#pragma omp atomic
					counterCovergedIteration = counterCovergedIteration+1;
				}
				#pragma  omp flush(counterCovergedIteration)
				#pragma  omp flush

				std::cout << std::endl << std::endl << "tid" <<  tid << " counterCovergedIteration -> " << counterCovergedIteration  << " / " << omp_get_num_threads() <<  std::endl << std::endl << std::endl;

				#pragma omp barrier

			} while(counterCovergedIteration != omp_get_num_threads());//!UWL.histogramConverged());
			//iteration coverged
			#pragma omp barrier

			#pragma omp single
					{
						counterCovergedIteration=0;
			#pragma  omp flush(counterCovergedIteration)
					}
			#pragma omp barrier
			#pragma  omp flush


			//run as long for each interation
			#pragma omp barrier
			//std::cout << std::endl << std::endl << "tid" <<  tid << " NextIterStart: counterCovergedIteration -> " << counterCovergedIteration  << " / " << omp_get_num_threads() <<  std::endl << std::endl << std::endl;


			//reset for new iteration
			UWL.doResetForNextIteration();

		}while( !(myIngredients.getModificationFactor() < std::exp(std::pow(10,-8)) ) );

		// all iteration converged and f < exp(10-8)
		#pragma omp barrier

		UWL.cleanup();



		/*
		TaskManager taskmanager;
		taskmanager.addUpdater(new UpdaterReadBfmFile<Ing>(infile,myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
		//here you can choose to use MoveLocalBcc instead. Careful though: no real tests made yet
		//(other than for latticeOccupation, valid bonds, frozen monomers...)
		//taskmanager.addUpdater(new UpdaterSimpleSimulator<Ing,MoveLocalSc>(myIngredients,save_interval));

		taskmanager.addUpdater(new UpdaterAdaptiveWangLandauSamplingNextNeighbor<Ing,MoveLocalSc>(myIngredients,
				save_interval,
				bias_update_interval, modFactor, max_mcs, minWinThread, maxWinThread),
				1);

		if (tid == 0){
			if(filedump)
				taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile,myIngredients));
		}

		taskmanager.initialize();
		taskmanager.run();
		taskmanager.cleanup();

		*/

		//if (tid == 0)
		{
			std::stringstream ssprefixWindow;
			ssprefixWindow << "_idxWin" << std::setw(2) << std::setfill('0') << tid;
			outfile=infile+ssprefixWindow.str()+"_final";

			AnalyzerWriteBfmFile<Ing> ABFM(outfile,myIngredients);
			ABFM.initialize();
			ABFM.execute();
			ABFM.cleanup();
		}

	}

	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

