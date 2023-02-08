

//#include <omp.h>

#include <cstring>
#include <sstream>      // std::stringstream, std::stringbuf
#include <cmath>

#include <cstdio>      //sprintf
#include <climits>     // INT_MAX

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

#include "FeatureWangLandauExtendedShellInteraction.h"
#include "UpdaterAdaptiveWangLandauSamplingNextNeighbor.h"
#include "ReadInHGLnDOS.h"

#include "mpi.h" // Include MPI header file containing the librarys API

int main(int argc, char* argv[])
{
	// Start up MPI
	MPI_Init (&argc, &argv);
	
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

		uint32_t min_statistic_entries=100;
		
		int numWalkerPerWindow = 2;

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
		| clara::Opt( [&min_statistic_entries](int const b)
			{
				if (b < 0)
				{
					return clara::ParserResult::runtimeError("min_statistic_entries must be greater than 0");
				}
				else
				{
					min_statistic_entries = b;
					return clara::ParserResult::ok(clara::ParseResultType::Matched);
				}
			}, "min_statistic_entries (=1000)")
			["--min-statistic"]
			("(required) <integer> Number of entries for calculating the average at specific energy." )
			.required()
		| clara::Opt(  numWalkerPerWindow, "walker per energy window (=2)" )
			["--walker"]
			("every energy window has this number of walker (=2)")
			.required()
		| clara::Help( showHelp );

		auto result = parser.parse( clara::Args( argc, argv ) );
		if( !result ) {
			std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
			exit(1);
		}
		else if(showHelp == true)
		{
			std::cout << "Simulator for the ScBFM with Ex.Vol and BondCheck and WangLandau in ES-shell" << std::endl
					<< "maximum number of connections per monomer is 8" << std::endl
					<< "Features used: FeatureBondset, FeatureAttributes, FeatureWangLandauExtendedShellInteraction<FeatureLattice<uint8_t> >" << std::endl
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
					<< "walker per energy window: " << numWalkerPerWindow << std::endl
					;
		}
	
	// ******************************************
	// Here goes the MPI stuff
	// ******************************************
	
	// See: DOI 10.1088/1742-6596/1012/1/012003
	
	// File handlers for I/O
	//FILE *file;
	FILE *stdoutlog;
	//FILE *wanderlog;
	//char filename[50];
	char stdoutlogname[128];
	
	// Set up local MPI communicators for replica-exchange (RE)
	// - each process belongs to two local groups/communicators
	// - each process has different local IDs in different communicators generally
	MPI_Status status;
	
	int numprocs; // total number of walkers
	// Get total number of processes
	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
	
	int myid; // my rank(ID) in the global communicator, MPI_COMM_WORLD
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	
	int multiple = numWalkerPerWindow;//2; // number of walkers having the same energy window
	
	// Error log for every process
	sprintf(stdoutlogname, "Error%04i.log", myid);
	
	// at the moment, the code works only for an _odd_ number of energy windows
	// (to make the RE in windows at the outside consistent)
	if ((numprocs/multiple)%2 == 0)
    {
      if (myid == 0) 
        {
			stdoutlog=fopen(stdoutlogname,"a");
			std::cerr << "ERROR: Even number of energy windows " << int(numprocs/multiple) << " requested. Please request an odd number of energy windows." << std::endl<< std::endl;
			fprintf(stdoutlog, "ERROR: Even number of energy windows (%d) requested. Please request an odd number of energy windows.\n\n", numprocs/multiple);
			fclose(stdoutlog);
          
        }

      MPI_Abort(MPI_COMM_WORLD, 1);
    }
	
	
	int comm_id; // ID for a communicator
	int mylocalid[2]; // my ID in local communicators
	
	// Get the group of processes in MPI_COMM_WORLD (i.e., all processors)
	MPI_Group world;
	MPI_Comm_group(MPI_COMM_WORLD, &world);
	
	// The followings are for defining a list of MPI local communicators
	stdoutlog=fopen(stdoutlogname,"a");
	int *ranks;  // an array to store global IDs
	ranks = (int*) malloc(2*multiple*sizeof(int)); // for walkers in a local comm.
	
	int numLocalComm = (numprocs/multiple) - 1; // number of local communicators
	
	MPI_Group *mpi_local_group; // an array to store local groups
	MPI_Comm *mpi_local_comm;  // an array to store local comm.
	
	mpi_local_group = (MPI_Group*) malloc(numLocalComm*sizeof(MPI_Group));
	mpi_local_comm = (MPI_Comm*) malloc(numLocalComm*sizeof(MPI_Comm));
	
	for (int i=0; i<numLocalComm; i++) // i: counter for local communicators
	{
		// For each walker in local communicator (group) i, calculate the global ID
		// (its rank in MPI_COMM_WORLD) and put them in the ‘ranks’ array
		
		for (int j=0; j<2*multiple; j++) // j: counter for walkers in a local comm.
		{
			ranks[j] = i*multiple+j;
			
			if (myid==0) 
			{
			fprintf(stdoutlog,"Proc %3i: %i will be part of communicator/group %i\n",myid,ranks[j],i);
			}
		}
		
		MPI_Group_incl(world,2*multiple,ranks,&mpi_local_group[i]); // create local group
		MPI_Comm_create(MPI_COMM_WORLD,mpi_local_group[i],&mpi_local_comm[i]); // create communicator for that group
	}
	
	free(ranks);
	
	// get my local id (in my local communicators)
  if (myid<numprocs-multiple)  // Every processor except those in the last window
    {
      comm_id=2*(myid/(2*multiple));
      MPI_Comm_rank(mpi_local_comm[comm_id], &mylocalid[0]);
      fprintf(stdoutlog,"Proc %3i: I am part of communicator/group %i with local_id[0]=%i\n",myid,comm_id,mylocalid[0]);
    }
  else
    {
      mylocalid[0]=INT_MAX; // just to give it a value
      fprintf(stdoutlog,"Proc %3i: got local_id[0]=%i\n",myid,mylocalid[0]);
    }
    
	if (myid>=multiple)
    {
      comm_id=2*((myid-multiple)/(2*multiple))+1;
      MPI_Comm_rank(mpi_local_comm[comm_id], &mylocalid[1]);
      fprintf(stdoutlog,"Proc %3i: I am part of communicator/group %i with local_id[1]=%i\n",myid,comm_id,mylocalid[1]);
    }
  else
    {
      mylocalid[1]=INT_MAX; // just to give it a value
      fprintf(stdoutlog,"Proc %3i: got local_id[1]=%i\n",myid,mylocalid[1]);
    }

    fprintf(stdoutlog,"Proc %3i: Start WL iteration\n",myid);
	fclose(stdoutlog);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
	
	int nthreads, tid;

	// FeatureExcludedVolume<> is equivalent to FeatureExcludedVolume<FeatureLattice<bool> >
	//typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureFixedMonomers,FeatureAttributes,FeatureExcludedVolumeSc<>) Features;
	//typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes, FeatureNNInteractionSc< FeatureLattice >) Features;
	//typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes, FeatureWangLandauNextNeighbor< FeatureLattice >) Features;
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes< >, FeatureWangLandauExtendedShellInteraction< FeatureLattice >) Features;

	//typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureWangLandau, FeatureAttributes) Features;

	typedef ConfigureSystem<VectorInt3,Features, 8> Config;
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


//#pragma omp parallel private(nthreads, tid) shared(mol, energyState, energyWinStart, energyWinEnd, lnDOSenergyOld, lnDOSenergyNew, rnd_tid, acceptExchange, counterCovergedIteration)
	{
		/* Obtain thread number */
		/* Obtain thread number */
		tid = myid;//omp_get_thread_num();

		stdoutlog=fopen(stdoutlogname,"a");
		fprintf(stdoutlog,"Hello World from thread = %3i out of %3i\n",myid, numprocs);
		fclose(stdoutlog);

		/* Only master thread does this */
		if (tid == 0)
		{
			stdoutlog=fopen(stdoutlogname,"a");
			nthreads = numprocs;//omp_get_num_threads();
			fprintf(stdoutlog,"Number of all processes = %3i\n",nthreads);
			fclose(stdoutlog);
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

			int counter = 0;
			for(size_t n=0;n<in.getHGLnDOS().getVectorValues().size();n++){

				if(in.getHGLnDOS().getVectorValues()[n].ReturnN() != 0)
				{
					myIngredients.modifyHGLnDOS().addValue(in.getHGLnDOS().getCenterOfBin(n), in.getHGLnDOS().getVectorValues()[n].ReturnM1());
					counter++;
				}
			}

			// check if something has been read-in
			if(counter == 0)
				throw std::runtime_error("File HGLnDOS is empty. EXITing...\n");
		}

		// run the simulation and gather the information

		/* if lengthIncrease == 0
		double lengthWindow = (maxWin-minWin)/(omp_get_num_threads()-(omp_get_num_threads()-1)*overlap);
		double minWinThread = lengthWindow*tid*(1.0-overlap)+minWin;
		double maxWinThread = minWinThread+lengthWindow;
		*/

		//double lengthFirst=(maxWin-minWin)/(overlap+(1.0-overlap)*(std::pow((1.0+lengthIncrease), 1.0*omp_get_num_threads())-1.0)/lengthIncrease );
		double lengthFirst=(maxWin-minWin)/(overlap+(1.0-overlap)*(std::pow((1.0+lengthIncrease), 1.0*(numprocs/multiple))-1.0)/lengthIncrease );

		int tid_window = tid/multiple;
		double minWinThread = minWin+lengthFirst*(std::pow((1.0+lengthIncrease), tid_window)-1.0)/lengthIncrease  - overlap*lengthFirst*(std::pow((1.0+lengthIncrease), (tid_window+1.0))-(1.0+lengthIncrease))/lengthIncrease;
		double maxWinThread = minWinThread+lengthFirst*std::pow((1.0+lengthIncrease), 1.0*tid_window);

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
			ssprefixWindowBFM << "_idxWin" << std::setw(4) << std::setfill('0') << tid;

			UpdaterReadBfmFile<Ing> UR(std::string(infile + ssprefixWindowBFM.str() + ".bfm"),myIngredients,UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE);
			UR.initialize();
			UR.execute();
			UR.cleanup();

			myIngredients.setName(infile);
			myIngredients.modifyMolecules().setAge(0); // reset the clock within the file
		}

		UpdaterAdaptiveWangLandauSamplingNextNeighbor<Ing,MoveLocalSc> UWL(myIngredients,
						save_interval,
						bias_update_interval, flatness, modFactor, max_mcs, minWinThread, maxWinThread, tid, min_statistic_entries);




		UWL.initialize();

		std::cout << "Thread " << tid << " with window [ " <<   minWinThread << " ; " << maxWinThread << " ] " << std::endl;

		// if some BFM-files are read in its highly likely that they already in window
		//if(readinBFMinWin == true)
		{
			UWL.checkForFirstWindowAppearance();
		}
		// run as long to reach desired window
		do
		{
			UWL.execute();
		} while(!myIngredients.isEnergyInWindow());

		myIngredients.modifyMolecules().setAge(0); // reset the clock within the file
		
		 MPI_Barrier(MPI_COMM_WORLD);
		 
		 int tryleft=0, tryright=0, exchangeleft=0, exchangeright=0;
		 
		 
		int hasFullStatisics_slowest = 0; // the slowest simulation determines the convergence
		int hasFullStatisics_recent  = 0; // the slowest simulation determines the convergence
		 
		 
		 do
		{
			//std::cout << "rsync all threads iterartion " << myid << std::endl;
//MPI_Barrier(MPI_COMM_WORLD);
// #pragma omp barrier
			// run the one iterartion until histogram converged

			
			
			//do // run until all walker WITHIN ONE windows converged for their iteration
			{
				//std::cout << "rsync all threads windows conversion " << myid << std::endl;
MPI_Barrier(MPI_COMM_WORLD);
//#pragma omp barrier
				// output configuration
				///if(filedump)
				//{
					//std::stringstream ss;
					//ss << infile << "_" << tid << ".bfm";

					//AnalyzerWriteBfmFile<Ing> ABFM(ss.str(),myIngredients);
					//ABFM.initialize();
					//ABFM.execute();
					//ABFM.cleanup();
				//}/

				//for(int count = 0; count < 1; count++)
				UWL.execute();

				//#pragma omp barrier
				// workaround for broadcasting the swap direction to all processes
				MPI_Barrier(MPI_COMM_WORLD);
				int swap_direction=rng.r250_rand32()%2;
				MPI_Bcast(&swap_direction, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				
				//exchange configurations

				// THIS IS THE MASTER RE / SWAP FUNCTION
				//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				//int replica_exchange(int swap_direction, int index_akt)
				{
					//energyState[0] =  myIngredients.getInternalEnergyCurrentConfiguration(myIngredients);
					//energyWinStart[0]= myIngredients.getMinWin();
					//energyWinEnd[0]= myIngredients.getMaxWin();
					//lnDOSenergyOld[0] = myIngredients.getHGLnDOS().getCountAt(energyState[0]);
					
					double E_new; // energy by other process
					double E_old; // current energy
					
					
					// frac refers to local exchange probability
					// wk is combined exchange probability
					double myfrac,otherfrac,randx,wk;
					
					int change=0; // boolean: 0 for not exchanging, 1 for exchanging
					int swap_partner=-1; // id of swap partner (receive_buffer)
					
					// 0 -> to the right; 1 -> to the left
					/*int swap_direction=rng.r250_rand32()%2;*///(myIngredients.getMolecules().getAge()/save_interval)%2;// swap_direction%2; // comes actually as number of swap attempt
					
					// everyone has to find its swap-partner
					
					int *pairs; // array containing the partners of each process (send_buffer)
					pairs = (int*)malloc(2*multiple*sizeof(int));
					
					if (mylocalid[swap_direction]==0) // 'head-node' in the energy window determines pairs of flippartners
					{
						int choose_from=multiple; // number of free partners in higher window of communicator
						int select; // storage for random number
						
						int *libre; // list of free partners from higher window in communicator
						libre = (int*)malloc(multiple*sizeof(int));
						
						for (int i=0;i<multiple;i++) libre[i]=multiple+i; // initialise
						
						// idea: processes from the lower window choose someone from the higher window at random
						// of course, the chosen walker can't have an exchange partner yet
						for (int i=0;i<multiple;i++) // loop over processes in the lower window
						{
							select=rng.r250_rand32()%choose_from;
							pairs[i]=libre[select];
							pairs[libre[select]]=i; // the 'vice-versa pair'
							// update list
							choose_from--;
							for (int j=select;j<choose_from;j++)
								libre[j]=libre[j+1];
						}
						
						//{
							//stdoutlog=fopen(stdoutlogname,"a");
							//fprintf(stdoutlog,"Proc %3i: Drew the following swap partners:\n",myid);
							//for (int i=0;i<2*multiple;i++)
								//fprintf(stdoutlog,"Proc %3i: %i -- %i (local ids in communicator)\n",myid,i,pairs[i]);
							
							//fprintf(stdoutlog,"Proc %3i: tryleft: %i, exchangeleft %i (Akzeptanzleft:%.2lf) <--> tryright: %i, exchangeright %i (Akzeptanzright:%.2lf)\n",myid,tryleft,exchangeleft,(double)exchangeleft/(double)tryleft,tryright,exchangeright,(double)exchangeright/(double)tryright);						      fclose(stdoutlog);
						//}
						
						free(libre);
					}
					
					// at this point, every walker has a swap partner assigned, now they must be communicated
					if ((swap_direction==0)&&(myid<(numprocs-multiple))) // the walkers from the last node should not swap
					{
						comm_id=2*(myid/(2*multiple)); // ! all integer, the '/' is a (div) ! Not the same as myid/multiple !
						MPI_Scatter(pairs,1,MPI_INT,&swap_partner,1,MPI_INT,0,mpi_local_comm[comm_id]);
					}
					
					if ((swap_direction==1)&&(myid>=multiple)) // the walkers from the zero-node should not swap
					{
						comm_id=((myid-multiple)/(2*multiple))*2+1; // ! all integer, the '/' is a (div) ! See above
						MPI_Scatter(pairs,1,MPI_INT,&swap_partner,1,MPI_INT,0,mpi_local_comm[comm_id]);
					}
					
					std::cout << "rsync all threads " << myid << std::endl;
					MPI_Barrier(MPI_COMM_WORLD);
					
					free(pairs);
					
					if (swap_partner!=-1) // i.e. if there is a swap-partner for me (if I am at a boundary, I might not have a swap partner this time)
					{
						// statistics
						if (swap_partner>mylocalid[swap_direction]) tryright++;
						else tryleft++;
						
						// safety cross check
						E_old=myIngredients.getInternalEnergyCurrentConfiguration(myIngredients);
						
						
						E_new = myIngredients.getInternalEnergyCurrentConfiguration(myIngredients);
						//i_new= myIngredients.getHGLnDOS().getBinNo(E_new);
						
						
						// get histogram index from my swap partner
						//MPI_Sendrecv_replace(&i_new,1,MPI_INT,swap_partner,1,swap_partner,1,mpi_local_comm[comm_id],&status);
						MPI_Sendrecv_replace(&E_new,1,MPI_DOUBLE,swap_partner,1,swap_partner,1,mpi_local_comm[comm_id],&status);
						
						//if (Ecur+(2*numberspins)!=index_akt)
						//{
							//stdoutlog=fopen(stdoutlogname,"a");
							//fprintf(stdoutlog,"Proc %3i, replica_exchange(): E_old=%f with lnDOS=%f, E_new=%f with lnDOS=%f from %3i on swap direction %i at time %i.\n",myid,E_old,myIngredients.getHGLnDOS().getCountAt(E_old), E_new, myIngredients.getHGLnDOS().getCountAt(E_new), swap_partner, swap_direction, myIngredients.getMolecules().getAge());
							//fclose(stdoutlog);
							////MPI_Abort(MPI_COMM_WORLD,1);
						//}
						
						
						if ((E_new>myIngredients.getMaxWin())||(E_new<myIngredients.getMinWin())) // energyranges must overlap!
						{
							myfrac=-1.0;
						}
						else
						{
							// calculate my part of the exchange probability
							//myfrac=exp(lngE[index_akt]-lngE[i_new]); // g(myE)/g(otherE)
							
							myfrac=std::exp(myIngredients.getHGLnDOS().getCountAt(E_old) - myIngredients.getHGLnDOS().getCountAt(E_new)); //lnDOSenergyOld[0]+lnDOSenergyOld[1]-lnDOSenergyNew[0]-lnDOSenergyNew[1];

						}
						
						if (mylocalid[swap_direction]<multiple) // I am receiver and calculator
						{
							// get my partners part of the exchange probability
							MPI_Recv(&otherfrac,1,MPI_DOUBLE,swap_partner,2,mpi_local_comm[comm_id],&status);
							
							// calculate combined exchange probability and roll the dice
							if ((myfrac>0.0)&&(otherfrac>0.0))
							{
								//randx=(1.0*rng.r250_drand()/(RAND_MAX+1.0));
								randx=rng.r250_drand();
								wk=myfrac*otherfrac;
								if (randx<wk) change=1;
							}
							
							// tell my swap partner whether to exchange or not
							MPI_Send(&change,1,MPI_INT,swap_partner,3,mpi_local_comm[comm_id]);
						}
						else // I just send my part of exchange probability and await decision
						{
							MPI_Send(&myfrac,1,MPI_DOUBLE,swap_partner,2,mpi_local_comm[comm_id]);
							MPI_Recv(&change,1,MPI_INT,swap_partner,3,mpi_local_comm[comm_id],&status);
						}
						
						// if decision was made to exchange configurations
						if (change==1)
						{
							int *tmp_coordinates; // array containing the molecules coordinates process (send_buffer)
							size_t tmp_coordinates_size = myIngredients.getMolecules().size();
							tmp_coordinates = (int*)malloc(3*tmp_coordinates_size*sizeof(int));
							
							for(int i = 0; i < tmp_coordinates_size; i++)
							{
								tmp_coordinates[i*3+0] =  myIngredients.getMolecules()[i].getX();
								tmp_coordinates[i*3+1] =  myIngredients.getMolecules()[i].getY();
								tmp_coordinates[i*3+2] =  myIngredients.getMolecules()[i].getZ();
							}
							
							/*{
							stdoutlog=fopen(stdoutlogname,"a");
							fprintf(stdoutlog,"Proc %3i, replica_exchange(): E_old=%f with lnDOS=%f, E_new=%f with lnDOS=%f from %3i on swap direction %i at time %i.\n",myid,E_old,myIngredients.getHGLnDOS().getCountAt(E_old), E_new, myIngredients.getHGLnDOS().getCountAt(E_new), swap_partner, swap_direction, myIngredients.getMolecules().getAge());
							fprintf(stdoutlog,"Proc %3i, replica_exchange(): Mono0(%3i,%3i,%3i) and MonoEnd(%3i,%3i,%3i)\n",myid,tmp_coordinates[0], tmp_coordinates[1], tmp_coordinates[2], tmp_coordinates[3*(tmp_coordinates_size-1)], tmp_coordinates[3*(tmp_coordinates_size-1)+1], tmp_coordinates[3*(tmp_coordinates_size-1)+2]);
							fclose(stdoutlog);
							//MPI_Abort(MPI_COMM_WORLD,1);
							}
							*/
							
							// exchange conformations (incl. the 3 'special' polymer)
							//MPI_Sendrecv_replace(&latticepoint[0],numberspins+2+1,MPI_INT,swap_partner,1,swap_partner,1,mpi_local_comm[comm_id],&status);
							MPI_Sendrecv_replace(&tmp_coordinates[0],3*tmp_coordinates_size,MPI_INT,swap_partner,1,swap_partner,1,mpi_local_comm[comm_id],&status);
							
							for(int i = 0; i < tmp_coordinates_size; i++)
							{
								myIngredients.modifyMolecules()[i].setX(tmp_coordinates[i*3+0]);
								myIngredients.modifyMolecules()[i].setY(tmp_coordinates[i*3+1]);
								myIngredients.modifyMolecules()[i].setZ(tmp_coordinates[i*3+2]);
							}
							
							myIngredients.synchronize();
							
							/*{
							stdoutlog=fopen(stdoutlogname,"a");
							fprintf(stdoutlog,"Proc %3i, replica_exchange(): E_old=%f with lnDOS=%f, E_new=%f with lnDOS=%f from %3i on swap direction %i at time %i.\n",myid,E_old,myIngredients.getHGLnDOS().getCountAt(E_old), E_new, myIngredients.getHGLnDOS().getCountAt(E_new), swap_partner, swap_direction, myIngredients.getMolecules().getAge());
							fprintf(stdoutlog,"Proc %3i, replica_exchange(): Mono0(%3i,%3i,%3i) and MonoEnd(%3i,%3i,%3i)\n",myid,tmp_coordinates[0], tmp_coordinates[1], tmp_coordinates[2], tmp_coordinates[3*(tmp_coordinates_size-1)], tmp_coordinates[3*(tmp_coordinates_size-1)+1], tmp_coordinates[3*(tmp_coordinates_size-1)+2]);
							fclose(stdoutlog);
							//MPI_Abort(MPI_COMM_WORLD,1);
							}
							*/
							free(tmp_coordinates);
							
							
							// statistics
							if (swap_partner>mylocalid[swap_direction]) exchangeright++;
							else exchangeleft++;
						}
						//else
						//{
							//myIngredients.rejectMove(myIngredients); // exchange was rejected
						//}
					}
					
					//return(change);

				} // end of replica exchange
				//std::cout << "rsync all threads after RE" << myid << std::endl;
				MPI_Barrier(MPI_COMM_WORLD);
				
				//flat = 0;
				
				
				
				//if(UWL.histogramConverged() && UWL.isFirstConverged())
				/*if(UWL.hasSufficientStatistics())
				{
					UWL.outputConvergedIteration();

					
				}*/
				
				// get the recent decision of sufficient statistics
				hasFullStatisics_recent = UWL.hasSufficientStatistics() ? 1 : 0;
				
					
				// communicate the LOWEST progess to ALL process: 0->not all have good statistics
				MPI_Allreduce(&hasFullStatisics_recent,&hasFullStatisics_slowest,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
					
				
			
				
			} //while(flat != 1);//omp_get_num_threads());//!UWL.histogramConverged());
			
			
			
		// termination for ALL processes (NOT windows) only if ALL process reached the enough statistics
		}while( hasFullStatisics_slowest != 1 );//std::exp(std::pow(10,-8)) ) );
		
//#pragma omp barrier
		// here all threads are sync and in therir desried window
		// run the simulations
		// output the final logDOS
		UWL.cleanup();
		
		{
			std::stringstream ssprefixWindow;
			ssprefixWindow << "_idxWin" << std::setw(4) << std::setfill('0') << tid;
			outfile=infile+ssprefixWindow.str()+"_final";

			AnalyzerWriteBfmFile<Ing> ABFM(outfile,myIngredients);
			ABFM.initialize();
			ABFM.execute();
			ABFM.cleanup();
		}

		

	}

	}
	catch(std::exception& err){std::cerr<<err.what();}
	
	MPI_Finalize(); // Finish up MPI
	
	
	return 0;
  
}


