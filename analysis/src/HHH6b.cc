/////////////////////////////////////////////////////////////////////////////////////////////

// C++
#include <fstream>
#include <random>
#include <functional>

#include "run.h"
#include "analysis.h"
#include "samples.h"
#include "detector.h"

#include "HepMC/IO_Exception.h"

using namespace Pythia8;

void printEvent(finalState ev){
    // For debugging
    for (int j = 0; j < ev.size(); j++){
        std::cout << "Jet " << j << ":" << std::endl;
        std::cout << "PDG: " << ev[j].user_index() << std::endl;
        for (int m = 0; m < ev[j].four_mom().size(); m++){
            std::cout << "p_" << m << ": " << ev[j].four_mom()[m] << std::endl;
            }
        }
    }

////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::uint32_t> initSeeds( runCard const& run, sampleCard const& sample, int const& subsample )
{
  const uint32_t primary = run.runseed;
  const uint32_t samplename = std::hash<std::string>()(sample.samplename);
  const uint32_t runname = std::hash<std::string>()(run.runname);
  const uint32_t subsampleID = subsample;

  std::seed_seq seq{primary,samplename,runname,subsampleID}; 
  std::vector<std::uint32_t> seeds(4);
  seq.generate(seeds.begin(), seeds.end());
  return seeds;
}
 
////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] ) 
{  
  if (argc != 4)
  {
    cerr << "Error: Wrong number of arguments!"<<endl;
    cerr << "Usage: HHH6b <run card> <sample card> <subsample>" <<endl;
    exit(-1);
  }

  // Read run card
  const std::string runfile = std::string(argv[1]);
  const runCard run(runfile);

  // Read sample card 
  const std::string samplefile = std::string(argv[2]);
  const sampleCard sample(samplefile);

  // Determine subsample constants
  const int subsample = atoi(argv[3]);
  const int sampleStart = subsample*run.sub_samplesize; // start point of the subsample

  const std::vector<std::uint32_t> seeds = initSeeds(run, sample, subsample);
  const uint32_t pythiaSeed   = ((double)seeds[0]/pow(2,32))*9E8; // Pythia seeds must be < 9E8
  const uint32_t pileupSeed   = ((double)seeds[1]/pow(2,32))*9E8; // Pythia seeds must be < 9E8
  const uint32_t detectorSeed = seeds[2];
  const uint32_t analysisSeed = seeds[3];

  cout << "Processing sample: " <<sample.samplename<< ", subsample: "<<subsample <<std::endl;
  cout  <<"  RNG Seeds - Shower:   "<<pythiaSeed  <<std::endl
        <<"            - PU:       "<<pileupSeed <<std::endl
        <<"            - Detector: "<<detectorSeed <<std::endl
        <<"            - Analysis: "<<analysisSeed <<std::endl;

  // Initialise Pythia and HepMC
  Pythia pythiaRun(std::string(PYTHIADIR)); // Pythia input
  std::ifstream hepmc_is( sample.eventpath.c_str() );                   // HepMC input

  // Initialise the event sample and weight normalisation
  double weight_norm = 0;
  if (!sample.hepmc) InitPythia(run, sample, pythiaSeed, pythiaRun, weight_norm );
  else InitHepMC( run, sample, weight_norm);
  weight_norm *= 1000*sample.xsec_norm; // Includes pb->fb conversion

  // Initialse Analyses and detector simulation
  vector<Analysis*> analyses;
  InitAnalyses(analyses, run, sample, subsample);
  Detector detector(run, sample, pileupSeed, detectorSeed);

  // Skip to subsample x
  cout << "Skipping to startpoint: " << sampleStart <<endl;
  for (int iEvent = 0; iEvent < sampleStart; ++iEvent) 
  {
    try{
        double dum; finalState dum2;
        if (!sample.hepmc) get_final_state_particles(pythiaRun, dum2, dum);
            else get_final_state_particles(hepmc_is,  dum2, dum);
        }
    catch(int ex){
        // Bad values
        continue;
        }
  }
    
  ////////////////////////////////////////////////////////////////////////////////
  //                                                                            //  
  //                                EVENT LOOP                                  //
  //                                                                            //
  ////////////////////////////////////////////////////////////////////////////////
  
  cout << "*************** Analysis Begins ***************" <<endl;
  const int targetSize = min(run.sub_samplesize, sample.nevt_sample - sampleStart);
  cout << "Analysing: " << targetSize <<" events"<<endl;
  int negev = 0;
  int nanev = 0;
  for (int iEvent = 0; iEvent < targetSize; ++iEvent) 
  {
    finalState ifs, fs; // The event final state
    
    double event_weight=0;
    
    if (!sample.hepmc) get_final_state_particles(pythiaRun, ifs, event_weight);
    else try{
        get_final_state_particles(hepmc_is,  ifs, event_weight);
    }
    catch (HepMC::IO_Exception &e) { // nan values appearing - disregard event
        std::cout << "WARNING: invalid kinematical quantity (NaN), event discarded" << iEvent << std::endl;
        nanev++;
        continue;
    }
    catch(int ex){
        std::cout << "WARNING: negative invariant mass squared, event discarded" << iEvent << std::endl;
        negev++;
        continue;
        }
        
    // Perform detector simulation
    detector.Simulate(ifs,fs);
        
    // Run over analyses
    event_weight *= weight_norm;
    
    if (ifs.size() != 0){
        for (size_t i=0; i<analyses.size(); i++){
          analyses[i]->Analyse(sample.is_signal, event_weight, fs);
        }
 //       if ((iEvent+1) % 100 == 0 )
          cout << iEvent+1 <<" events analysed"<<endl;
    }
  }

    if (negev > 0) std::cout << "WARNING: discarded " << negev << " events due to negative invariant masses." << std::endl;
    if (nanev > 0) std::cout << "WARNING: discarded " << nanev << " events due to a NaN kinematic quantity." << std::endl;
    
  // Clean up
  for (size_t i=0; i<analyses.size(); i++)
    delete analyses[i];
  hepmc_is.close();

  // End of the main program
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
