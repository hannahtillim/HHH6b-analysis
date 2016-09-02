#include "triHiggsAnalysis.h"
#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"

#include "utils.h"
#include "run.h"

#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/Filter.hh"

#include "fastjet/contrib/SoftKiller.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/contrib/Nsubjettiness.hh"

#include <algorithm>

using std::vector;
using namespace fastjet::contrib;
using namespace fastjet;

// LST energy correlations
const EnergyCorrelatorC2 C2(2, EnergyCorrelator::pt_R);
const EnergyCorrelatorD2 D2(2, EnergyCorrelator::pt_R);

// tau2/tau1 NSubjettiness
const NsubjettinessRatio tau21(2,1, KT_Axes(), UnnormalizedMeasure(1));



const int nAnalysis = 4;
const int nCuts = 4;
const std::string aString[nAnalysis] = {"_res", "_interR", "_interB", "_boost"};
const std::string cString[nCuts] = {"_c0", "_c1", "_c2", "_c3"};

triHiggsAnalysis::triHiggsAnalysis(runCard const& run, sampleCard const& sample, int const& subsample) : Analysis("triHiggsAnalysis", run, sample, subsample){

    // Histogram settings
    const int nbins = 50;

    // TriHiggs system
    const double pt_HHH_min = 0;
    const double pt_HHH_max = 1000;
    
    const double m_HHH_min = 0.;
    const double m_HHH_max = 1000.;

    // Individual Higgses
    const double DeltaRmin = 0;
    const double DeltaRmax = 5;

    const double DeltaPhimin = -3.2;
    const double DeltaPhimax = 3.2;

    const double DeltaEtamin = -2.5;
    const double DeltaEtamax = 2.5;
  
    const double m_min = 0.;
    const double m_max = 180.; 
 
    const double pt_min = 0.;
    const double pt_max = 900.;
  
    const double eta_min = -6.;
    const double eta_max = +6.;
  
    const double phi_min = -3.15;
    const double phi_max = +3.15;
    
    for (int j = 0; j < nAnalysis; j++){
        for (int i = 0; i < nCuts; i++){
            const std::string suffix = aString[j] + cString[i];
        
            BookHistogram(new YODA::Histo1D(nbins, pt_HHH_min, pt_HHH_max), "pT_HHH" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, m_HHH_min, m_HHH_max), "mHHH" + suffix);
            
            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H0" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H1" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, pt_min, pt_max), "pt_H2" + suffix);
            
            BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H0" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H1" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, m_min, m_max), "m_H2" + suffix);
            
            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_H0" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_H1" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, eta_min, eta_max), "eta_H2" + suffix);
            
            BookHistogram(new YODA::Histo1D(nbins, 0, 200), "split12_fj1" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, 0, 200), "split12_fj2" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, 0, 200), "split12_fj3" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, 0, 200), "split12_fj" + suffix);  

            BookHistogram(new YODA::Histo1D(nbins, 0, 1), "tau21_fj1" + suffix); 
            BookHistogram(new YODA::Histo1D(nbins, 0, 1), "tau21_fj2" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, 0, 1), "tau21_fj3" + suffix);
            BookHistogram(new YODA::Histo1D(nbins, 0, 1), "tau21_fj" + suffix);      

            BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "C2_fj1" + suffix);  
            BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "C2_fj2" + suffix);
            BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "C2_fj3" + suffix);
            BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "C2_fj" + suffix);   

            BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "D2_fj1" + suffix);  
            BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "D2_fj2" + suffix);
            BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "D2_fj3" + suffix);
            BookHistogram(new YODA::Histo1D(nbins*2, 0, 1), "D2_fj" + suffix); 
        }
    
        //  Bookkeeping
    
        BookHistogram(new YODA::Histo1D(nCuts, 0, nCuts), "CF" + aString[j]);
        
        BookHistogram(new YODA::Histo1D(nbins, 0, 1000), "b0_pT" + aString[j]);
        BookHistogram(new YODA::Histo1D(nbins, 0, 1000), "b1_pT" + aString[j]);
        BookHistogram(new YODA::Histo1D(nbins, 0, 1000), "b2_pT" + aString[j]);
        BookHistogram(new YODA::Histo1D(nbins, 0, 1000), "b3_pT" + aString[j]);
        BookHistogram(new YODA::Histo1D(nbins, 0, 1000), "b4_pT" + aString[j]);
        BookHistogram(new YODA::Histo1D(nbins, 0, 1000), "b5_pT" + aString[j]);
        
        BookHistogram(new YODA::Histo1D(10, 0, 100), "pT_cuts_6j" + aString[j]);
        BookHistogram(new YODA::Histo1D(10, 0, 100), "pT_cuts_6b" + aString[j]);
        
        BookHistogram(new YODA::Histo1D(nbins, 0, 5), "deltaR_02" + aString[j]);
        BookHistogram(new YODA::Histo1D(nbins, 0, 5), "deltaR_12" + aString[j]);
    
    }
    
    BookHistogram(new YODA::Histo1D(nbins, 0, 180), "lR_higgs");
    BookHistogram(new YODA::Histo1D(nbins, 0, 180), "lR_higgs_MW");



    // ********************* Ntuple definition **********************
    
    const std::string tupleSpec = "# signal source weight pt_H0 pt_H1 pt_H2 pt_HHH m_H0 m_H1 m_H2 m_HHH dR_H0H1 dR_H1H2 dR_H0H2 dPhi_H0H1 dPhi_H0H2 dPhi_H1H2 dEta_H0H1 dEta_H1H2 dEta_H0H2 chi_H0H1 chi_H2H1 chi_H0H2 pt_H0_sub0 pt_H0_sub1 pt_H1_sub0 pt_H1_sub1 pt_H2_sub0 pt_H2_sub1";

    const std::string root = "." + GetRoot() + GetSample() + "/";
    std::stringstream suffix; suffix << "." << GetSubSample() <<".dat";

    const std::string resDir = root+"resNTuple"+suffix.str();
    const std::string intRDir = root+"intRNTuple"+suffix.str();
    const std::string intBDir = root+"intBNTuple"+suffix.str();
    const std::string bstDir = root+"bstNTuple"+suffix.str();

    resNTuple.open(resDir.c_str());
    intRNTuple.open(intRDir.c_str());
    intBNTuple.open(intBDir.c_str());
    bstNTuple.open(bstDir.c_str());

    resNTuple << tupleSpec <<std::endl;
    intRNTuple << tupleSpec <<" split12_fj tau21_fj C2_fj D2_fj"<<std::endl;
    intBNTuple << tupleSpec <<" split12_fj1 split12_fj2 tau21_fj1 tau21_fj2 C2_fj1 C2_fj2 D2_fj1 D2_fj2"<<std::endl;
    bstNTuple << tupleSpec <<" split12_fj1 split12_fj2 split12_fj3 tau21_fj1 tau21_fj2 tau21_fj3 C2_fj1 C2_fj2 C2_fj3 D2_fj1 D2_fj2 D2_fj3"<<std::endl;
}

// Check if jets are mass-drop tagged
static vector<PseudoJet> MDtagJets( vector<PseudoJet> const& injets)
{
  // Mass-drop tagger
  double const mu = 0.67;
  double const ycut = 0.09;

  const fastjet::JetDefinition CA10(fastjet::cambridge_algorithm, 1.0);
  const fastjet::MassDropTagger md_tagger(mu, ycut);
  vector<PseudoJet> MDTJets;
  for (size_t i=0; i<injets.size(); i++)
  {
    const fastjet::ClusterSequence cs_sub( injets[i].constituents(), CA10);
    vector<PseudoJet> ca_jets = sorted_by_pt( cs_sub.inclusive_jets() );
    const PseudoJet ca_jet = ca_jets[0];
    const PseudoJet tagged_jet = md_tagger(ca_jet);
    if ( tagged_jet != 0 )
      MDTJets.push_back(injets[i]);
  }
  return MDTJets;
}

// Returns the two hardest GA subjets for each input largeR jet
static vector< vector<PseudoJet> > getSubJets( vector<PseudoJet> const& largeRJets, vector<PseudoJet> const& trackJets )
{
  vector< vector<PseudoJet> > largeRsubJets;
  for ( PseudoJet jet : largeRJets )
  {
    vector<PseudoJet> subJets;
    get_assoc_trkjets( jet, trackJets, subJets, false);
    largeRsubJets.push_back( SelectorNHardest(2)(sorted_by_pt(subJets)) );
  }
  return largeRsubJets;
}

// Small-R B-tagging
// Uses fastjet .user_index() method with PDG numbering system - checks contents directly,
// returns vector of btagTypes
static vector<btagType> BTagging( const vector<PseudoJet>& jets_vec )
{
  vector<btagType> btag_vec;
  for( auto jet : jets_vec )
  {
    btagType type = NTAG;
    const vector<PseudoJet>& jet_constituents = SelectorPtMin(15)(jet.constituents());
    for( auto constituent : jet_constituents )
    {
      const int userid= constituent.user_index();
      if(abs(userid) == 5) type = BTAG;
      if(abs(userid) == 4 && type != BTAG ) type = CTAG;
      if( type == NTAG ) type = LTAG;
    }
    btag_vec.push_back(type);
  }
  return btag_vec;
}

// params: number of tags required, actual contents
static double btagprob(int const& nTag, int const& nB, int const& nC, int const& nL){
    
    const double btag_prob = 1.0;
    const double btag_mistag = 0;
    const double ctag_prob = 0;
    
    double totalprob = 0;
    
    // Loop over all possible classification combinations
    for (int iB=0; iB <= std::min(nTag, nB); iB++){
        for (int iC=0; iC <= std::min(nC, nTag-iB); iC++){
            for (int iL=0; iL <= std::min(nL, nTag-iB-iC); iL++){
                const double bProb = pow(btag_prob, iB) * pow(1.0-btag_prob, nB-iB);
                const double cProb = pow(ctag_prob, iC) * pow(1.0-ctag_prob, nC-iC);
                const double lProb = pow(btag_mistag, iL) * pow(1.0-btag_mistag, nL-iL);
                
                const int permutationTags = iB+iC+iL;
                if (permutationTags == nTag) totalprob += bProb*cProb*lProb;
                
                
                }
            }
        }
        
    return totalprob;
    }

void pT_sort(std::pair<fastjet::PseudoJet, fastjet::PseudoJet> *dij){
    if (dij->first.pt() < dij->second.pt()){
        std::swap(dij->first, dij->second);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//                                  COMBINATORIC FUNCTIONS                                 //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////
bool overlapping_indices(vector<int>* indexvec){
// Checks if any indices in a vector are the same
    for (int o = 0; o < indexvec->size(); o++){
        for (int p = o + 1; p < indexvec->size(); p++){
            if ((*indexvec)[o] == (*indexvec)[p]){
                return true;
            }
        }
    }
    return false;
}

std::vector< std::pair<int, int> > Xchoose2(int x){
    // Returns vector of possible integer combinations of X consecutive numbers starting from 0
    std::vector< std::pair<int, int> > intvec = {};
    
    for (int i = 0; i < x; i++){
        for (int j = i + 1; j < x; j++){
            intvec.push_back(std::make_pair(i, j));
            }
        }
    return intvec;
    }

static std::vector<std::vector< std::pair<fastjet::PseudoJet, fastjet::PseudoJet> > > Jet_Combos(int higgses, 
                                std::vector<fastjet::PseudoJet> jets, std::vector<btagType> sRbt, 
                                std::vector< std::vector<btagType> >* reshuf_btags){
    // Returns a vector of possible combinations of jets which could form a group of Higgs bosons
    // i.e. a vector (complete set) OF
                // vectors (possible arrangements) OF
                            // pairs (dijets i.e. from H->bb) OF
                                    // jets (i.e. b jet)
    
    // Also fills a matching vector of vectors of btagTypes - i.e. the corresponding tag for each jet, in order
    
    typedef std::pair<int,int> intpair;
    typedef std::pair<fastjet::PseudoJet, fastjet::PseudoJet> dijet;
    typedef std::vector<dijet> hset;
    
    const int pool = jets.size();
    
    // Raw possible index combinations
    std::vector<intpair> _pchoose2 = Xchoose2(pool);
    int pchoose2 = _pchoose2.size();
    
    // Possible trios of dijets to be filled- i.e. decaying Higgs candidates
    std::vector<hset> hcand = {};
    

    // ******************* FINDING POSSIBLE COMBINATIONS ******************* // 
    
    if (higgses == 3){
        // FOR RESOLVED
        for (int i = 0; i < pchoose2; i++){
            
            for (int j = i + 1; j < pchoose2; j++){
                
                for (int k = j + 1; k < pchoose2; k++){

                    // Get indices from pair list above
                    std::vector<int> jet_indices = {_pchoose2[i].first, _pchoose2[j].first, _pchoose2[k].first, 
                                                _pchoose2[i].second, _pchoose2[j].second, _pchoose2[k].second};

                    // Check for overlapping jet indices - 2 dijets can't share same jet
                    bool overlap = overlapping_indices(&jet_indices);
              
                    if (!overlap){
                        hset temp = {dijet(jets[_pchoose2[i].first], jets[_pchoose2[i].second]),
                                    dijet(jets[_pchoose2[j].first], jets[_pchoose2[j].second]),
                                    dijet(jets[_pchoose2[k].first], jets[_pchoose2[k].second])};
                        hcand.push_back(temp);
                        
                        std::vector<btagType> btags = {sRbt[_pchoose2[i].first], sRbt[_pchoose2[i].second],
                                    sRbt[_pchoose2[j].first], sRbt[_pchoose2[j].second],
                                    sRbt[_pchoose2[k].first], sRbt[_pchoose2[k].second]};
                        reshuf_btags->push_back(btags);
                    }
                }
            }
        }
    }
    
    else if (higgses == 2){
        // FOR INTER-RESOLVED
        
        for (int i = 0; i < pchoose2; i++){
            
            for (int j = i + 1; j < pchoose2; j++){

                std::vector<int> jet_indices = {_pchoose2[i].first, _pchoose2[j].first, 
                                            _pchoose2[i].second, _pchoose2[j].second};

                // Check for overlapping jet indices - 2 dijets can't share same jet
                bool overlap = overlapping_indices(&jet_indices);
 
                if (!overlap){
                    hset temp = {dijet(jets[_pchoose2[i].first], jets[_pchoose2[i].second]),
                                    dijet(jets[_pchoose2[j].first], jets[_pchoose2[j].second])};
                    hcand.push_back(temp);
 
                    std::vector<btagType> btags = {sRbt[_pchoose2[i].first], sRbt[_pchoose2[i].second],
                                    sRbt[_pchoose2[j].first], sRbt[_pchoose2[j].second]};
                    reshuf_btags->push_back(btags);
                    }
                }
            }
        }
        
    else if (higgses == 1){
        // FOR INTER-BOOSTED
        for (int i = 0; i < pchoose2; i++){
            hset temp = {dijet(jets[_pchoose2[i].first], jets[_pchoose2[i].second])};
            hcand.push_back(temp);
            
            std::vector<btagType> btags = {sRbt[_pchoose2[i].first], sRbt[_pchoose2[i].second]};
            reshuf_btags->push_back(btags);
            }
        };
    
    return hcand;
}
    




/////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                             //
//                                      ANALYSIS FUNCTIONS                                     //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////

void triHiggsAnalysis::Analyse(bool const& signal, double const& weightnorm, finalState const& ifs){
    
    
    Analysis::Analyse(signal, weightnorm, ifs);
    
    const finalState fs = ifs;
    // *************************************** General cut selectors ************************************* //

    const Selector LR_kinematics = SelectorNHardest(3) * ( SelectorAbsRapMax(2.0) && SelectorPtMin(250.0) );
    const Selector SR_kinematics = SelectorAbsRapMax(2.5) && SelectorPtMin(20.0);
    const Selector TR_kinematics = SelectorAbsRapMax(2.5) && SelectorPtMin(50.0);
    
    // ******************************************** Clustering ********************************************//
    
    // Small-R
    const double smallR = 0.4;
    const fastjet::JetDefinition akt_sR(fastjet::antikt_algorithm, smallR);         // Algorithm to be used w param R
    const fastjet::ClusterSequence cs_akt_sR(fs, akt_sR);                          // Actual clustering
    const vector<PseudoJet> sRjets = sorted_by_pt(SR_kinematics(cs_akt_sR.inclusive_jets()));      // Sorted vec of pseudojets
    const vector<btagType> tagType_sR(BTagging(sRjets));                            // Corresponding btag info
    
    // Track Jets
    const double GAjetR = 0.3; // Boosted subjet radius for ghost-association
    const fastjet::JetDefinition jd_subjets(fastjet::antikt_algorithm, GAjetR);
    const fastjet::ClusterSequence cs_subjets(ifs, jd_subjets);
    const vector<PseudoJet> trackJets = sorted_by_pt( TR_kinematics(cs_subjets.inclusive_jets() ) );

    // Cluster large-R jets
    const double BoostJetR=1.0;
    const fastjet::JetDefinition akt_boost(fastjet::antikt_algorithm, BoostJetR);
    const fastjet::ClusterSequence cs_akt_bst(ifs, akt_boost);
    const vector<PseudoJet> largeRJets_noTrim = LR_kinematics(cs_akt_bst.inclusive_jets()); 
    const vector<PseudoJet> largeRJets_Trim = largeRJets_noTrim; //subtractPU ? trimJets(largeRJets_noTrim): largeRJets_noTrim;
    const vector<PseudoJet> largeRJets = sorted_by_pt( MDtagJets(largeRJets_Trim) );
    const vector< vector<PseudoJet> > largeRsubJets = getSubJets(largeRJets, trackJets);
    vector< vector<btagType> > tagType_LR;
    for ( auto subjets : largeRsubJets )
        tagType_LR.push_back( BTagging(subjets) );
    
    FillHistogram("CF_res", weightnorm, 0.1);
    FillHistogram("CF_interR", weightnorm, 0.1);
    FillHistogram("CF_interB", weightnorm, 0.1);
    FillHistogram("CF_boost", weightnorm, 0.1);

    ResolvedAnalysis(sRjets, tagType_sR, signal, weightnorm);
    IntermediateAnalysis_Rtype(largeRJets, sRjets, largeRsubJets, tagType_LR, tagType_sR, signal, weightnorm);
    IntermediateAnalysis_Btype(largeRJets, sRjets, largeRsubJets, tagType_LR, tagType_sR, signal, weightnorm);
    BoostedAnalysis(largeRJets, largeRsubJets, tagType_LR, signal, weightnorm);
    
    return;
}


void triHiggsAnalysis::HiggsFill(fastjet::PseudoJet const& H0, fastjet::PseudoJet const& H1, fastjet::PseudoJet const& H2, 
                std::string const& analysis, size_t const& cut, double const& weight){
                    
    const fastjet::PseudoJet triHiggs = H0 + H1 + H2;
    
    const std::string suffix = "_" + analysis + cString[cut];
    
    FillHistogram("CF_" + analysis, weight, cut + 0.1);
    
    FillHistogram("mHHH" + suffix, weight, triHiggs.m());
    FillHistogram("pT_HHH" + suffix, weight, triHiggs.pt());
    
    FillHistogram("pt_H0" + suffix, weight, H0.pt());
    FillHistogram("pt_H1" + suffix, weight, H1.pt());
    FillHistogram("pt_H2" + suffix, weight, H2.pt());

    FillHistogram("m_H0" + suffix, weight, H0.m());
    FillHistogram("m_H1" + suffix, weight, H1.m());
    FillHistogram("m_H2" + suffix, weight, H2.m());
  
    FillHistogram("eta_H0" + suffix, weight, H0.eta());
    FillHistogram("eta_H1" + suffix, weight, H1.eta());
    FillHistogram("eta_H2" + suffix, weight, H2.eta());
}


void triHiggsAnalysis::BoostFill( fastjet::PseudoJet const& H0,
                                          fastjet::PseudoJet const& H1,
                                          fastjet::PseudoJet const& H2,
                                          std::string const& analysis, 
                                          size_t const& cut, 
                                          double const& weight )
{
      // Histo fill suffix
      const std::string suffix = "_" + analysis + cString[cut];

      // Splitting scales
      const double split12_fj1 = SplittingScales( H0 );
      const double split12_fj2 = SplittingScales( H1 );
      const double split12_fj3 = SplittingScales( H2 );
      

      // 2-subjettiness / 1-subjettiness
      const double tau21_fj1 = tau21( H0 );
      const double tau21_fj2 = tau21( H1 );
      const double tau21_fj3 = tau21( H2 );

      // C2 energy correlation double-ratio
      const double C2_fj1 = C2(H0);
      const double C2_fj2 = C2(H1);
      const double C2_fj3 = C2(H2);

      // D2 energy correlation double-ratio
      const double D2_fj1 = D2(H0);
      const double D2_fj2 = D2(H1);
      const double D2_fj3 = D2(H2);

      FillHistogram("split12_fj1" + suffix, weight, split12_fj1);
      FillHistogram("split12_fj2" + suffix, weight, split12_fj2);
      FillHistogram("split12_fj3" + suffix, weight, split12_fj3);

      FillHistogram("tau21_fj1" + suffix, weight, tau21_fj1);
      FillHistogram("tau21_fj2" + suffix, weight, tau21_fj2);
      FillHistogram("tau21_fj3" + suffix, weight, tau21_fj3);

      FillHistogram("C2_fj1" + suffix, weight, C2_fj1);
      FillHistogram("C2_fj2" + suffix, weight, C2_fj2);
      FillHistogram("C2_fj3" + suffix, weight, C2_fj3);

      FillHistogram("D2_fj1" + suffix, weight, D2_fj1);
      FillHistogram("D2_fj2" + suffix, weight, D2_fj2);
      FillHistogram("D2_fj3" + suffix, weight, D2_fj3);
      
}

void triHiggsAnalysis::InterBFill( fastjet::PseudoJet const& H0,
                                          fastjet::PseudoJet const& H1,
                                          std::string const& analysis, 
                                          size_t const& cut, 
                                          double const& weight )
{
      // Histo fill suffix
      const std::string suffix = "_" + analysis + cString[cut];

      // Splitting scales
      const double split12_fj1 = SplittingScales( H0 );
      const double split12_fj2 = SplittingScales( H1 );
      

      // 2-subjettiness / 1-subjettiness
      const double tau21_fj1 = tau21( H0 );
      const double tau21_fj2 = tau21( H1 );

      // C2 energy correlation double-ratio
      const double C2_fj1 = C2(H0);
      const double C2_fj2 = C2(H1);

      // D2 energy correlation double-ratio
      const double D2_fj1 = D2(H0);
      const double D2_fj2 = D2(H1);

      FillHistogram("split12_fj1" + suffix, weight, split12_fj1);
      FillHistogram("split12_fj2" + suffix, weight, split12_fj2);

      FillHistogram("tau21_fj1" + suffix, weight, tau21_fj1);
      FillHistogram("tau21_fj2" + suffix, weight, tau21_fj2);

      FillHistogram("C2_fj1" + suffix, weight, C2_fj1);
      FillHistogram("C2_fj2" + suffix, weight, C2_fj2);

      FillHistogram("D2_fj1" + suffix, weight, D2_fj1);
      FillHistogram("D2_fj2" + suffix, weight, D2_fj2);
      
}

void triHiggsAnalysis::InterRFill( fastjet::PseudoJet const& H0,
                                          std::string const& analysis, 
                                          size_t const& cut, 
                                          double const& weight )
{
      // Histo fill suffix
      const std::string suffix = "_" + analysis + cString[cut];

      // Splitting scales
      const double split12_fj1 = SplittingScales( H0 );
      

      // 2-subjettiness / 1-subjettiness
      const double tau21_fj1 = tau21( H0 );

      // C2 energy correlation double-ratio
      const double C2_fj1 = C2(H0);

      // D2 energy correlation double-ratio
      const double D2_fj1 = D2(H0);

      FillHistogram("split12_fj" + suffix, weight, split12_fj1);

      FillHistogram("tau21_fj1" + suffix, weight, tau21_fj1);

      FillHistogram("C2_fj1" + suffix, weight, C2_fj1);

      FillHistogram("D2_fj1" + suffix, weight, D2_fj1);
}

double triHiggsAnalysis::ResolvedAnalysis(const vector<PseudoJet>& jets, const vector<btagType>& btags, const bool& signal, const double& weight){
    
    if (jets.size() >= 6){
        
        //*************** Bookkeeping *****************//
        FillHistogram("b0_pT_res", weight, jets[0].pt());
        FillHistogram("b1_pT_res", weight, jets[1].pt());
        FillHistogram("b2_pT_res", weight, jets[2].pt());
        FillHistogram("b3_pT_res", weight, jets[3].pt());
        FillHistogram("b4_pT_res", weight, jets[4].pt());
        FillHistogram("b5_pT_res", weight, jets[5].pt());
        for (int pt = 0; pt < 100; pt+=10){
            if (jets[5].pt() >= pt){
                FillHistogram("pT_cuts_6j_res", weight, pt+5);
                }
            }
        //*********************************************//
        
        typedef std::pair<fastjet::PseudoJet, fastjet::PseudoJet> dijet;
        typedef std::vector<dijet> htrio;

        
        vector< vector<btagType> > _btags_reshuf;
        vector< vector<btagType> >* btags_reshuf = &_btags_reshuf;
        vector< vector<dijet> > hcand = Jet_Combos(3, jets, btags, btags_reshuf);

        
        // ********************************** Choosing the best candidate **********************************//

        double min_hcand_variance = INFINITY;
        double min_hcand_variance_index;
        for (int i = 0; i < hcand.size(); i++){
            double mass1 = (hcand[i][0].first + hcand[i][0].second).m();
            double mass2 = (hcand[i][1].first + hcand[i][1].second).m();
            double mass3 = (hcand[i][2].first + hcand[i][2].second).m();
            double variance = (mass1*mass1 + mass2*mass2 + mass3*mass3)/3 - pow(((mass1 + mass2 + mass3)/3),2);
            if (variance < min_hcand_variance){
                min_hcand_variance = variance;
                min_hcand_variance_index = i;
                }
        }
        htrio triHiggs_cand = hcand[min_hcand_variance_index];

        dijet* h1 = &triHiggs_cand[0]; dijet* h2 = &triHiggs_cand[1]; dijet* h3 = &triHiggs_cand[2]; 
 
        vector<fastjet::PseudoJet> sorted_higgses = {h1->first + h1->second, h2->first + h2->second, h3->first + h3->second};
        
        // Sort by pT
        sorted_higgses = sorted_by_pt(sorted_higgses);
        pT_sort(h1);
        pT_sort(h2);
        pT_sort(h3);
        
        // ***************************************** B-Tagging ******************************************//
        
        // Probability that these 6 will all be b-tagged
        const int nB = std::count(btags.begin(), btags.begin() + 6, BTAG); 
        const int nC = std::count(btags.begin(), btags.begin() + 6, CTAG);
        const int nL = std::count(btags.begin(), btags.begin() + 6, LTAG);
        const double sel_eff = btagprob(6, nB, nC, nL);
        const double sel_wgt = sel_eff * weight;
        
        for (int pt = 0; pt < 100; pt+=10){
            if (jets[5].pt() >= pt){
                FillHistogram("pT_cuts_6b_res", sel_wgt, pt+5);
                }
        }
        if (sel_wgt == 0) return 0;
        
        HiggsFill(sorted_higgses[0], sorted_higgses[1], sorted_higgses[2], "res", 1, sel_wgt);          // (c1)

        // *************************************** Mass window *****************************************//
        const double massdiff1 = fabs(sorted_higgses[0].m() - 125.0);
        const double massdiff2 = fabs(sorted_higgses[1].m() - 125.0);
        const double massdiff3 = fabs(sorted_higgses[2].m() - 125.0);
        const double masswindow = 40.0;
        
        if (massdiff1 <= masswindow && massdiff2 <= masswindow && massdiff3 <= masswindow){
            HiggsFill(sorted_higgses[0], sorted_higgses[1], sorted_higgses[2], "res", 2, sel_wgt);      // (c2)
            
            fastjet::PseudoJet triHiggs = sorted_higgses[0] + sorted_higgses[1] + sorted_higgses[2];
            
            // Bookkeeping
            FillHistogram("deltaR_02_res", sel_wgt, sorted_higgses[0].delta_R(sorted_higgses[2]));
            FillHistogram("deltaR_12_res", sel_wgt, sorted_higgses[1].delta_R(sorted_higgses[2]));
            
            
            resNTuple << signal <<"\t"<<GetSample()<<"\t"<<sel_wgt << "\t"
          << sorted_higgses[0].pt() << "\t"
          << sorted_higgses[1].pt() << "\t"
          << sorted_higgses[2].pt() << "\t"
          << triHiggs.pt() << "\t"
          << sorted_higgses[0].m() << "\t"
          << sorted_higgses[1].m() << "\t"
          << sorted_higgses[2].m() << "\t"
          << triHiggs.m() << "\t"
          << sorted_higgses[0].delta_R(sorted_higgses[1]) << "\t"
          << sorted_higgses[1].delta_R(sorted_higgses[2]) << "\t"
          << sorted_higgses[0].delta_R(sorted_higgses[2]) << "\t"
          << getDPhi(sorted_higgses[0].phi(), sorted_higgses[1].phi()) << "\t"
          << getDPhi(sorted_higgses[0].phi(), sorted_higgses[2].phi()) << "\t"
          << getDPhi(sorted_higgses[1].phi(), sorted_higgses[2].phi()) << "\t"
          << fabs( sorted_higgses[0].eta() - sorted_higgses[1].eta())  << "\t"
          << fabs( sorted_higgses[1].eta() - sorted_higgses[2].eta())  << "\t"
          << fabs( sorted_higgses[0].eta() - sorted_higgses[2].eta())  << "\t"
          << Chi( sorted_higgses[0], sorted_higgses[1])  << "\t"
          << Chi( sorted_higgses[2], sorted_higgses[1])  << "\t"
          << Chi( sorted_higgses[0], sorted_higgses[2])  << "\t"
          << h1->first.pt() << "\t"
          << h1->second.pt() << "\t"
          << h2->first.pt() << "\t"
          << h2->second.pt() << "\t"
          << h3->first.pt() << "\t"
          << h3->second.pt() << "\t"
          << std::endl;
        }
        else HiggsFill(sorted_higgses[0], sorted_higgses[1], sorted_higgses[2], "res", 3, sel_wgt);     // (c3)
        return sel_wgt;
        
        }
    return 0;
}


double triHiggsAnalysis::IntermediateAnalysis_Rtype( const std::vector<PseudoJet>& largeRJets,
                                                     const std::vector<PseudoJet>& smallRJets,
                                                     const std::vector< std::vector<PseudoJet> >& largeRsubJets,
                                                     const std::vector< std::vector<btagType> > largeRbtags,
                                                     const std::vector<btagType>& smallRbtags,
                                                     const bool& signal, const double& event_weight )
{
    if (largeRJets.size() == 1 && largeRbtags[0].size() >= 2){
            
            typedef std::pair<fastjet::PseudoJet, fastjet::PseudoJet> dijet;
            
            //Bookkeeping
            FillHistogram("lR_higgs", event_weight, largeRJets[0].m());
            
            // Pushing smallR jets sufficiently far from largeR Higgs
            std::vector<fastjet::PseudoJet> sRjets;
            double R_range = 1.1;
            
            for (int i = 0; i < smallRJets.size(); i++ ){
                if (largeRJets[0].delta_R(smallRJets[i]) > R_range){
                    sRjets.push_back(smallRJets[i]);
                }
            }
            if (sRjets.size() < 4) return 0;
            
            // Possible combinations of the remaining sR jets into dijets (Higgses)
            std::vector< std::vector<btagType> > _sRbtags_reshuf;
            std::vector< std::vector<btagType> >* sRbtags_reshuf = &_sRbtags_reshuf;
            std::vector< std::vector<dijet> > sR_possibles = Jet_Combos(2, sRjets, smallRbtags, sRbtags_reshuf);
     

            // ********************************** Choosing the best candidate **********************************//
            double min_hcand_variance = INFINITY;
            double min_hcand_variance_index = 0;
            for (int i = 0; i < sR_possibles.size(); i++){
                double mass1 = largeRJets[0].m();
                double mass2 = (sR_possibles[i][0].first + sR_possibles[i][0].second).m();
                double mass3 = (sR_possibles[i][1].first + sR_possibles[i][1].second).m();
                double variance = (mass1*mass1 + mass2*mass2 + mass3*mass3)/3 - pow(((mass1 + mass2 + mass3)/3),2);
                if (variance < min_hcand_variance){
                    min_hcand_variance = variance;
                    min_hcand_variance_index = i;
                    }
            }
            std::vector<dijet> sR_best = sR_possibles[min_hcand_variance_index];
      
            std::vector<btagType> sR_best_btags = (*sRbtags_reshuf)[min_hcand_variance_index];
            
            
            dijet* h2 = &sR_best[0];
            dijet* h3 = &sR_best[1];
            vector<fastjet::PseudoJet> sorted_higgses = {largeRJets[0], h2->first + h2->second, h3->first + h3->second};
            
            // Sort by pT
            sorted_higgses = sorted_by_pt(sorted_higgses);
            pT_sort(h2);
            pT_sort(h3);

            // ****************************************** B-Tagging *******************************************//
            
            // Small-R jets
            int nB = std::count(sR_best_btags.begin(), sR_best_btags.begin() + 4, BTAG);
            int nC = std::count(sR_best_btags.begin(), sR_best_btags.begin() + 4, CTAG);
            int nL = std::count(sR_best_btags.begin(), sR_best_btags.begin() + 4, LTAG);
            
            // Large-R jet - test leading and subleading subjets
            for (int i = 0; i < 2; i++){
                if (largeRbtags[0][i] == BTAG) nB++;
                else if(largeRbtags[0][i] == CTAG) nC++;
                else nL++;
            }
            
            const double sel_eff = btagprob(6, nB, nC, nL);
            const double sel_wgt = sel_eff * event_weight;
            if (sel_wgt == 0) return 0;
            
            HiggsFill(sorted_higgses[0], sorted_higgses[1], sorted_higgses[2], "interR", 1, sel_wgt);       // (c1)
            InterRFill(largeRJets[0], "interR", 1, sel_wgt);
            
            // *************************************** Mass window *****************************************//
            const double massdiff1 = fabs(sorted_higgses[0].m() - 125.0);
            const double massdiff2 = fabs(sorted_higgses[1].m() - 125.0);
            const double massdiff3 = fabs(sorted_higgses[2].m() - 125.0);
            const double masswindow = 40.0;
      
            if (massdiff1 <= masswindow && massdiff2 <= masswindow && massdiff3 <= masswindow){
                HiggsFill(sorted_higgses[0], sorted_higgses[1], sorted_higgses[2], "interR", 2, sel_wgt);   // (c2)
                InterRFill(largeRJets[0], "interR", 2, sel_wgt);
                
                //Bookkeeping
     //           FillHistogram("deltaR_02_interR", sel_wgt, sorted_higgses[0].delta_R(sorted_higgses[2]));
      //          FillHistogram("deltaR_12_interR", sel_wgt, sorted_higgses[1].delta_R(sorted_higgses[2]));
       //         FillHistogram("lR_higgs_MW", event_weight, largeRJets[0].m());
                
                
                fastjet::PseudoJet triHiggs = sorted_higgses[0] + sorted_higgses[1] + sorted_higgses[2];
                
                // Calculate some substructure variables
                const std::vector<double> split12_vec = SplittingScales( largeRJets );
                const double tau21_fj1 = tau21( largeRJets[0] );

                // C2 energy correlation double-ratio
                const double C2_fj1 = C2(largeRJets[0]);

                // D2 energy correlation double-ratio
                const double D2_fj1 = D2(largeRJets[0]);
                
                
                intRNTuple << signal <<"\t"<<GetSample()<<"\t"<<sel_wgt << "\t"
                << sorted_higgses[0].pt() << "\t"
                << sorted_higgses[1].pt() << "\t"
                << sorted_higgses[2].pt() << "\t"
                << triHiggs.pt() << "\t"
                << sorted_higgses[0].m() << "\t"
                << sorted_higgses[1].m() << "\t"
                << sorted_higgses[2].m() << "\t"
                << triHiggs.m() << "\t"
                << sorted_higgses[0].delta_R(sorted_higgses[1]) << "\t"
                << sorted_higgses[1].delta_R(sorted_higgses[2]) << "\t"
                << sorted_higgses[0].delta_R(sorted_higgses[2]) << "\t"
                << getDPhi(sorted_higgses[0].phi(), sorted_higgses[1].phi()) << "\t"
                << getDPhi(sorted_higgses[0].phi(), sorted_higgses[2].phi()) << "\t"
                << getDPhi(sorted_higgses[1].phi(), sorted_higgses[2].phi()) << "\t"
                << fabs( sorted_higgses[0].eta() - sorted_higgses[1].eta())  << "\t"
                << fabs( sorted_higgses[1].eta() - sorted_higgses[2].eta())  << "\t"
                << fabs( sorted_higgses[0].eta() - sorted_higgses[2].eta())  << "\t"
                << Chi( sorted_higgses[0], sorted_higgses[1])  << "\t"
                << Chi( sorted_higgses[2], sorted_higgses[1])  << "\t"
                << Chi( sorted_higgses[0], sorted_higgses[2])  << "\t"
                << largeRsubJets[0][0].pt() << "\t"
                << largeRsubJets[0][1].pt() << "\t"
                << h2->first.pt() << "\t"
                << h2->second.pt() << "\t"
                << h3->first.pt() << "\t"
                << h3->second.pt() << "\t"
                << split12_vec[0] << "\t"
                << tau21_fj1 << "\t"
                << C2_fj1 << "\t"
                << D2_fj1 << "\t"
              << std::endl;
                }
            else HiggsFill(sorted_higgses[0], sorted_higgses[1], sorted_higgses[2], "interR", 3, sel_wgt);  // (c3)
            return sel_wgt;
        
    }
    return 0;
    
}

double triHiggsAnalysis::IntermediateAnalysis_Btype( const std::vector<PseudoJet>& largeRJets,
                                                     const std::vector<PseudoJet>& smallRJets,
                                                     const std::vector< std::vector<PseudoJet> >& largeRsubJets,
                                                     const std::vector< std::vector<btagType> > largeRbtags,
                                                     const std::vector<btagType>& smallRbtags,
                                                     const bool& signal, const double& event_weight )
{
//    std::cout << "IB begin" << std::endl;
    if (largeRJets.size() == 2){
        
        if (largeRbtags.size() < 2 || largeRbtags[0].size() < 2 || largeRbtags[1].size() < 2) return 0;
        
        
        typedef std::pair<fastjet::PseudoJet, fastjet::PseudoJet> dijet;
        
        std::vector<fastjet::PseudoJet> sRjets;
        
        // Pushing smallR jets sufficiently far from largeR Higgses
        double R_range = 1.1;
        for (int i = 0; i < smallRJets.size(); i++ ){
            if (largeRJets[0].delta_R(smallRJets[i]) > R_range && largeRJets[1].delta_R(smallRJets[i]) > R_range){
                sRjets.push_back(smallRJets[i]);
            }
        }
        if (sRjets.size() < 2) return 0;
        
        // Possible combinations of the remaining sR jets into dijets (Higgses)
        std::vector< std::vector<btagType> > _sRbtags_reshuf;
        std::vector< std::vector<btagType> >* sRbtags_reshuf = &_sRbtags_reshuf;
        std::vector< std::vector<dijet> > sR_possibles = Jet_Combos(1, sRjets, smallRbtags, sRbtags_reshuf);
        
        
        // ********************************** Choosing the best candidate **********************************//
        double min_hcand_variance = INFINITY;
        double min_hcand_variance_index;
        for (int i = 0; i < sR_possibles.size(); i++){
            double mass1 = largeRJets[0].m();
            double mass2 = largeRJets[1].m();
            double mass3 = (sR_possibles[i][0].first + sR_possibles[i][0].second).m();
            double variance = (mass1*mass1 + mass2*mass2 + mass3*mass3)/3 - pow(((mass1 + mass2 + mass3)/3),2);
            if (variance < min_hcand_variance){
                min_hcand_variance = variance;
                min_hcand_variance_index = i;
                }
        }
        vector<dijet> sR_best = sR_possibles[min_hcand_variance_index];
    
        vector<btagType> sR_best_btags = (*sRbtags_reshuf)[min_hcand_variance_index];
        
        dijet* h3 = &sR_best[0];
        
        vector<fastjet::PseudoJet> sorted_higgses = {largeRJets[0], largeRJets[1], h3->first + h3->second};
    
        // Sort by pT
        sorted_higgses = sorted_by_pt(sorted_higgses);
        pT_sort(h3);

        // ****************************************** B-Tagging *******************************************//
        // Small-R jets
        int nB = std::count(sR_best_btags.begin(), sR_best_btags.begin() + 2, BTAG);
        int nC = std::count(sR_best_btags.begin(), sR_best_btags.begin() + 2, CTAG);
        int nL = std::count(sR_best_btags.begin(), sR_best_btags.begin() + 2, LTAG);
        
        // Large-R jets - test leading and subleading subjets
        for (int j = 0; j < 2; j++){
            for (int i = 0; i < 2; i++){
                if (largeRbtags[j][i] == BTAG) nB++;
                else if(largeRbtags[j][i] == CTAG) nC++;
                else nL++;
            }
        }
        
        const double sel_eff = btagprob(6, nB, nC, nL);       // Prob that all will be b-tagged
        const double sel_wgt = sel_eff * event_weight;
        if (sel_wgt == 0) return 0;
        
        HiggsFill(sorted_higgses[0], sorted_higgses[1], sorted_higgses[2], "interB", 1, sel_wgt);       // (c1)
        InterBFill(largeRJets[0], largeRJets[1], "interB", 1, sel_wgt);
        
        
        // *************************************** Mass window *****************************************//
   
        const double massdiff1 = fabs(sorted_higgses[0].m() - 125.0);
        const double massdiff2 = fabs(sorted_higgses[1].m() - 125.0);
        const double massdiff3 = fabs(sorted_higgses[2].m() - 125.0);
        const double masswindow = 40.0;
        
        if (massdiff1 <= masswindow && massdiff2 <= masswindow && massdiff3 <= masswindow){
        
            HiggsFill(sorted_higgses[0], sorted_higgses[1], sorted_higgses[2], "interB", 2, sel_wgt);   // (c2)
     
            InterBFill(largeRJets[0], largeRJets[1], "interB", 2, sel_wgt);
            // Bookkeeping
 //           FillHistogram("deltaR_02_interB", sel_wgt, sorted_higgses[0].delta_R(sorted_higgses[2]));
  //          FillHistogram("deltaR_12_interB", sel_wgt, sorted_higgses[1].delta_R(sorted_higgses[2]));
   //         FillHistogram("lR_higgs_MW", event_weight, largeRJets[0].m());
            
            
            fastjet::PseudoJet triHiggs = sorted_higgses[0] + sorted_higgses[1] + sorted_higgses[2];
            
            // Calculate some substructure variables
            const std::vector<double> split12_vec = SplittingScales( largeRJets );
            const double tau21_fj1 = tau21( largeRJets[0] );
            const double tau21_fj2 = tau21( largeRJets[1] );

            // C2 energy correlation double-ratio
            const double C2_fj1 = C2(largeRJets[0]);
            const double C2_fj2 = C2(largeRJets[1]);

            // D2 energy correlation double-ratio
            const double D2_fj1 = D2(largeRJets[0]);
            const double D2_fj2 = D2(largeRJets[1]);
            
            intBNTuple << signal <<"\t"<<GetSample()<<"\t"<<sel_wgt << "\t"
             << sorted_higgses[0].pt() << "\t"
            << sorted_higgses[1].pt() << "\t"
            << sorted_higgses[2].pt() << "\t"
            << triHiggs.pt() << "\t"
            << sorted_higgses[0].m() << "\t"
            << sorted_higgses[1].m() << "\t"
            << sorted_higgses[2].m() << "\t"
            << triHiggs.m() << "\t"
            << sorted_higgses[0].delta_R(sorted_higgses[1]) << "\t"
            << sorted_higgses[1].delta_R(sorted_higgses[2]) << "\t"
            << sorted_higgses[0].delta_R(sorted_higgses[2]) << "\t"
            << getDPhi(sorted_higgses[0].phi(), sorted_higgses[1].phi()) << "\t"
            << getDPhi(sorted_higgses[0].phi(), sorted_higgses[2].phi()) << "\t"
            << getDPhi(sorted_higgses[1].phi(), sorted_higgses[2].phi()) << "\t"
            << fabs( sorted_higgses[0].eta() - sorted_higgses[1].eta())  << "\t"
            << fabs( sorted_higgses[1].eta() - sorted_higgses[2].eta())  << "\t"
            << fabs( sorted_higgses[0].eta() - sorted_higgses[2].eta())  << "\t"
            << Chi( sorted_higgses[0], sorted_higgses[1])  << "\t"
            << Chi( sorted_higgses[1], sorted_higgses[2])  << "\t"
            << Chi( sorted_higgses[0], sorted_higgses[2])  << "\t"
            << largeRsubJets[0][0].pt() << "\t"
            << largeRsubJets[0][1].pt() << "\t"
            << largeRsubJets[1][0].pt() << "\t"
            << largeRsubJets[1][1].pt() << "\t"
            << h3->first.pt() << "\t"
            << h3->second.pt() << "\t"
            << split12_vec[0] << "\t"
            << split12_vec[1] << "\t"
            << tau21_fj1 << "\t"
            << tau21_fj2 << "\t"
            << C2_fj1 << "\t"
            << C2_fj2 << "\t"
            << D2_fj1 << "\t"
            << D2_fj2 << "\t"
            << std::endl;
            std::cout << "done" << std::endl;
            }
        else HiggsFill(sorted_higgses[0], sorted_higgses[1], sorted_higgses[2], "interB", 3, sel_wgt);  // (c3)
        
        return sel_wgt;
        
    }
    return 0;
}

double triHiggsAnalysis::BoostedAnalysis( const std::vector<PseudoJet>& largeRJets,
                                                     const std::vector< std::vector<PseudoJet> >& largeRsubJets,
                                                     const std::vector< std::vector<btagType> > largeRbtags,
                                                     const bool& signal, const double& event_weight )
{
    if (largeRJets.size() == 3){
        
        std::vector<fastjet::PseudoJet> sorted_higgses = sorted_by_pt(largeRJets);
        
        // ****************************************** B-Tagging *******************************************//
        // Test leading and subleading subjets for each jet
        int nB = 0; int nC = 0; int nL = 0;
        for (int j = 0; j < 3; j++){
            for (int i = 0; i < 2; i++){
                if (largeRbtags[j][i] == BTAG) nB++;
                else if(largeRbtags[j][i] == CTAG) nC++;
                else nL++;
            }
        }
        
        const double sel_eff = btagprob(6, nB, nC, nL);       // Prob that all will be b-tagged
        const double sel_wgt = sel_eff * event_weight;
        if (sel_wgt == 0) return 0;
        
        
        HiggsFill(sorted_higgses[0], sorted_higgses[1], sorted_higgses[2], "boost", 1, sel_wgt);       // (c1)
        
        
        // *************************************** Mass window *****************************************//
   
        const double massdiff1 = fabs(sorted_higgses[0].m() - 125.0);
        const double massdiff2 = fabs(sorted_higgses[1].m() - 125.0);
        const double massdiff3 = fabs(sorted_higgses[2].m() - 125.0);
        const double masswindow = 40.0;
  
        if (massdiff1 <= masswindow && massdiff2 <= masswindow && massdiff3 <= masswindow){
            HiggsFill(sorted_higgses[0], sorted_higgses[1], sorted_higgses[2], "boost", 2, sel_wgt);   // (c2)
            
            
            fastjet::PseudoJet triHiggs = sorted_higgses[0] + sorted_higgses[1] + sorted_higgses[2];
            
            // Calculate some substructure variables
            const std::vector<double> split12_vec = SplittingScales( largeRJets );
            const double tau21_fj1 = tau21( largeRJets[0] );
            const double tau21_fj2 = tau21( largeRJets[1] );
            const double tau21_fj3 = tau21( largeRJets[2] );

            // C2 energy correlation double-ratio
            const double C2_fj1 = C2(largeRJets[0]);
            const double C2_fj2 = C2(largeRJets[1]);
            const double C2_fj3 = C2(largeRJets[2]);

            // D2 energy correlation double-ratio
            const double D2_fj1 = D2(largeRJets[0]);
            const double D2_fj2 = D2(largeRJets[1]);
            const double D2_fj3 = D2(largeRJets[2]);
            
            bstNTuple << signal <<"\t"<<GetSample()<<"\t"<<sel_wgt << "\t"
             << sorted_higgses[0].pt() << "\t"
            << sorted_higgses[1].pt() << "\t"
            << sorted_higgses[2].pt() << "\t"
            << triHiggs.pt() << "\t"
            << sorted_higgses[0].m() << "\t"
            << sorted_higgses[1].m() << "\t"
            << sorted_higgses[2].m() << "\t"
            << triHiggs.m() << "\t"
            << sorted_higgses[0].delta_R(sorted_higgses[1]) << "\t"
            << sorted_higgses[1].delta_R(sorted_higgses[2]) << "\t"
            << sorted_higgses[0].delta_R(sorted_higgses[2]) << "\t"
            << getDPhi(sorted_higgses[0].phi(), sorted_higgses[1].phi()) << "\t"
            << getDPhi(sorted_higgses[0].phi(), sorted_higgses[2].phi()) << "\t"
            << getDPhi(sorted_higgses[1].phi(), sorted_higgses[2].phi()) << "\t"
            << fabs( sorted_higgses[0].eta() - sorted_higgses[1].eta())  << "\t"
            << fabs( sorted_higgses[1].eta() - sorted_higgses[2].eta())  << "\t"
            << fabs( sorted_higgses[0].eta() - sorted_higgses[2].eta())  << "\t"
            << Chi( sorted_higgses[0], sorted_higgses[1])  << "\t"
            << Chi( sorted_higgses[1], sorted_higgses[2])  << "\t"
            << Chi( sorted_higgses[0], sorted_higgses[2])  << "\t"
            << largeRsubJets[0][0].pt() << "\t"
            << largeRsubJets[0][1].pt() << "\t"
            << largeRsubJets[1][0].pt() << "\t"
            << largeRsubJets[1][1].pt() << "\t"
            << largeRsubJets[2][0].pt() << "\t"
            << largeRsubJets[2][1].pt() << "\t"
            << split12_vec[0] << "\t"
            << split12_vec[1] << "\t"
            << split12_vec[2] << "\t"
            << tau21_fj1 << "\t"
            << tau21_fj2 << "\t"
            << tau21_fj3 << "\t"
            << C2_fj1 << "\t"
            << C2_fj2 << "\t"
            << C2_fj3 << "\t"
            << D2_fj1 << "\t"
            << D2_fj2 << "\t"
            << D2_fj3 << "\t"
            <<std::endl;
            }
        else HiggsFill(sorted_higgses[0], sorted_higgses[1], sorted_higgses[2], "boost", 3, sel_wgt);  // (c3)
        
        return sel_wgt;
        
    }
    
    return 0;
}
