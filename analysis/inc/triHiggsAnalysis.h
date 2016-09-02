#pragma once
#include "analysis.h"
#include "utils.h"

#include <iostream>
#include <vector>


class triHiggsAnalysis : public Analysis
{
public:
    triHiggsAnalysis(runCard const& run, sampleCard const& sample, int const& subsample);
    
    void Analyse(bool const& signal, double const& weight_norm, finalState const&);
    
    

private:

    double ResolvedAnalysis(const std::vector<fastjet::PseudoJet>& smallRJets, 
                            const std::vector<btagType>& btags, 
                            const bool& signal, const double& event_weight);
    
    double IntermediateAnalysis_Rtype( const std::vector<fastjet::PseudoJet>& largeRJets,
                                                     const std::vector<fastjet::PseudoJet>& smallRJets,
                                                     const std::vector< std::vector<fastjet::PseudoJet> >& largeRsubJets,
                                                     const std::vector< std::vector<btagType> > largeRbtags,
                                                     const std::vector<btagType>& smallRbtags,
                                                     const bool& signal, const double& event_weight );
    double IntermediateAnalysis_Btype( const std::vector<fastjet::PseudoJet>& largeRJets,
                                                     const std::vector<fastjet::PseudoJet>& smallRJets,
                                                     const std::vector< std::vector<fastjet::PseudoJet> >& largeRsubJets,
                                                     const std::vector< std::vector<btagType> > largeRbtags,
                                                     const std::vector<btagType>& smallRbtags,
                                                     const bool& signal, const double& event_weight );
    double BoostedAnalysis( const std::vector<fastjet::PseudoJet>& largeRJets,
                                                     const std::vector< std::vector<fastjet::PseudoJet> >& largeRsubJets,
                                                     const std::vector< std::vector<btagType> > largeRbtags,
                                                     const bool& signal, const double& event_weight );
                                                     
    void HiggsFill(fastjet::PseudoJet const& H0, fastjet::PseudoJet const& H1, fastjet::PseudoJet const& H2, 
                std::string const& analysis, size_t const& cut, double const& weight);
                
    void BoostFill( fastjet::PseudoJet const& H0,
                                          fastjet::PseudoJet const& H1,
                                          fastjet::PseudoJet const& H2,
                                          std::string const& analysis, 
                                          size_t const& cut, 
                                          double const& weight );
    void InterBFill( fastjet::PseudoJet const& H0,
                                          fastjet::PseudoJet const& H1,
                                          std::string const& analysis, 
                                          size_t const& cut, 
                                          double const& weight );
    void InterRFill( fastjet::PseudoJet const& H0,
                                          std::string const& analysis, 
                                          size_t const& cut, 
                                          double const& weight );
                
    std::ofstream resNTuple;
    std::ofstream intRNTuple;
    std::ofstream intBNTuple;
    std::ofstream bstNTuple;
};

