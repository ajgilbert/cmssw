// -*- C++ -*-
//
// Package:    ConversionTrackRefUpdater
// Class:      ConversionTrackRefUpdater
//
/**\class ConversionTrackRefUpdater ConversionTrackRefUpdater.cc TauAnalysis/ConversionTrackRefUpdater/src/ConversionTrackRefUpdater.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Tomasz Maciej Frueboes
//         Created:  Fri Apr  9 12:15:56 CEST 2010
// $Id: ConversionTrackRefUpdater.cc,v 1.1 2012/03/01 17:03:28 fruboes Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/EgammaTrackReco/interface/ConversionTrack.h"
#include "DataFormats/EgammaTrackReco/interface/ConversionTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"

#include <DataFormats/Math/interface/deltaR.h>
//
// class decleration
//
using namespace reco;


class ConversionTrackRefUpdater : public edm::EDProducer {
  public:
    explicit ConversionTrackRefUpdater(const edm::ParameterSet&);
    ~ConversionTrackRefUpdater();

  private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    edm::InputTag input_;
    edm::InputTag srcTracks_;
    edm::InputTag targetTracks_;
};

//
// constructors and destructor
//
ConversionTrackRefUpdater::ConversionTrackRefUpdater(const edm::ParameterSet& iConfig) :
    input_(iConfig.getParameter<edm::InputTag>("input")),
    srcTracks_(iConfig.getParameter<edm::InputTag>("srcTracks")),
    targetTracks_(iConfig.getParameter<edm::InputTag>("targetTracks")) {
  produces<ConversionTrackCollection>();
}


ConversionTrackRefUpdater::~ConversionTrackRefUpdater() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ConversionTrackRefUpdater::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   auto_ptr<ConversionTrackCollection> output(new ConversionTrackCollection);

   edm::Handle<reco::TrackCollection> hTargetTracks;
   iEvent.getByLabel(targetTracks_, hTargetTracks);
   
   edm::Handle<reco::TrackCollection> hSrcTracks;
   iEvent.getByLabel(srcTracks_, hSrcTracks);

   edm::Handle<ConversionTrackCollection> hConversions;
   iEvent.getByLabel(input_, hConversions);

   for (unsigned i = 0; i < hConversions->size(); ++i) {
     ConversionTrack newConvTrack(hConversions->at(i));
     TrackBaseRef currentRef = hConversions->at(i).trackRef();
     if (currentRef.isNull()) {
       output->push_back(newConvTrack);
       continue;
     }
     if (currentRef.id() == hSrcTracks.id()) {
       size_t newIndex = -1;
       bool found = false;
       for (size_t i = 0; i < hTargetTracks->size(); ++i) {
         if (deltaR(currentRef->momentum(), hTargetTracks->at(i).momentum()) < 0.001){
           newIndex = i;
           found = true;
           break;
         }
       }
       if (found) {
         TrackBaseRef trackRef(TrackRef(hTargetTracks, newIndex));
         ConversionTrack tmpConvTrack(trackRef);
         tmpConvTrack.setTrajRef(newConvTrack.trajRef());
         tmpConvTrack.setIsTrackerOnly(newConvTrack.isTrackerOnly());
         tmpConvTrack.setIsArbitratedEcalSeeded(newConvTrack.isArbitratedEcalSeeded());
         tmpConvTrack.setIsArbitratedMerged(newConvTrack.isArbitratedMerged());
         tmpConvTrack.setIsArbitratedMergedEcalGeneral(newConvTrack.isArbitratedMergedEcalGeneral());
         output->push_back(tmpConvTrack);
       } else {
         std::cout << "XXXX whoops! Cannot set track ref!!" << std::endl;
       }
     } else {
      output->push_back(newConvTrack);
     }
   }
   iEvent.put(output);
}

// ------------ method called once each job just before starting event loop  ------------
void
ConversionTrackRefUpdater::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ConversionTrackRefUpdater::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ConversionTrackRefUpdater);
