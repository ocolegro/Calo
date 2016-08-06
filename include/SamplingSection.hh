#ifndef _samplingsection_hh_
#define _samplingsection_hh_

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Colour.hh"
#include "G4Track.hh"
//#include "EventAction.hh"
#include "G4RunManager.hh"

#include <iomanip>
#include <vector>
#include "G4SiHit.hh"
//class EventAction;
class SamplingSection {
public:
	//CTOR
	SamplingSection(std::vector<std::pair <G4double,std::string>> iEle);


			//DTOR
			~SamplingSection() {};

			inline void setNumberOfSectors(const unsigned nSec) {
				n_sectors = nSec;
				sublayer_vol.clear();
				for (unsigned ie(0); ie<n_elements*n_sectors; ++ie) {
					sublayer_vol.push_back(0);
				}
			};

			//
			std::pair<G4bool,G4bool> add(G4double depositRawE,G4VPhysicalVolume* vol,G4Track *lTrack,const G4ThreeVector & position);

			inline bool isSensitiveElement(const unsigned & aEle) {
				if (aEle < n_elements &&
						(ele_name[aEle] == "Si" || ele_name[aEle] == "Scintillator")
				) return true;
				return false;
			};

			inline unsigned getSensitiveLayerIndex(std::string astr) {
				if (astr.find("_")== astr.npos) return 0;
				size_t pos = astr.find("phys");
				if (pos != astr.npos && pos>1) {
					unsigned idx = 0;
					std::istringstream(astr.substr(pos-1,1))>>idx;
					return idx;
				}
				return 0;
			};

			inline G4Colour g4Colour(const unsigned & aEle) {
				if (isSensitiveElement(aEle)) return G4Colour::Red();
				if (ele_name[aEle] == "Cu") return G4Colour::Black();
				if (isAbsorberElement(aEle)) return G4Colour::Gray();
				if (ele_name[aEle] == "PCB") return G4Colour::Blue();
				if (ele_name[aEle] == "Air") return G4Colour::Cyan();
				return G4Colour::Yellow();
			};

			inline bool isAbsorberElement(const unsigned & aEle) {
				if (aEle < n_elements &&
						(
								ele_name[aEle] == "Pb" || ele_name[aEle] == "Cu" ||
								ele_name[aEle] == "W" || ele_name[aEle] == "Brass" ||
								ele_name[aEle] == "Fe" || ele_name[aEle] == "Steel" ||
								ele_name[aEle] == "SSteel" || ele_name[aEle] == "Al" ||
								ele_name[aEle] == "WCu" || ele_name[aEle] == "NeutMod"
						)
				) return true;
				return false;
			};

			//reset
			inline void resetCounters()
			{
				sublayer_RawDep.clear();
				sublayer_RawDep.resize(n_elements,0);

				for (unsigned idx(0); idx<n_sens_elements; ++idx) {
					if (sens_HitVec[idx].size() > sens_HitVec_size_max) {
						sens_HitVec_size_max = 2*sens_HitVec[idx].size();
					}
					sens_HitVec[idx].clear();
					sens_HitVec[idx].reserve(sens_HitVec_size_max);
				}

			}

			inline unsigned getLayer(std::string astr) {
					size_t pos = astr.find("phys");
					unsigned num = 0;
					if (astr.find("_")== astr.npos) {
						std::string truncated = astr.substr(0,pos);
						size_t last_index = truncated.find_last_not_of("0123456789");
						num = std::atoi(truncated.substr(last_index+1).c_str());

					}
					else{
						pos = astr.find("_");
						std::string truncated = astr.substr(0,pos);
						size_t last_index = truncated.find_last_not_of("0123456789");
						num = std::atoi(truncated.substr(last_index+1).c_str());
					}

					return num;

				};


			const G4SiHitVec & getSiHitVec(const unsigned & idx) const;

			G4double getTotalSensE();

			//
			void report(bool header=false);

			//members
			unsigned n_elements;
			unsigned n_sectors;
			unsigned n_sens_elements;

			//EventAction *eventAction_;

			std::vector<G4double> sublayer_thick;
			std::vector<std::string> ele_name;
			std::vector<G4double> sublayer_X0;
			std::vector<G4double> sublayer_dEdx;
			std::vector<G4double> sublayer_L0;
			std::vector<G4double> sublayer_RawDep;
			std::vector<G4double> sublayer_PrimaryDep;
			std::vector<G4double> sublayer_dl;
			std::vector<G4VPhysicalVolume*> sublayer_vol;

			std::vector<G4SiHitVec> sens_HitVec;
			G4double Total_thick;
			bool hasScintillator;
			unsigned sens_HitVec_size_max;

		};

#endif
