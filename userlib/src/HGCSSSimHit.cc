#include "HGCSSSimHit.hh"

#include <iomanip>
#include <cmath>
#include <stdlib.h>

HGCSSSimHit::HGCSSSimHit(const G4SiHit & aSiHit,
		TH2Poly* map, const float) {
	double x = aSiHit.hit_x;
	double y = aSiHit.hit_y;
	assert(map);

	nGammas_ = 0;
	nElectrons_ = 0;
	nMuons_ = 0;
	nNeutrons_ = 0;
	nProtons_ = 0;
	nHadrons_ = 0;

	if (abs(aSiHit.pdgId) == 22)
		nGammas_++;
	else if (abs(aSiHit.pdgId) == 11)
		nElectrons_++;
	else if (abs(aSiHit.pdgId) == 13)
		nMuons_++;
	else if (abs(aSiHit.pdgId) == 2112)
		nNeutrons_++;
	else if (abs(aSiHit.pdgId) == 2212)
		nProtons_++;
	else
		nHadrons_++;

	cellid_ = map->FindBin(x, y);
	if (cellid_ > 1e5)
		std::cout << "The cell id = " << cellid_ << ", (x,y) = " << "(" << x << ", " << y << ")" << std::endl;
	layer_ = aSiHit.layer;
	energy_ = aSiHit.energyDep;
	trackIDMainParent_ = aSiHit.trackId;
	KEMainParent_ = aSiHit.parentKE;
	pdgIDMainParent_ = aSiHit.pdgId;
	eDepMainParent_ = aSiHit.energyDep ;
}

void HGCSSSimHit::Add(const G4SiHit & aSiHit) {

	if (abs(aSiHit.pdgId) == 22)
		nGammas_++;
	else if (abs(aSiHit.pdgId) == 11)
		nElectrons_++;
	else if (abs(aSiHit.pdgId) == 13)
		nMuons_++;
	else if (abs(aSiHit.pdgId) == 2112)
		nNeutrons_++;
	else if (abs(aSiHit.pdgId) == 2212)
		nProtons_++;
	else
		nHadrons_++;

	energy_ += aSiHit.energyDep;

	if (aSiHit.energyDep > eDepMainParent_) {
		trackIDMainParent_ = aSiHit.parentId;
		eDepMainParent_ = aSiHit.energyDep;
		pdgIDMainParent_ = aSiHit.pdgId;
		KEMainParent_ = aSiHit.parentKE;

	}

}



std::pair<double, double> HGCSSSimHit::get_xy(const bool isScintillator,
		const HGCSSGeometryConversion & aGeom) const {
	if (isScintillator)
		return aGeom.squareGeom.find(cellid_)->second;
	else
		return aGeom.hexaGeom.find(cellid_)->second;

}

ROOT::Math::XYZPoint HGCSSSimHit::position(const bool isScintillator,
		const HGCSSGeometryConversion & aGeom) const {
	std::pair<double, double> xy = get_xy(isScintillator, aGeom);
	return ROOT::Math::XYZPoint(xy.first / 10., xy.second / 10., zpos_ / 10.);
}

double HGCSSSimHit::theta(const bool isScintillator,
		const HGCSSGeometryConversion & aGeom) const {
	return 2 * atan(exp(-1. * eta(isScintillator, aGeom)));
}

double HGCSSSimHit::eta(const bool isScintillator,
		const HGCSSGeometryConversion & aGeom) const {
	return position(isScintillator, aGeom).eta();
}

double HGCSSSimHit::phi(const bool isScintillator,
		const HGCSSGeometryConversion & aGeom) const {
	return position(isScintillator, aGeom).phi();
}

void HGCSSSimHit::Print(std::ostream & aOs) const {
	aOs << "====================================" << std::endl << " = Layer "
			<< layer() << " siLayer " << silayer() << " cellid " << cellid_
			<< std::endl << " = Energy " << energy_ << " time " << time_
			<< std::endl << " = g " << nGammas_ << " e " << nElectrons_
			<< " mu " << nMuons_ << " neutron " << nNeutrons_ << " proton "
			<< nProtons_ << " had " << nHadrons_ << std::endl
			<< " = main parent: trackID " << trackIDMainParent_ << " efrac "
			<< mainParentEfrac() << std::endl
			<< "====================================" << std::endl;

}

