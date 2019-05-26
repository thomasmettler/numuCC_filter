#ifndef MCPARTICLEHELPER_H
#define MCPARTICLEHELPER_H

#include "larcorealg/Geometry/TPCGeo.h"
#include "larcore/Geometry/Geometry.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"

#include "GeometryHelper.h"
#include "TVector3.h"

namespace constants
{
// The corners of the bottom CRT panel
const float BX1 = -130.;
const float BY = -258.;
const float BZ1 = 275.;
const float BX2 = 400.;
const float BZ2 = 800.;
} // namespace constants

struct MCParticleInfo
{
    uint process; // std::string
    uint end_process;
    bool startInside;
    bool endInside;
    bool partInside;  // This means that the track is crossing, starts/ends inside, is completely inside.
    float startX_tpc; // if the track starts outside but crosses, this will be stored here.
    float startY_tpc;
    float startZ_tpc;
    float startX_sce; // spacecharge corrected version of the tpc edge or the inside start end point.
    float startY_sce;
    float startZ_sce;
    float endX_tpc;
    float endY_tpc;
    float endZ_tpc;
    float endX_sce;
    float endY_sce;
    float endZ_sce;
    float length;
    float lengthTPC;
    float length_sce;
    float scatterdist; // Difference between the true end point and the start+dir*length point
};

struct CRTcrossing
{
    bool CRT_cross;
    float crossX;
    float crossY;
    float crossZ;
    float crossE;
    float crossT;
};

class MCParticleHelper
{
  public:
    MCParticleHelper();
    ~MCParticleHelper() = default;

    void Configure(bool is_mcc9);

    std::set<std::string> getProcesses() { return string_process; };
    MCParticleInfo fillMCP(simb::MCParticle const &mcparticle);
    CRTcrossing isCrossing(simb::MCParticle const &mcparticle);
    bool inTPC(const TVector3 x) const;
    void crossTPC(MCParticleInfo &this_mcp, const simb::MCTrajectory &traj);

  private:
    geo::BoxBoundedGeo theTpcGeo;
    art::ServiceHandle<geo::Geometry> geo_service;
    GeometryHelper geoHelper;
    bool m_is_mcc9;

    std::set<std::string> string_process; // This variable counts the different processes invloved.
    std::map<std::string, uint> map_process =
        {
            {"0", 0},
            {"CoulombScat", 1},
            {"Decay", 2},
            {"annihil", 3},
            {"compt", 4},
            {"conv", 5},
            {"dInelastic", 6},
            {"eBrem", 7},
            {"eIoni", 8},
            {"hBertiniCaptureAtRest", 9},
            {"hIoni", 10},
            {"hadElastic", 11},
            {"muBrems", 12},
            {"muIoni", 13},
            {"muMinusCaptureAtRest", 14},
            {"muPairProd", 15},
            {"muonNuclear", 16},
            {"nCapture", 17},
            {"neutronInelastic", 18},
            {"phot", 19},
            {"photonNuclear", 20},
            {"pi+Inelastic", 21},
            {"pi-Inelastic", 22},
            {"primary", 23},
            {"CoupledTransportation", 24},
            {"FastScintillation", 25},
            {"protonInelastic", 26},
            {"anti_protonInelastic", 27},
            {"electronNuclear", 28}};
};

MCParticleHelper::MCParticleHelper()
{
    //// Set up the boxboundedgeo:
    geo::TPCGeo const &thisTPC = geo_service->TPC();
    theTpcGeo = thisTPC.ActiveBoundingBox();
    std::cout << std::endl;
    std::cout << "[MCParticleHelper constructor] ActivePTC volume max: " << theTpcGeo.MaxX() << ", " << theTpcGeo.MaxY() << ", " << theTpcGeo.MaxZ() << std::endl;
    std::cout << "[MCParticleHelper constructor] ActivePTC volume min: " << theTpcGeo.MinX() << ", " << theTpcGeo.MinY() << ", " << theTpcGeo.MinZ() << std::endl;
    //// Check if spacecharge correction is working.
    auto const &SCE(*lar::providerFrom<spacecharge::SpaceChargeService>());
    auto scecorr = SCE.GetPosOffsets(geo::Point_t(theTpcGeo.MinX(), theTpcGeo.MinY(), theTpcGeo.MinZ()));
    double xOffset = scecorr.X();
    double yOffset = scecorr.Y();
    double zOffset = scecorr.Z();
    std::cout << "[MCParticleHelper constructor] Spacecharge correction at lower TPC corner: " << xOffset << ", " << yOffset << ", " << zOffset << std::endl;
    scecorr = SCE.GetPosOffsets(geo::Point_t(theTpcGeo.MaxX(), theTpcGeo.MaxY(), theTpcGeo.MaxZ()));
    xOffset = scecorr.X();
    yOffset = scecorr.Y();
    zOffset = scecorr.Z();
    std::cout << "[MCParticleHelper constructor] Spacecharge correction at uper TPC corner: " << xOffset << ", " << yOffset << ", " << zOffset << std::endl;
    scecorr = SCE.GetPosOffsets(geo::Point_t(0., -116.5, 0.));
    xOffset = scecorr.X();
    yOffset = scecorr.Y();
    zOffset = scecorr.Z();
    std::cout << "[MCParticleHelper constructor] Spacecharge correction at (0., -116.5, 0.): " << xOffset << ", " << yOffset << ", " << zOffset << std::endl;
    scecorr = SCE.GetPosOffsets(geo::Point_t(256.35, 116.5, 1036.8));
    xOffset = scecorr.X();
    yOffset = scecorr.Y();
    zOffset = scecorr.Z();
    std::cout << "[MCParticleHelper constructor] Spacecharge correction at (256.35., 116.5, 1036.8): " << xOffset << ", " << yOffset << ", " << zOffset << std::endl;
    std::cout << std::endl;
}

void MCParticleHelper::Configure(bool is_mcc9)
{
    m_is_mcc9 = is_mcc9;
    std::cout << "[MCParticleHelper configuration] File is mcc9? " << m_is_mcc9 << std::endl;
}

MCParticleInfo MCParticleHelper::fillMCP(simb::MCParticle const &mcparticle)
{
    MCParticleInfo this_mcp;

    string_process.insert(mcparticle.Process());
    if (map_process.find(mcparticle.Process()) != map_process.end())
    {
        this_mcp.process = map_process[mcparticle.Process()];
    }
    else
    {
        this_mcp.process = 999;
        std::cout << "[CosmicStudies::analyze] New MC interaction process found: " << mcparticle.Process() << std::endl;
    }

    string_process.insert(mcparticle.EndProcess());
    if (map_process.find(mcparticle.EndProcess()) != map_process.end())
    {
        this_mcp.end_process = map_process[mcparticle.EndProcess()];
    }
    else
    {
        this_mcp.end_process = 999;
        std::cout << "[CosmicStudies::analyze] New MC interaction end_process found!" << std::endl;
    }

    TVector3 mc_start = mcparticle.Position().Vect();
    TVector3 mc_end = mcparticle.EndPosition().Vect();
    this_mcp.length = (mc_start - mc_end).Mag(); // This is the total length, not the length in the detector!
    this_mcp.startInside = inTPC(mc_start);
    this_mcp.endInside = inTPC(mc_end);

    if (this_mcp.startInside)
    {
        this_mcp.startX_tpc = mc_start.X();
        this_mcp.startY_tpc = mc_start.Y();
        this_mcp.startZ_tpc = mc_start.Z();
    }
    if (this_mcp.endInside)
    {
        this_mcp.endX_tpc = mc_end.X();
        this_mcp.endY_tpc = mc_end.Y();
        this_mcp.endZ_tpc = mc_end.Z();
    }
    if (!this_mcp.startInside || !this_mcp.endInside)
    {
        crossTPC(this_mcp, mcparticle.Trajectory());
    }
    else
    {
        this_mcp.partInside = true;
    }

    /*
    this_mcp.partInside = true;
    //Find the section that is inside the tpc
    if (!this_mcp.startInside || !this_mcp.endInside)
    {
        TVector3 mc_startdir = mcparticle.Momentum().Vect();
        std::vector<TVector3> intersections = theTpcGeo.GetIntersections(mc_start, mc_startdir);
        uint num_intersections = intersections.size();

        //Particles completely passes the TPC without entering
        if (num_intersections == 0)
        {
            this_mcp.partInside = false;
            this_mcp.lengthTPC = 0;
            this_mcp.length_sce = 0;
        }

        // Particle started inside the TPC and is exiting
        else if (num_intersections == 1)
        {
            this_mcp.startX_tpc = mc_start.X();
            this_mcp.startY_tpc = mc_start.Y();
            this_mcp.startZ_tpc = mc_start.Z();
            this_mcp.endX_tpc = intersections[0].X();
            this_mcp.endY_tpc = intersections[0].Y();
            this_mcp.endZ_tpc = intersections[0].Z();
        }

        // Particle crosses TPC or particle stops in TPC or particle stops before TPC
        else if (num_intersections == 2)
        {
            float len_start_cross0 = (mc_start - intersections[0]).Mag();
            float len_start_cross1 = (mc_start - intersections[1]).Mag();

            // Particle is crossing the TPC
            if (std::max(len_start_cross0, len_start_cross1) < this_mcp.length)
            {
                this_mcp.startX_tpc = intersections[0].X();
                this_mcp.startY_tpc = intersections[0].Y();
                this_mcp.startZ_tpc = intersections[0].Z();
                this_mcp.endX_tpc = intersections[1].X();
                this_mcp.endY_tpc = intersections[1].Y();
                this_mcp.endZ_tpc = intersections[1].Z();
            }
            // Particle stops before entering
            else if (std::min(len_start_cross0, len_start_cross1) > this_mcp.length)
            {
                this_mcp.partInside = false;
                this_mcp.lengthTPC = 0;
                this_mcp.length_sce = 0;
            }
            // Particle stops inside
            else
            {
                this_mcp.endX_tpc = mc_end.X();
                this_mcp.endY_tpc = mc_end.Y();
                this_mcp.endZ_tpc = mc_end.Z();
                if (len_start_cross0 < this_mcp.length)
                {
                    this_mcp.startX_tpc = intersections[0].X();
                    this_mcp.startY_tpc = intersections[0].Y();
                    this_mcp.startZ_tpc = intersections[0].Z();
                }
                else
                {
                    this_mcp.startX_tpc = intersections[1].X();
                    this_mcp.startY_tpc = intersections[1].Y();
                    this_mcp.startZ_tpc = intersections[1].Z();
                }
            }
        }
    }
    else //Start and end are inside
    {
        this_mcp.startX_tpc = mc_start.X();
        this_mcp.startY_tpc = mc_start.Y();
        this_mcp.startZ_tpc = mc_start.Z();
        this_mcp.endX_tpc = mc_end.X();
        this_mcp.endY_tpc = mc_end.Y();
        this_mcp.endZ_tpc = mc_end.Z();
    }
    */

    if (this_mcp.partInside)
    {
        std::vector<float> start_tpc = {this_mcp.startX_tpc, this_mcp.startY_tpc, this_mcp.startZ_tpc};
        std::vector<float> end_tpc = {this_mcp.endX_tpc, this_mcp.endY_tpc, this_mcp.endZ_tpc};
        this_mcp.lengthTPC = geoHelper.distance(start_tpc, end_tpc);

        //Correct the inside tpcpoints for spacecharge

        auto const &SCE(*lar::providerFrom<spacecharge::SpaceChargeService>());
        auto sce_start = SCE.GetPosOffsets(geo::Point_t(this_mcp.startX_tpc, this_mcp.startY_tpc, this_mcp.startZ_tpc));
        auto sce_end = SCE.GetPosOffsets(geo::Point_t(this_mcp.endX_tpc, this_mcp.endY_tpc, this_mcp.endZ_tpc));

        this_mcp.startY_sce = this_mcp.startY_tpc + sce_start.Y();
        this_mcp.startZ_sce = this_mcp.startZ_tpc + sce_start.Z();
        this_mcp.endY_sce = this_mcp.endY_tpc + sce_end.Y();
        this_mcp.endZ_sce = this_mcp.endZ_tpc + sce_end.Z();

        if (m_is_mcc9)
        {
            auto const &detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
            auto const &detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
            float g4Ticks = detClocks->TPCG4Time2Tick(mcparticle.T()) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
            float xtimeoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);
            this_mcp.startX_sce = (this_mcp.startX_tpc + xtimeoffset + sce_start.X()) * (1.114 / 1.098) + 0.6;
            this_mcp.endX_sce = (this_mcp.endX_tpc + xtimeoffset + sce_end.X()) * (1.114 / 1.098) + 0.6;
        }
        else
        {
            this_mcp.startX_sce = this_mcp.startX_tpc - sce_start.X() + 0.7;
            this_mcp.endX_sce = this_mcp.endX_tpc - sce_end.X() + 0.7;
        }

        std::vector<float> start_sce = {this_mcp.startX_sce, this_mcp.startY_sce, this_mcp.startZ_tpc};
        std::vector<float> end_sce = {this_mcp.endX_sce, this_mcp.endY_sce, this_mcp.endZ_sce};
        this_mcp.length_sce = geoHelper.distance(start_sce, end_sce);
    }
    else // not part inside
    {
        this_mcp.lengthTPC = 0;
        this_mcp.length_sce = 0;
    }

    // Calculate the scatterdist:
    this_mcp.scatterdist = (mcparticle.Position().Vect() + mcparticle.Momentum().Vect().Unit() * this_mcp.length - mcparticle.EndPosition().Vect()).Mag();

    return this_mcp;
}

CRTcrossing MCParticleHelper::isCrossing(simb::MCParticle const &mcparticle)
{
    CRTcrossing this_cross;
    this_cross.CRT_cross = false;

    const simb::MCTrajectory &traj = mcparticle.Trajectory();
    TVector3 pt1;
    TVector3 pt2;
    float t;

    for (size_t i = 1; i < traj.size(); i++)
    {
        pt1 = traj.Position(i - 1).Vect();
        pt2 = traj.Position(i).Vect();

        t = (constants::BY - pt1.Y()) / (pt2.Y() - pt1.Y());
        // if t < 0 or > 1 then the intersection is beyond the segment
        if ((t > 0) && (t <= 1))
        {
            this_cross.crossX = pt1.X() + (pt2.X() - pt1.X()) * t;
            this_cross.crossZ = pt1.Z() + (pt2.Z() - pt1.Z()) * t;
            if ((this_cross.crossX > constants::BX1) && (this_cross.crossX < constants::BX2) && (this_cross.crossZ > constants::BZ1) && (this_cross.crossZ < constants::BZ2))
            {
                this_cross.crossY = pt1.Y() + (pt2.Y() - pt1.Y()) * t;
                this_cross.crossE = traj.E(t - 1);
                this_cross.crossT = traj.T(t - 1);
                this_cross.CRT_cross = true;
                return this_cross;
            } // if they intersec
        }     // if t is between 0 and 1
    }
    return this_cross;
}

bool MCParticleHelper::inTPC(const TVector3 x) const
{
    bool is_x = x.X() >= theTpcGeo.MinX()  && x.X() <= theTpcGeo.MaxX() ;
    bool is_y = x.Y() >= theTpcGeo.MinY()  && x.Y() <= theTpcGeo.MaxY() ;
    bool is_z = x.Z() >= theTpcGeo.MinZ()  && x.Z() <= theTpcGeo.MaxZ() ;
    return is_x && is_y && is_z;
}

void MCParticleHelper::crossTPC(MCParticleInfo &this_mcp, const simb::MCTrajectory &traj)
{
    this_mcp.partInside = this_mcp.startInside;

    TVector3 pt1 = traj.Position(0).Vect();
    bool pt1_inTPC =inTPC(pt1);

    for (size_t i = 1; i < traj.size(); i++)
    {
        TVector3 pt2 = traj.Position(i).Vect();
        bool pt2_inTPC =inTPC(pt2);

        if (pt1_inTPC != pt2_inTPC) // Crossing the TPC!
        {
            // std::cout << "[MCParticleHelper::crossTPC] " << i << "\t" << pt2.X() << "\t" << pt2.Y() << "\t" << pt2.Z() << std::endl;
            TVector3 pt1_dir = pt2-pt1;
            std::vector<TVector3> intersections = theTpcGeo.GetIntersections(pt1, pt1_dir);
            if (intersections.size() < 2)
            {
                std::cout << "[MCParticleHelper::crossTPC] crossing expected but only one found!" << std::endl;
                std::cout << "[MCParticleHelper::crossTPC] " << pt1.X() << "\t" << pt1.Y() << "\t" << pt1.Z() << std::endl;
                std::cout << "[MCParticleHelper::crossTPC] " << pt2.X() << "\t" << pt2.Y() << "\t" << pt2.Z() << std::endl;
            }
            if (!this_mcp.partInside)             // Entering the TPC
            {
                TVector3 mc_cross = intersections[0]; // Take the first crossing point.
                this_mcp.partInside = true;
                this_mcp.startX_tpc = mc_cross.X();
                this_mcp.startY_tpc = mc_cross.Y();
                this_mcp.startZ_tpc = mc_cross.Z();
                if (this_mcp.endInside)
                {
                    return; // The end point is already filled in.
                }
            }
            else // Exiting the TPC
            {
                TVector3 mc_cross = intersections[1]; // Take the first crossing point.
                this_mcp.endX_tpc = mc_cross.X();
                this_mcp.endY_tpc = mc_cross.Y();
                this_mcp.endZ_tpc = mc_cross.Z();
                // std::cout << "[MCParticleHelper::crossTPC] " << mc_cross.X() << "\t" << mc_cross.Y() << "\t" << mc_cross.Z() << std::endl;
                return;
            }
        }
        pt1 = pt2;
        pt1_inTPC = pt2_inTPC;
    }
}

#endif