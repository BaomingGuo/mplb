/**
 * Copyright 2019 United Kingdom Research and Innovation
 *
 * Authors: See AUTHORS
 *
 * Contact: [jianping.meng@stfc.ac.uk and/or jpmeng@gmail.com]
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice
 *    this list of conditions and the following disclaimer in the documentation
 *    and or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * ANDANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef FREESURFACE2D_KERNEL_INC
#define FREESURFACE2D_KERNEL_INC
#include "ops_lib_core.h"

#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "model.h"
#include "type.h"
#include "flowfield_host_device.h"
#include "OWC2d_host_device.h"
#include "point_position.h"
#define VertexType_Interface 2000
#define VertexType_Empty 2001
#define VertexType_NewInterface 2003
#define VertexType_Interface2Empty 2004
#define VertexType_Interface2Fluid 2005
#define VertexType_SlipWall 2006
#define VertexType_Empty2Solid 2007
#define VertexType_Solid2Interface 2008
#define VertexType_Solid2Fluid 2009
const int iMax {3027-1};
int jMax {0};
const int Finaldepth{80};
const int BoundaryScheme_HalfwayBounceBack{100};
const int BoundaryScheme_HalfwayBounceForward{101};
int jofInletInterface{jMax};
const Real VOFofInletInterface{0.8};
const Real DistanceforUIBBWall{1-(VOFofInletInterface-0.5)};
bool MultiLayerConvertedSolid{false};
bool ConsistentInitialize{false};
#define FreeMpi_2D
enum VertexQuality { No_Fluid = 1, No_Empty = 2, Standard = 3, Solid = -1 };
enum VertexPosition { In_Chamber = 1, Out_Chamber = 2, ImmerseSolid = -1};
#ifdef OPS_2D
    //  #define CIRCLE
    Real Radii{0.5};
    Real IntialCenter[2]{1.6025, 0.8755};
    Real RotationCenter[2]{1.6025, 0.8755};
    Real TotalHydroForce[6]{};
    Real TotalAirForce[6]{};
    Real VelocityofRoCenter[3]{0,0,0};
    Real RotationCenterOld[2]{1.6025, 0.8755};
    Real VelocityofRoCenterOld[3]{0,0,0};
    Real TotalHydroForceOld[6]{};
#endif

#ifdef OPS_3D
    Real RotationCenter[3]{};
    Real TotalHydroForce[6]{};
    Real VelocityofRoCenter[6]{};
#endif

// define the square solid boundary and vertices
#ifdef OPS_2D
    //  2********1
    //  *        *
    //  *        * 
    //  *        *
    //  3********0
    #define SQUARE
    const int dim{2};
    const int Polygon {4};
    const int NumofPol{2};
    const int NumofChamber{0};
    const int BasicUnit{2};//two vertex in 2D make a line
    Real SolidVertex[]{ 2.0025, 0.4975, 2.0025, 1.2535, 
                        1.2025, 1.2535, 1.2025, 0.4975,

                        1.8025, 0.4975, 1.8025, 1.2015, 
                        1.4025, 1.2015, 1.4025, 0.4975};
    Real SolidPlaneNormal[]{ 1,  0, 0,  1,
                            -1,  0, 0, -1,

                             1,  0, 0,  1,
                            -1,  0, 0, -1};
    Real SolidVertexOld[NumofPol * Polygon * dim]{0};
    Real SolidPlaneNormalOld[NumofPol * Polygon * dim]{0};
    const Real GeometryPar[]{0.8, 0.756, 0.4, 0.704};
#endif


#define OWCOn
    Real DeltaRho{0};  // prepare for chamber pressure
    Real DeltaRhoOld{0};
#ifdef OWCOn
    Real VOChamber{0.08};
    Real VOChamberOld{0.08};
    Real POld{0};
#endif

    void KerSetInitialMacroVars(ACC<Real>& rho, ACC<Real>& u, ACC<Real>& v,
                                ACC<Real>& mass, ACC<Real>& vof, 
                                const ACC<int>& nodeType,
                                const int* idx,
                                const ACC<Real>& acceleration,
                                const ACC<Real>& Porosity) {
        const int vt{nodeType(0, 0)};
        rho(0, 0) =
            1 + CS * TimeStep() * (idx[1] - (jMax +VOFofInletInterface-0.5) ) * acceleration(1, 0, 0);
        //+MeshSize*(j-Finaldepth)*acceleration(1,0,0);//+(jMax-0.5)*MeshSize*Gravity;
        u(0, 0) = 0;
        v(0, 0) = 0;
        // const int i{idx[0]};
        //  Initialise the mass and vof here according to x,y
        if ((int)VertexType::Fluid == vt) {
            vof(0, 0) = 1.0;
        } else if (VertexType_Interface == vt) {
            vof(0, 0) = VOFofInletInterface;
        } else {
            vof(0, 0) = 0.0;
        }
        mass(0, 0) = vof(0, 0) * rho(0, 0)*Porosity(0,0);
        // porosity is key for sovling the vof of interface cell
}

void KerSetInitialNodeType(ACC<int>& nodeType, 
                           const ACC<Real>& coordinates,const int* idx) {
    const int j{idx[1]};
    const int i{idx[0]};
    nodeType(0, 0) = VertexType_Empty;
    // dam breaking from left boundary
    if ((j < jMax && j > 0)) {
        nodeType(0, 0) = (int)VertexType::Fluid;
    } else if ((jMax == j)) {
        nodeType(0, 0) = VertexType_Interface;
    } else {
        nodeType(0, 0) = VertexType_Empty;
    }
}

void KerInitialiseBGK2ndfree(ACC<Real>& f, const ACC<int>& nodeType,
                         const ACC<Real>& Rho, const ACC<Real>& U,
                         const ACC<Real>& V, const int* lattIdx,
                         const ACC<Real>& acceleration) {
#ifdef OPS_2D
    int vt = nodeType(0, 0);
    if (vt != (int)VertexType::ImmersedSolid && 
        VertexType_Empty != vt) {
        Real rho{Rho(0, 0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)-acceleration(1,0,0)*TimeStep()/2};
        const Real T{1};
        const int polyOrder{2};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            f(xiIdx, 0, 0) = CalcBGKFeq(xiIdx, rho, u, v, T, polyOrder);
#ifdef CPU
            const Real res{f(xiIdx, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function %f becomes invalid at the "
                    "lattice %i\n",
                    res, xiIdx);
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}

void KerUpdateMacroBodyForce(ACC<Real>& acceleration,
                             const ACC<Real>& coordinates, const int* idx) {
    // here we put some constant for gravity to the following two lines
    acceleration(1, 0, 0) = 0;//-1.072e-4;
    acceleration(0, 0, 0) = 0;
}

void KerUpdateMacroBodyForcetimeVOF(ACC<Real>& acceleration,
                             const ACC<Real>& coordinates, 
                             const int* idx,const ACC<Real>& vof) {
    // here we put some constant for gravity to the following two lines
    acceleration(1, 0, 0) = acceleration(1, 0, 0)*vof(0,0);
    acceleration(0, 0, 0) = acceleration(0, 0, 0)*vof(0,0);
}


// we generate wave and absorb reflected wave by momentum source
void KerGenerateWave(ACC<Real>& acceleration,const ACC<Real>& coordinates, 
                     const int* idx,const Real* time) {
    const Real iCoordinate{coordinates(0,0,0)};
    Real SourceWidth{WaveLength};
    // we always assume: SourceWidth = WaveLength,
    // wave absorption length = 2*WaveLength
    // 5 from half SourceWidth and sponge width (2*L)
    const Real CenterofGArea = BlockCoor[0] + DampWidth + 0.5*SourceWidth/L0;
    if (iCoordinate >= CenterofGArea - 0.5*SourceWidth/L0 &&
        iCoordinate <= CenterofGArea + 0.5*SourceWidth/L0) {
        Real alpha{-0.38955};
        Real WaveNumber = 2 * PI / WaveLength;
        Real WaveFrequency = 2 * PI / WavePeriod;
        Real beta = 20 / (SourceWidth * SourceWidth);
        Real I = sqrt(PI / beta) * exp(-WaveNumber * WaveNumber / (4 * beta));
        Real SourceStrength = 2 * WaveAmplitude *
                              (pow(WaveFrequency, 2) -
                               (alpha + 1.0 / 3.0) * 9.81 * pow(WaveNumber, 4) *
                                   pow(WaterDepth, 3)) /
                              (WaveFrequency * I * WaveNumber *
                               (1 - alpha * pow(WaveNumber * WaterDepth, 2)));
        Real Distance = (iCoordinate - CenterofGArea)*L0;
        acceleration(0, 0, 0) += L0/pow(SoundofWater,2) * 9.81 * 2 * beta * Distance *
                                 exp(-beta * Distance * Distance) *
                                 SourceStrength / WaveFrequency *
                                 sin(WaveFrequency * (*time) * (L0/ SoundofWater));
    }  // generate waves
}

void KerAbsorbWave(ACC<Real>& acceleration, const ACC<Real>& coordinates,
                   const int* idx,
                   const ACC<Real>& U, const ACC<Real>& V) {
    // we generate wave and absorb reflected wave by momentum source
    const Real iCoordinate{coordinates(0,0,0)};
    if (iCoordinate <= BlockCoor[0] + DampWidth|| 
        iCoordinate >= BlockCoor[1] - DampWidth ) {
        // || iCoordinate >= *endCoordinate - DampWidth
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        //Real WaveAmplitude{0.01};
        Real LinearConstant{LinearConst};  // unit 1/T
        Real NonLinearConst{0 * L0};                 // unit 1/L                             // empirical coefficient
        // we always assume: SourceWidth = WaveLength,
        // wave absorption length = 2*WaveLength
        Real StartPt{0};
        if (iCoordinate <= BlockCoor[0] + DampWidth) {
            StartPt = BlockCoor[0] + DampWidth;
        } else if (iCoordinate >= BlockCoor[1] - DampWidth) {
            StartPt = BlockCoor[1] - DampWidth;
        }
        Real DistanceScale = (iCoordinate - StartPt) / DampWidth;
        acceleration(0, 0, 0) +=
            -(LinearConstant * u + NonLinearConst * u * sqrt(u * u + v * v)) *
            (exp(pow(DistanceScale, ns)) - 1) / (exp(1) - 1);
        acceleration(1, 0, 0) +=
            -(LinearConstant * v + NonLinearConst * v * sqrt(u * u + v * v)) *
            (exp(pow(DistanceScale, ns)) - 1) / (exp(1) - 1);
    }
}

void KerCalcDensityfree(ACC<Real>& Rho, const ACC<Real>& f,
                    const ACC<int>& nodeType, const int* lattIdx) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    if (vt != (int)VertexType::ImmersedSolid && 
        VertexType_Empty != vt) {
        Real rho{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            rho += f(xiIdx, 0, 0);
        }
#ifdef CPU
        if (isnan(rho) || rho <= 0 || isinf(rho)) {
            ops_printf(
                "Error! Density %f becomes invalid！Something "
                "wrong...",
                rho);
            assert(!(isnan(rho) || rho <= 0 || isinf(rho)));
        }
#endif
        Rho(0, 0) = rho;
    }
#endif  // OPS_2D
}

void KerCalcUForcefree(ACC<Real>& U, const ACC<Real>& f, const ACC<int>& nodeType,
                   const ACC<Real>& coordinates, const ACC<Real>& acceleration,
                   const ACC<Real>& Rho, const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    const Real x{coordinates(0, 0, 0)};
    const Real y{coordinates(1, 0, 0)};
    const int vt {nodeType(0, 0)};
    if (vt != (int) VertexType::ImmersedSolid && 
        VertexType_Empty != vt) {
        Real u{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            u += CS * XI[xiIdx * LATTDIM] * f(xiIdx, 0, 0);
        }
        u /= Rho(0, 0);
        if ((int)VertexType::Fluid == vt || (int) VertexType::MDPeriodic == vt ||
            VertexType_Interface == vt||
            (int)VertexType::Inlet == vt||
            (int)VertexType::OutLet == vt) {
            u += ((*dt) * acceleration(0, 0, 0) / 2);
        }
#ifdef CPU
        if (isnan(u) || isinf(u)) {
            ops_printf(
                "Error! Velocity U=%f becomes invalid! Maybe something wrong "
                "at x=%f y=%f\n",
                u, x, y);
            assert(!(isnan(u) || isinf(u)));
        }
#endif
        U(0, 0) = u;
    }
#endif  // OPS_2D
}

void KerCalcVForcefree(ACC<Real>& V, const ACC<Real>& f, const ACC<int>& nodeType,
                   const ACC<Real>& coordinates, const ACC<Real>& acceleration,
                   const ACC<Real>& Rho, const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    const Real x{coordinates(0, 0, 0)};
    const Real y{coordinates(1, 0, 0)};
    const int vt {nodeType(0, 0)};
    if (vt != (int) VertexType::ImmersedSolid && 
        VertexType_Empty != vt) {
        Real v{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            v += CS * XI[xiIdx * LATTDIM + 1] * f(xiIdx, 0, 0);
        }
        v /= Rho(0, 0);
        if ((int)VertexType::Fluid == vt || (int) VertexType::MDPeriodic == vt||
            VertexType_Interface == vt||
            (int)VertexType::Inlet == vt||
            (int)VertexType::OutLet == vt) {
            v += ((*dt) * acceleration(1, 0, 0) / 2);
        }
#ifdef CPU
        if (isnan(v) || isinf(v)) {
            ops_printf(
                "Error! Velocity V=%f becomes invalid! Maybe something wrong "
                "at x=%f y=%f\n",
                v, x, y);
            assert(!(isnan(v) || isinf(v)));
        }
#endif
        V(0, 0) = v;
    }
#endif  // OPS_2D
}

void KerCalcBodyForce1STfree(ACC<Real>& fStage, const ACC<Real>& acceleration,
                         const ACC<Real>& Rho, const ACC<int>& nodeType,
                         const int* lattIdx, const ACC<Real>& U,const ACC<Real>& V) {
#ifdef OPS_2D

    const int vt {nodeType(0, 0)};
    if (vt == (int)VertexType::Fluid || 
        vt == (int)VertexType::MDPeriodic||
        VertexType_Interface == vt||
        (int)VertexType::Inlet == vt||
        (int)VertexType::OutLet == vt) {
        Real rho{Rho(0, 0)};
        Real g[]{acceleration(0, 0, 0), acceleration(1, 0, 0)};
        Real u{U(0,0)};
        Real v{V(0,0)};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real bodyForce{CalcBodyForcefree(xiIndex, rho, g,u,v)};
#ifdef CPU
            if (isnan(bodyForce) || isinf(bodyForce)) {
                ops_printf(
                    "Error! Body force  %f becomes invalid  at  the lattice "
                    "%i\n",
                    bodyForce, xiIndex);
                assert(!(isnan(bodyForce) || isinf(bodyForce)));
            }
#endif
            fStage(xiIndex, 0, 0) = bodyForce;
        }
    }
#endif  // OPS_2D
}

//only add "VertexType_Interface == vt" into if statement
void KerStreamfree(ACC<Real>& f, const ACC<Real>& fStage, const ACC<int>& nodeType,
               const ACC<int>& geometry, const int* lattIdx) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0);
    for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
        int cx = (int)XI[xiIndex * LATTDIM];
        int cy = (int)XI[xiIndex * LATTDIM + 1];

        if (vt == (int)VertexType::Fluid || 
            vt == (int)VertexType::MDPeriodic ||
            vt == (int)VertexType::VirtualBoundary||
            VertexType_Interface == vt) {
            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
            continue;
        }
        // Block boundary and immersed boundary
        if (vt != (int)VertexType::ImmersedSolid &&
            vt != (int)VertexType::Fluid &&
            vt != (int)VertexType::VirtualBoundary &&
            vt != (int)VertexType::MDPeriodic && 
            VertexType_Interface != vt &&
            (int)VertexType::Wall != vt && 
            VertexType_Empty != vt) {
            // this is reserved in case some boundary conditions
            // might not allow stream
            bool streamRequired{false};

            if (streamRequired) {
                if ((cx == 0) && (cy == 0)) {
                    f(xiIndex, 0, 0) = fStage(xiIndex, 0, 0);
                    continue;
                }
            }
            switch (vg) {
                case VG_IP:
                    // (cx=0 means stream is implemented at i=0,so here we
                    //  disable the step at boundary)
                    if (streamRequired) {
                        if (cx <= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cx < 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IM:
                    if (streamRequired) {
                        if (cx >= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cx > 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_JP:
                    if (streamRequired) {
                        if (cy <= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy < 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_JM:
                    if (streamRequired) {
                        if (cy >= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy > 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IPJP_I:
                    if (streamRequired) {
                        if (cy <= 0 && cx <= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy < 0 && cx < 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IPJM_I:
                    if (streamRequired) {
                        if (cy >= 0 && cx <= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy > 0 && cx < 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IMJP_I:
                    if (streamRequired) {
                        if (cy <= 0 && cx >= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy < 0 && cx > 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IMJM_I:
                    if (streamRequired) {
                        if (cy >= 0 && cx >= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy > 0 && cx > 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IPJP_O:
                    if (streamRequired) {
                        if (cy <= 0 || cx <= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy < 0 || cx < 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IPJM_O:
                    if (streamRequired) {
                        if (cy >= 0 || cx <= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy > 0 || cx < 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IMJP_O:
                    if (streamRequired) {
                        if (cy <= 0 || cx >= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy < 0 || cx > 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                case VG_IMJM_O:
                    if (streamRequired) {
                        if (cy >= 0 || cx >= 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    } else {
                        if (cy > 0 || cx > 0) {
                            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
                        }
                    }
                    break;
                default:
                    break;
            }
        }
    }
#endif  // OPS_2D
}

void KerCollideBGKFreeSurface(ACC<Real>& fStage, const ACC<Real>& f,
                             const ACC<Real>& coordinates,
                             const ACC<int>& nodeType, const ACC<Real>& Rho,
                             const ACC<Real>& U, const ACC<Real>& V,
                             const Real* tauRef, const Real* dt,
                             const int* lattIdx,ACC<Real>&Taucopy,
                             const int* idx, const int* blocksizeX) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != (int)VertexType::ImmersedSolid);
    if (collisionRequired) {
        Real rho{Rho(0, 0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        const Real T{1};
        const int polyOrder{2};
        Real tau = (*tauRef);
        Real fstream[lattIdx[1]+1]{0};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            fstream[xiIndex] = f(xiIndex,0,0);
        }
        tau = Smagorinsky_LES(rho, u, v, tau, fstream, lattIdx);
        // Turbulence model

        Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real feq{CalcBGKFeq(xiIndex, rho, u, v, T, polyOrder)};
            if (vt == (int)VertexType::Fluid || 
                vt == (int)VertexType::MDPeriodic||
                VertexType_Interface== vt||
                (int)VertexType::Inlet == vt||
                (int)VertexType::OutLet == vt) {
                fStage(xiIndex, 0, 0) =
                    feq + (1 - dtOvertauPlusdt) * (f(xiIndex, 0, 0) - feq) +
                    tau * dtOvertauPlusdt * fStage(xiIndex, 0, 0);
            } else {
                fStage(xiIndex, 0, 0) =
                    feq + (1 - dtOvertauPlusdt) * (f(xiIndex, 0, 0) - feq);
            }
#ifdef CPU
            const Real res{fStage(xiIndex, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where feq=%e and rho=%e u=%e v=%e at "
                    "x=%e y=%e\n",
                    res, xiIndex, feq, rho, u, v, coordinates(0, 0, 0),
                    coordinates(1, 0, 0));
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}

// MRT model
Real TransformationM[9][9]{
                { 1,  1,  1,  1,  1,  1,  1,  1,  1},
                {-4,  2, -1,  2, -1,  2, -1,  2, -1},
                { 4,  1, -2,  1, -2,  1, -2,  1, -2},
                { 0, -1, -1, -1,  0,  1,  1,  1,  0},
                { 0, -1,  2, -1,  0,  1, -2,  1,  0},
                { 0,  1,  0, -1, -1, -1,  0,  1,  1},
                { 0,  1,  0, -1,  2, -1,  0,  1, -2},
                { 0,  0,  1,  0, -1,  0,  1,  0, -1},
                { 0, -1,  0,  1,  0, -1,  0,  1,  0}};
Real InverseM [9][9] {
             {1.0/9, -1.0/9 ,  1.0/9 ,  0    ,  0     ,  0    ,  0     ,  0    ,  0    },
             {1.0/9,  1.0/18,  1.0/36, -1.0/6, -1.0/12,  1.0/6,  1.0/12,  0    , -1.0/4},
             {1.0/9, -1.0/36, -1.0/18, -1.0/6,  1.0/6 ,  0    ,  0     ,  1.0/4,  0    },
             {1.0/9,  1.0/18,  1.0/36, -1.0/6, -1.0/12, -1.0/6, -1.0/12,  0    ,  1.0/4},
             {1.0/9, -1.0/36, -1.0/18,  0    ,  0     , -1.0/6,  1.0/6 , -1.0/4,  0    },
             {1.0/9,  1.0/18,  1.0/36,  1.0/6,  1.0/12, -1.0/6, -1.0/12,  0    , -1.0/4},
             {1.0/9, -1.0/36, -1.0/18,  1.0/6, -1.0/6 ,  0    ,  0     ,  1.0/4,  0    },
             {1.0/9,  1.0/18,  1.0/36,  1.0/6,  1.0/12,  1.0/6,  1.0/12,  0    ,  1.0/4},
             {1.0/9, -1.0/36, -1.0/18,  0    ,  0     ,  1.0/6, -1.0/6 , -1.0/4,  0    }};
// Real TransformationM[9][9]{
//                 { 1,  1,  1,  1,  1,  1,  1,  1,  1},
//                 {-4,  2, -1,  2, -1,  2, -1,  2, -1},
//                 { 4,  1, -2,  1, -2,  1, -2,  1, -2},
//                 { 0, -1, -1, -1,  0,  1,  1,  1,  0},
//                 { 0, -1,  2, -1,  0,  1, -2,  1,  0},
//                 { 0,  1,  0, -1, -1, -1,  0,  1,  1},
//                 { 0,  1,  0, -1,  2, -1,  0,  1, -2},
//                 { 0,  0,  1,  0, -1,  0,  1,  0, -1},
//                 { 0, -1,  0,  1,  0, -1,  0,  1,  0}};
void KerCollideMRTFreeSurface(ACC<Real>& fStage, const ACC<Real>& f,
                              const ACC<Real>& coordinates,
                              const ACC<int>& nodeType, const ACC<Real>& Rho,
                              const ACC<Real>& U, const ACC<Real>& V,
                              const Real* tauRef, const Real* dt,
                              const int* lattIdx,
                              const ACC<Real>& acceleration,const int* idx,
                              ACC<Real>&Tau){
#ifdef OPS_2D
    const Real iCoordinate{coordinates(0,0,0)};
    const int vt{nodeType(0, 0)};
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != (int)VertexType::ImmersedSolid && 
                              VertexType_Empty != vt &&
                              vt != (int)VertexType::Wall);
    if (!collisionRequired){
        Tau(0,0) = *tauRef;
    }
    if (collisionRequired) {
        const Real rho{Rho(0, 0)};
        const Real u{U(0, 0)};
        const Real v{V(0, 0)};
        const Real Fx{acceleration(0, 0, 0)};
        const Real Fy{acceleration(1, 0, 0)};
        Real tau = (*tauRef);
        Real Mf[lattIdx[1] + 1]{0};
        Real MfStage[lattIdx[1] + 1]{0};
        //Turbulence model
        Real fstream[lattIdx[1]+1]{0};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            fstream[xiIndex] = f(xiIndex,0,0);
        }
        //moment for f
        for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                Mf[i] += TransformationM[i][xiIndex] * f(xiIndex, 0, 0);
            }
        }

        Real LinearAbsorb{0};
        if (iCoordinate <= BlockCoor[0]+DampWidth||
            iCoordinate >= BlockCoor[1]-DampWidth) { 
            //iCoordinate >= BlockCoor[1]-DampWidth
            // Real NonLinearConst{10 * L0};       // unit 1/L
            Real StartPt{0};
            if (iCoordinate <= BlockCoor[0]+DampWidth) {
                StartPt = BlockCoor[0]+DampWidth;
            } else if (iCoordinate >= BlockCoor[1]-DampWidth) {
                StartPt = BlockCoor[1]-DampWidth;
            }
            Real DistanceScale = (iCoordinate - StartPt) / DampWidth;
            LinearAbsorb =TimeStep()*
                (-LinearConst * (exp(pow(DistanceScale, ns)) - 1) / (exp(1) - 1));
            // ops_printf("Updating the halos...%e\n", LinearConst);
        }
        tau = Smagorinsky_LES(rho, u, v, tau, fstream, lattIdx);
        Tau(0, 0) = tau;
        Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        // relaxation moment
        Real RelaxationRate[lattIdx[1]+1]{1, 
                                        0.3, 
                                        1, 
                                        1, 1.5, 
                                        1, 1.5, 
                                        dtOvertauPlusdt, 
                                        dtOvertauPlusdt};
         Real AbsorbWave[lattIdx[1]+1] {1,
                        1,
                        1,
                        1 + LinearAbsorb * 1/RelaxationRate[3],
                        1 + LinearAbsorb * 1/RelaxationRate[4],
                        1 + LinearAbsorb * 1/RelaxationRate[5],
                        1 + LinearAbsorb * 1/RelaxationRate[6],
                        1 + LinearAbsorb * 2/RelaxationRate[7],
                        1 + LinearAbsorb * 2/RelaxationRate[8]};
        //moment space
         Real EquilibriumMoment[lattIdx[1]+1]{rho, 
                                        rho * (-2 + (u * u + v * v) * 
                                        (1 + LinearAbsorb * 2/RelaxationRate[1])),
                                        rho * (1 - (u * u + v * v) * 
                                        (1 + LinearAbsorb * 2/RelaxationRate[2])),
                                        rho * u, -rho * u,
                                        rho * v, -rho * v,
                                        rho * (u * u - v * v), 
                                        rho * u * v};
         Real ForceMoment[lattIdx[1]+1]   {0,   
                                    rho * 2 * (Fx * u + Fy * v), 
                                    rho * (-2 * (Fx * u + Fy * v)), 
                                    rho * Fx, -rho * Fx, 
                                    rho * Fy, -rho * Fy, 
                                    rho * 2 * (Fx * u - Fy * v), 
                                    rho * (Fy * u + Fx * v)};
        //  Real EquilibriumMoment[lattIdx[1]+1]{rho, 
        //                                 rho * ((u * u + v * v))* 
        //                                 (1 + LinearAbsorb * 2/RelaxationRate[1]),
        //                                 0,
        //                                 rho * u, 0,
        //                                 rho * v, 0,
        //                                 rho * (u * u - v * v), 
        //                                 rho * u * v};
        //  Real ForceMoment[lattIdx[1]+1]   {0,   
        //                             rho * 2 * (Fx * u + Fy * v), 
        //                             0, 
        //                             rho * Fx, 0, 
        //                             rho * Fy, 0, 
        //                             rho * 2 * (Fx * u - Fy * v), 
        //                             rho * (Fy * u + Fx * v)};
        // MRT collision
        for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
            if (vt == (int)VertexType::Fluid ||
                vt == (int)VertexType::MDPeriodic ||
                VertexType_Interface == vt || 
                (int)VertexType::Inlet == vt ||
                (int)VertexType::OutLet == vt) {
                MfStage[i] =
                    Mf[i] - RelaxationRate[i] * 
                    (Mf[i] - AbsorbWave[i] * EquilibriumMoment[i]) +
                    (1 - RelaxationRate[i] / 2) * (*dt) * ForceMoment[i];
            } else {
                MfStage[i] =
                    Mf[i] - RelaxationRate[i] * 
                    (Mf[i] - AbsorbWave[i] * EquilibriumMoment[i]);
            }
        }
        // inverse to population space
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            Real InversefStage{0};
            for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
                InversefStage += InverseM[xiIndex][i] * MfStage[i];
            }
            fStage(xiIndex, 0, 0) = InversefStage;
#ifdef CPU
            const Real res{fStage(xiIndex, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where rho=%e u=%e v=%e at "
                    "x=%e y=%e\n",
                    res, xiIndex, rho, u, v, coordinates(0, 0, 0),
                    coordinates(1, 0, 0));
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}

void KerCalcNormalDirection( ACC<Real>& NormalDirections,
                        const ACC<Real>& VolumeFraction, 
                        const ACC<int>& nodeType,
                        const ACC<int>& geometryProperty,
                        const int* lattIdx) {
    const int vt{nodeType(0, 0)};
    NormalDirections(0, 0, 0) = 0.0;
    NormalDirections(1, 0, 0) = 0.0;
    if (VertexType_Interface == vt || VertexType_NewInterface == vt) {
        Real NormalDirection[2]{};
        int vgneighbors[3][3]{};
        Real vofneighbors[3][3]{};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            vgneighbors[cx + 1][cy + 1] = geometryProperty(cx, cy);
            vofneighbors[cx + 1][cy + 1] = VolumeFraction(cx, cy);
            if ((int)VertexType::ImmersedSolid == nodeType(cx, cy)) {
                vgneighbors[cx + 1][cy + 1] = (int)VG_ImmersedSolid;
                vofneighbors[cx + 1][cy + 1] = VolumeFraction(0, 0);
            }
        }
        CalcNormalDirection(NormalDirection, vgneighbors, vofneighbors,
                            lattIdx);
        NormalDirections(0, 0, 0) = NormalDirection[0];
        NormalDirections(1, 0, 0) = NormalDirection[1];
    }
}

//define the node is No_Fluid or No_Empty to prepare for mass exchange
void KerUpdatNodeQuality(ACC<int>& Quality,const ACC<int>& nodeType,  
                      const int* lattIdx,
                      const ACC<Real>& Porosity,
                      const ACC<Real>& NormalDirections,
                      const ACC<Real>& VOF,
                      const ACC<Real>& coordinates,
                      const ACC<int>& SegmentNo) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    const Real x{coordinates(0, 0, 0)};
    // const Real y{coordinates(1, 0, 0)};
    // const Real vof{VOF(0,0)};
    Quality(0,0) = (int)Standard;
    if (VertexType_Interface == vt) {
        int NumofFluid{0};
        int NumofEmpty{0};
        int NumofSolid{0};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor=nodeType(cx, cy);
            if ((int)VertexType::Fluid == vtneighbor) {
                NumofFluid++;
            } else if (VertexType_Empty == vtneighbor) {
                NumofEmpty++;
            } else if ((int)VertexType::ImmersedSolid == vtneighbor){
                NumofSolid++;
                // Real Intersect[]{-100, -100};
                // Real Distance{-100};
                Real NormalDirDotLatticeVel{-100};
                // Real IntersectFree[]{-100, -100};
                // Real NormalDirection[SpaceDim()]{NormalDirections(0, 0, 0),
                //                                  NormalDirections(1, 0, 0)};
                NormalDirDotLatticeVel = NormalDirections(0, 0, 0) * cx +
                                         NormalDirections(1, 0, 0) * cy;
// #ifdef CIRCLE
//                 CalcIntersectIBB(cx, cy, x, y, RotationCenter, Radii,
//                                  Intersect, Distance);
// #endif //CIRCLE                
// #ifdef SQUARE
//                 const Real P[]{x, y};
//                 CalcIntersectIBB(SegmentNo(0, 0), CS * TimeStep() * cx,
//                                  CS * TimeStep() * cy, P, SolidVertex,
//                                  SolidPlaneNormal, Polygon, NumofPol, Intersect,
//                                  Distance);
// #endif  // SQUARE
                // bool Subvof = false;
//                 CalcIntersectFree(NormalDirection, vof, cx, cy, x, y,
//                                   IntersectFree, Subvof);
                if (NormalDirDotLatticeVel > 0 + 1e-14) {
                    // if (vof > 0.5) {
                    //     const Real VectorInterF[SpaceDim()]{
                    //         IntersectFree[0] * CS * TimeStep() - x,
                    //         IntersectFree[1] * CS * TimeStep() - y,
                    //     };
                    //     const Real VectorInterS[SpaceDim()]{
                    //         Intersect[0] - x,
                    //         Intersect[1] - y,
                    //     };
                    //     const Real RelationofVector{
                    //         VectorInterF[0] * VectorInterS[0] +
                    //         VectorInterF[1] * VectorInterS[1]};
                    //     if (RelationofVector <= 0 - 1e-15) {
                    //         ops_printf(
                    //             "Error, the two intersect vectors"
                    //             "are opposite to each other..."
                    //             "IntersectFreeX =%e \n",
                    //             IntersectFree[0]);
                    //         assert(RelationofVector > 0 - 1e-15);
                    //     } else {                           
                    //         if (abs(VectorInterF[0]) < abs(VectorInterS[0]) ||
                    //             abs(VectorInterF[1]) < abs(VectorInterS[1])) {
                    //             NumofEmpty++;
                    //         }
                    //     }
                    // } else if (vof <= 0.5) {
                    //     NumofEmpty++;
                    // }
                    NumofEmpty++;
                } else if (NormalDirDotLatticeVel < 0 - 1e-14) {
                    NumofFluid++;
                } //else if (abs(NormalDirection[0]) <= 0 + 1e-14 &&
                //            abs(NormalDirection[1]) <= 0 + 1e-14) {
                //     // if (Porosity(0, 0) > 1 - 1e-15) {
                //         NumofEmpty++;
                //         NumofFluid++;
                //         ops_printf(
                //             "Warning, Normal direction is zero, nx =%e, ny=%e, "
                //             "x =%e, y =%e\n",
                //             NormalDirection[0], NormalDirection[1], x, y);
                //         break;
                //     // }
                // }
            }
   
        }
        //determine the quality of the interface cells:
        if ((0==NumofFluid)&&NumofEmpty>5) {
            Quality(0,0) = (int)No_Fluid ;
        } else if ((0==NumofEmpty)&&NumofFluid>5) {
            Quality(0,0)=(int)No_Empty;
        } else {Quality(0,0)=(int)Standard;}
        if (NumofSolid>0){
            // (x>(SolidVertex[4]+SolidVertex[6])/2-0.4&&
            //     x<(SolidVertex[12]+SolidVertex[14])/2+0.4)||
            // (x>(SolidVertex[8]+SolidVertex[10])/2-0.4&&
            //     x<(SolidVertex[0]+SolidVertex[2])/2+0.4)
            if (NumofFluid==0&&NumofEmpty>2) {
                Quality(0,0) = (int)No_Fluid ;
            } else if (NumofEmpty==0&&NumofFluid>2) {
                Quality(0,0)=(int)No_Empty;
            } else {Quality(0,0)=(int)Standard;}
        }
    }
#endif  // OPS_2D    
}

//mass exchange among interface, fluid cells
void KerMassExchange(ACC<Real>& Mass, const ACC<Real>& fStage,
                     const ACC<Real>& f,
                     const ACC<Real>& VolumeFraction, 
                     const ACC<int>& nodeType,
                     const ACC<int>& Quality, 
                     const ACC<int>& geometryProperty,
                     const int* idx, const int* jofInletInterface,
                     const int* lattIdx,
                     const ACC<Real>& porosity,
                     const ACC<Real>& NormalDirections ) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    VertexQuality vq = (VertexQuality)Quality(0, 0);
    // after interface type added into enum, revise it!
    if (VertexType_Interface == vt) {
        const int j{idx[1]};
        const Real vof {VolumeFraction(0, 0)} ;
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            Real deltam{0.0};
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            if (0 == cx && 0 == cy) {
                continue;
            }
            const int vtneighbor{nodeType(cx, cy)};
            const Real Area = 0.5 * (vof * porosity(0, 0) +
                                VolumeFraction(cx, cy) * porosity(cx, cy));
            switch (vtneighbor) {
                case VertexType_Interface: {
                    VertexQuality vqNeighbor = (VertexQuality)Quality(cx, cy);
                    // Force non-physical (or undesired) spatter and bubble
                    // to fill or empty first step:
                    if ((Standard == vq && Standard == vqNeighbor) ||
                        (No_Fluid == vq && No_Fluid == vqNeighbor) ||
                        (No_Empty == vq && No_Empty == vqNeighbor)) {
                        deltam = Area *
                                 (fStage(OPP[xiIndex], cx, cy) -fStage(xiIndex, 0, 0));
                    } else if ((Standard == vq && No_Empty == vqNeighbor) ||
                               (No_Fluid == vq && Standard == vqNeighbor) ||
                               (No_Fluid == vq && No_Empty == vqNeighbor)) {
                        deltam = Area *
                                 (-fStage(xiIndex, 0, 0));
                    } else if ((Standard == vq && No_Fluid == vqNeighbor) ||
                               (No_Empty == vq && Standard == vqNeighbor) ||
                               (No_Empty == vq && No_Fluid == vqNeighbor)) {
                        deltam = Area *
                                 (fStage(OPP[xiIndex], cx, cy));
                    }
                    break;
                }
                case (int)VertexType::Fluid: {
                    deltam =0.5 * (porosity(0, 0) + porosity(cx, cy))*
                        (fStage(OPP[xiIndex], cx, cy) - fStage(xiIndex, 0, 0));
                    break;
                }
                case VertexType_Empty: {
                    deltam = 0.0;
                    break;
                }
                case (int)VertexType::ImmersedSolid: {
                    // deltam = vof * porosity(0, 0) *
                    //     (f(OPP[xiIndex], 0, 0) -fStage(xiIndex, 0, 0));
                    // if (NormalDirections(0,0,0)*cx+NormalDirections(1,0,0)*cy>0){
                        deltam=0;
                    // }
                    break;
                } 
                case VertexType_SlipWall:{
                    const int vgneighbor = geometryProperty(cx, cy);
                    if (VG_IP == vgneighbor || VG_IM == vgneighbor) {
                        if (VertexType_Empty == nodeType(0, cy)) {
                            deltam = 0;
                        } else {
                            deltam = vof *porosity(0,0) * 
                                         (fStage(OPP[xiIndex], cx, cy) -
                                          fStage(xiIndex, 0, 0));
                        }
                    } else if (VG_JP == vgneighbor || VG_JM == vgneighbor) {
                        if (VertexType_Empty == nodeType(cx, 0)) {
                            deltam = 0;
                        } else {
                            deltam = vof *porosity(0,0) * 
                                         (fStage(OPP[xiIndex], cx, cy) -
                                          fStage(xiIndex, 0, 0));
                            // In fact, interface and fluid need to be addressed
                            // separately
                        }
                    } else {
                        deltam = 0.0;
                    }
                    break;
                }
                case (int)VertexType::Inlet: {
                    if (j + cy < *jofInletInterface) {
                        deltam = 0.5 * (porosity(0, 0) + porosity(cx, cy))*
                                 (fStage(OPP[xiIndex], cx, cy) -
                                  fStage(xiIndex, 0, 0));
                    } else if (j + cy == *jofInletInterface) {
                        deltam = vof * porosity(0, 0) *
                                 (fStage(OPP[xiIndex], cx, cy) -
                                  fStage(xiIndex, 0, 0));
                    }
                    break;
                } 
                case (int)VertexType::OutLet: {
                    // if ((int)VertexType::Fluid == nodeType(0, cy)) {
                    //     deltam = porosity(0, 0)*
                    //              (fStage(OPP[xiIndex], cx, cy) -
                    //               fStage(xiIndex, 0, 0));
                    // } else if (VertexType_Interface == nodeType(0, cy)) {
                    //     deltam = vof * porosity(0, 0) *
                    //              (fStage(OPP[xiIndex], cx, cy) -
                    //               fStage(xiIndex, 0, 0));
                    // }
                    deltam=0;
                    break;
                } 
                default: {
                    deltam = 0.0;
                    // exchange with boundary cells
                    break;
                }
            }
            Mass(0, 0) += deltam; 
        }
    }
#endif  // OPS_2D
}

void KerReconstructDF(ACC<Real>& f, const ACC<Real>& fStage,
                      const ACC<Real>& VolumeFraction, 
                      const ACC<int>& nodeType,
                      const ACC<int>& geometryProperty, 
                      const ACC<Real>& U, const ACC<Real>& V, 
                      const int* idx, const int* jofInletInterface,
                      const int* lattIdx,
                      const ACC<int>& Position,
                      const ACC<Real> &Tau,
                      const ACC<Real>& coordinates,
                     const ACC<Real>& NormalDirections) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    if (VertexType_Interface == vt) {
        const int j{idx[1]};
        Real rhoA {1.0} ;//the air pressure in free surface = 0
        const Real T{1.0};
        const int polyOrder{2};
        Real vof {VolumeFraction(0,0)};
        // reconstruct the missing DFs from empty cells
        // Warning:do not reconstruct DFs n*c>0 
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor {nodeType(cx, cy)};
            if(cx==0&&cy==0){
                continue;
            }
            if ( VertexType_Empty == vtneighbor ||
                ((int)VertexType::Inlet == vtneighbor &&
                 j + cy > *jofInletInterface) ||

                ((int)VertexType::OutLet == vtneighbor &&
                 VertexType_Empty == nodeType(0, cy))) {
                Real u = U(0, 0);
                Real v = V(0, 0);
                Real IntersectFree[]{-100, -100};
                // Real NormalDirection[SpaceDim()]{NormalDirections(0, 0, 0),
                //                                  NormalDirections(1, 0, 0)};
                // bool Subvof = false;
                // CalcIntersectFree(NormalDirection, vof, 
                //                 cx, cy, coordinates(0,0,0), coordinates(1,0,0), 
                //                 IntersectFree, Subvof);
                // if(abs(IntersectFree[1]-j)>0.5){
                    Real ny = NormalDirections(1,0,0);
                    IntersectFree[1] = (j+(abs(ny) <= 0 ? 0 : ny/abs(ny))*
                                            (vof - 0.5));
                // }
                // if (vof >= 0.5) {
                //     if ((int)VertexType::Fluid == nodeType(-cx, -cy) ||
                //     (VertexType_Interface == nodeType(-cx, -cy)&&
                //         VolumeFraction(-cx, -cy) >= 0.5 )) {
                //         Real Distance = (vof - 0.5);
                //         u = U(0, 0) * (1 + Distance) - U(-cx, -cy) * Distance;
                //         v = V(0, 0) * (1 + Distance) - V(-cx, -cy) * Distance;
                //     }
                // }else {
                //     if ((int)VertexType::Fluid == nodeType(-cx, -cy) ||
                //         (VertexType_Interface == nodeType(-cx, -cy) &&
                //         VolumeFraction(-cx, -cy) >= 0.5)) {
                //         if ((int)VertexType::Fluid == nodeType(-2*cx, -2*cy)) {
                //             Real Distance = (vof - 0.5+1);
                //             u = U(-cx, -cy) * (1 + Distance) - U(-2*cx, -2*cy) * Distance;
                //             v = V(-cx, -cy) * (1 + Distance) - V(-2*cx, -2*cy) * Distance;
                //         }else {
                //             if (vof <= 0.25){
                //                 u = U(-cx, -cy);
                //                 v = V(-cx, -cy);
                //             }
                //         }
                //     }
                // }
                if ((int)VertexType::Fluid == nodeType(-cx, -cy) ||
                    (VertexType_Interface == nodeType(-cx, -cy))) {
                    // if (vof >= 0.5) {
                    //     &&
                    //  VolumeFraction(-cx, -cy) >= 0.5
                        Real Distance = vof - 0.5 ;
                        u = U(0, 0) * (1 + Distance) - U(-cx, -cy) * Distance;
                        v = V(0, 0) * (1 + Distance) - V(-cx, -cy) * Distance;
                    // } else {
                    //     u = U(-cx, -cy);
                    //     v = V(-cx, -cy);
                    // }
                }
                rhoA = 1 + (IntersectFree[1] - (jMax+VOFofInletInterface-0.5)) * CS *
                              TimeStep() * Gravity;
                // rhoA =
                // 1+(j-0.5+(VolumeFraction(0,0)-0.5))*CS*TimeStep()*Gravity;
                // the pressure inside OWC chamber
                if ((int)VertexPosition::In_Chamber == Position(0, 0)) {
                    rhoA += DeltaRho;
                }
                // absorb wave in sponge region
                Real feq = CalcBGKFeqFA(xiIndex, rhoA, u, v, T, polyOrder,
                                Tau(0, 0), coordinates(0,0,0));
                Real feqopp =
                    CalcBGKFeqFA(OPP[xiIndex], rhoA, u, v, T, polyOrder,
                                Tau(0, 0), coordinates(0,0,0));
                f(OPP[xiIndex], 0, 0) = feq + feqopp - fStage(xiIndex, 0, 0);
            }
        }
    }
#endif  // OPS_2D
}

void KerPreUpdateFlag(ACC<int>& nodeType, const ACC<Real>& VolumeFraction,
                      const ACC<int>& Quality,
                      const int* idx, const int* jofInletInterface,
                      const int* lattIdx) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    const int vq {Quality(0, 0)};
    if (VertexType_Interface == vt) {
        const int j{idx[1]};
        const Real vof {VolumeFraction(0, 0)};
        const Real OFFSET {0.001} ;  // should be a global variable
        // Force non-physical (or undesired) spatter and bubble
        // to fill or empty second step:
        if ((vof < 0.0 - OFFSET) || ((int)No_Fluid == vq && vof < 0.1 )) {
            bool neighborisI2F{false};
            bool neighborisFluidatInlet{false};
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                int cx = (int)XI[xiIndex * LATTDIM];
                int cy = (int)XI[xiIndex * LATTDIM + 1];
                // const Real vofneighbor{VolumeFraction(cx, cy)};
                // if (vofneighbor > 1.0 + OFFSET ||
                //     ((int)No_Empty == Quality(cx, cy) && vofneighbor > 0.9)) {
                //     neighborisI2F = true;
                //     break;
                // }
                if ((int)VertexType::Inlet == nodeType(cx, cy) &&
                    j + cy < *jofInletInterface) {
                    neighborisFluidatInlet = true;
                    break;
                }
            }
            // do not change the flag of (interface next to the interface2fluid)
            // to interface2empty
            if (!(neighborisI2F||neighborisFluidatInlet)) {
                nodeType(0, 0) = VertexType_Interface2Empty;
            }
        } else if ((vof > 1.0 + OFFSET) || ((int)No_Empty == vq && vof > 0.9)) {
            bool neighborisEmptyatInlet{false};
            bool neighborisI2E{false};
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                int cx = (int)XI[xiIndex * LATTDIM];
                int cy = (int)XI[xiIndex * LATTDIM + 1];
                const Real vofneighbor{VolumeFraction(cx, cy)};
                if ((int)VertexType::Inlet == nodeType(cx, cy) &&
                    j + cy > *jofInletInterface) {
                    neighborisEmptyatInlet = true;
                    break;
                }
                if ((vofneighbor < 0.0 - OFFSET) ||
                    ((int)No_Fluid == Quality(cx, cy) && vofneighbor < 0.1 )) {
                    neighborisI2E = true;
                    break;
                }
            }
            if (!(neighborisI2E||neighborisEmptyatInlet)) {
            nodeType(0, 0) = VertexType_Interface2Fluid;
            }
        }
    }
#endif  // OPS_2D
}

void KerBuildNewInterface(ACC<int>& nodeType, 
                          const int* idx, 
                          const int* jofInletInterface,
                          const int* lattIdx) {
#ifdef OPS_2D
    const int vt{nodeType(0, 0)};
    const int j{idx[1]};
    if ((int)VertexType::Fluid == vt) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            if (VertexType_Interface2Empty == nodeType(cx, cy) ||
                ((int)VertexType::Inlet == nodeType(cx, cy) &&
                 j + cy > *jofInletInterface)) {
                nodeType(0, 0) = VertexType_Interface;
                break;
            }
        }
    } else if (VertexType_Empty == vt) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            if (VertexType_Interface2Fluid == nodeType(cx, cy) ||
                ((int)VertexType::Inlet == nodeType(cx, cy) &&
                 j + cy < *jofInletInterface)) {
                nodeType(0, 0) = VertexType_NewInterface;
                break;
            }
        }                         // local
    }
#endif  // OPS_2D
}


#ifdef FreeMpi_2D
void KerCalcTotalWeight(ACC<Real>& TotalWeightsWrong, ACC<Real>& TotalWeightsRight,
                        ACC<Real>& NormalDirections,
                        const ACC<Real>& VolumeFraction, const ACC<int>& nodeType,
                        const ACC<int>& geometryProperty, const ACC<int>& Quality, 
                        const int* idx, const int* jofInletInterface,
                        const int* lattIdx,const ACC<Real>& coordinates) {
#ifdef OPS_2D
    const int vt{nodeType(0, 0)};
    // NormalDirections(0, 0, 0) = 0.0;
    // NormalDirections(1, 0, 0) = 0.0;
    TotalWeightsWrong( 0, 0) = 0.0;
    TotalWeightsRight( 0, 0) = 0.0;
    if (VertexType_Interface2Fluid == vt || VertexType_Interface2Empty == vt) {
        const int j{idx[1]};
        Real NormalDirection[2]{};
        Real WeightforMass[lattIdx[1] + 1]{0.0} ;
        Real WeightforMassEven[lattIdx[1] + 1]{0.0} ;
        Real TotalWeight{0.0};
        int vgneighbors[3][3]{};
        Real vofneighbors[3][3]{};
        bool neighborisImmerseSolid{false};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            vgneighbors[cx + 1][cy + 1] = geometryProperty(cx, cy);
            vofneighbors[cx + 1][cy + 1] = VolumeFraction(cx, cy);
            if ((int)VertexType::ImmersedSolid == nodeType(cx, cy)) {
                vgneighbors[cx + 1][cy + 1] = (int)VG_ImmersedSolid;
                vofneighbors[cx + 1][cy + 1] = VolumeFraction(0, 0);
                neighborisImmerseSolid = true;
            }
        }
        CalcNormalDirection(NormalDirection, vgneighbors, vofneighbors,
                            lattIdx);
        NormalDirections(0, 0, 0) = NormalDirection[0];
        NormalDirections(1, 0, 0) = NormalDirection[1];
        // calculate normal direction along interface surface,
        if (VertexType_Interface2Fluid == vt) {
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                int cx = (int)XI[xiIndex * LATTDIM];
                int cy = (int)XI[xiIndex * LATTDIM + 1];
                const int vtneighbor{nodeType(cx, cy)};
                Real NormalDirDotLatticeVel =
                    NormalDirection[0] * cx + NormalDirection[1] * cy;
                if (VertexType_Interface == vtneighbor ||
                    VertexType_NewInterface == vtneighbor ||
                // the unfilled cell at Inlet can also be regarded as inerface
                    ((int)VertexType::Inlet == vtneighbor &&
                     j + cy == *jofInletInterface) ||
                // the unfilled cell at outlet might be inerface or new interface
                // as the flags at outlet are given as its neighboring values.
                    ((int)VertexType::OutLet == vtneighbor &&
                     (VertexType_Interface == nodeType(0, cy) ||
                      VertexType_NewInterface == nodeType(0, cy)))) {
                    if (NormalDirDotLatticeVel > 1e-8) {
                        WeightforMass[xiIndex] = NormalDirDotLatticeVel;
                    }
                    // WeightforMass[xiIndex] = WEIGHTS[xiIndex];
                    WeightforMassEven[xiIndex] = 1.0;
                }
                TotalWeight += WeightforMass[xiIndex];
            }
        }
        if (VertexType_Interface2Empty == vt) {
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                int cx = (int)XI[xiIndex * LATTDIM];
                int cy = (int)XI[xiIndex * LATTDIM + 1];
                const int vtneighbor{nodeType(cx, cy)};
                Real NormalDirDotLatticeVel =
                    NormalDirection[0] * cx + NormalDirection[1] * cy;
                if (VertexType_Interface == vtneighbor ||
                    VertexType_NewInterface == vtneighbor ||

                    ((int)VertexType::Inlet == vtneighbor &&
                     j + cy == *jofInletInterface) ||

                    ((int)VertexType::OutLet == vtneighbor &&
                     (VertexType_Interface == nodeType(0, cy) ||
                      VertexType_NewInterface == nodeType(0, cy)))) {
                    if (NormalDirDotLatticeVel < -1e-8) {
                        WeightforMass[xiIndex] = -NormalDirDotLatticeVel;
                    }
                    //  WeightforMass[xiIndex] = WEIGHTS[xiIndex];
                    WeightforMassEven[xiIndex] = 1.0;
                }
                TotalWeight += WeightforMass[xiIndex];
            }
        }
        // for KerAllocateMassMpi,
        // store the totalweight that is not strictly correct
        TotalWeightsWrong(0, 0) = TotalWeight;
        // EssentiallyEqual(&TotalWeight,&ZERO,EPS)
        if (TotalWeight < 1e-8) {
            TotalWeight = 0.0;
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                WeightforMass[xiIndex] = WeightforMassEven[xiIndex];
                TotalWeight += WeightforMass[xiIndex];
            }
        }
        if (TotalWeight < 1e-8 && ((int)No_Fluid != Quality(0, 0) &&
                                   (int)No_Empty != Quality(0, 0) &&
                                   (!neighborisImmerseSolid))) {
            ops_printf(
                "Error! TotalWeight =%e becomes invalid! Maybe updating flag "
                "is wrong "
                "at x=%f y=%f\n",
                TotalWeight, coordinates(0, 0, 0), coordinates(1, 0, 0));
            // Indeed the 'totalweight=0' always reflects wrong flag updating
            // assert(!(TotalWeight < 1e-8));
        }
        TotalWeightsRight(0, 0) = TotalWeight;
    }
#endif  // OPS_2D
}

void KerAllocateMassMpi(ACC<Real>& Mass, const ACC<int>& nodeType,
                        const ACC<Real>& Rho, 
                        const ACC<Real>& TotalWeightsWrong,
                        const ACC<Real>& TotalWeightsRight,
                        const ACC<Real>& NormalDirections, 
                        const int* lattIdx,
                        const ACC<Real>& VolumeFraction, 
                        const ACC<Real>& porosity,
                        const ACC<Real>& f) {
#ifdef OPS_2D
    const int vt = nodeType(0, 0);
    if (VertexType_Interface == vt || VertexType_NewInterface == vt) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor = nodeType(cx, cy);
            const Real vofneighbor = VolumeFraction(cx, cy);
            if ((VertexType_Interface2Fluid == vtneighbor && vofneighbor > 1 ) ||
                (VertexType_Interface2Empty == vtneighbor && vofneighbor < 0 )) {
                if (TotalWeightsRight(cx, cy) < 1e-8) {
                    continue;
                }
                Real NormalDirDotLatticeVel =
                    -(NormalDirections(0, cx, cy) * cx +
                      NormalDirections(1, cx, cy) * cy);
                // Note:the negative sign from that we stand at interface
                Real Weight{0.0};
                Real TotalWeight = TotalWeightsRight(cx, cy);//TotalWeights( cx, cy)
                Real ExcessMass{0.0};
                if (VertexType_Interface2Fluid == vtneighbor) {
                    if (NormalDirDotLatticeVel > 1e-8) {
                        Weight = NormalDirDotLatticeVel;
                    } else {
                        if (TotalWeightsWrong(cx, cy) < 1e-8) {
                            Weight = 1.0;
                        }
                    }
                    // Weight = WEIGHTS[xiIndex];
                    Real rho = 0;
                    for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
                        rho += f(xiIdx, cx, cy);
                    }
                    ExcessMass = Mass(cx, cy) - porosity(cx, cy) * rho;
                    // ops_printf("ExcessMass =%e, Weight =%e vofneighbor= %e\n",
                    // ExcessMass, Weight, vofneighbor);
                } else if (VertexType_Interface2Empty == vtneighbor) {
                    if (NormalDirDotLatticeVel < -1e-8){
                        Weight = -NormalDirDotLatticeVel;
                    }else{
                        if (TotalWeightsWrong(cx, cy) < 1e-8) {
                            Weight = 1.0;
                        }
                    }
                    // Weight = WEIGHTS[xiIndex];
                    ExcessMass = Mass(cx, cy);
                }
                Mass(0, 0) += ExcessMass * Weight/ TotalWeight;
            }
        }
    }
#endif  // OPS_2D
}
#endif //FreeMpi_2D

void KerCalcVolumeFraction(ACC<Real>& Mass, ACC<Real>& VolumeFraction,
                          ACC<int>& nodeType, const ACC<Real>& Rho,
                          const ACC<Real>& porosity,
                          const ACC<Real>& f,
                          const int* lattIdx) {
    const int vt{nodeType(0, 0)};
    Real rho = 0;
    for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
        rho += f(xiIdx, 0, 0);
    }
    if (VertexType_Interface == vt || 
        VertexType_NewInterface == vt) {
        nodeType(0, 0) = VertexType_Interface;
        VolumeFraction(0, 0) = Mass(0, 0) / (porosity(0, 0) * rho);
    } else if ((int)VertexType::Fluid == vt ||
               VertexType_Interface2Fluid == vt ||
               VertexType_Solid2Fluid == vt) {
        nodeType(0, 0) = (int)VertexType::Fluid;
        Mass(0, 0) = porosity(0,0)* rho;
        VolumeFraction(0, 0) = 1.0;
    } else if (VertexType_Empty == vt || 
               VertexType_Interface2Empty == vt) {
        nodeType(0, 0) = VertexType_Empty;
        Mass(0, 0) = 0.0;
        VolumeFraction(0, 0) = 0.0;
    } else {
        Mass(0, 0) = 0.0;
        VolumeFraction(0, 0) = 0.0;
    }
    assert(!(isnan(Mass(0, 0)) | isinf(Mass(0, 0))));
}

void KerHalfwayBounceBack(ACC<Real>& fStage, const ACC<int>& nodeType,
                          const int* lattIdx,const ACC<int>& geometryProperty,
                          const ACC<Real>& Rho) {
#ifdef OPS_2D
    // VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
    if ((int)VertexType::Wall == nodeType(0,0)) {
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            int cx = (int)XI[xiIdx * LATTDIM];
            int cy = (int)XI[xiIdx * LATTDIM + 1];
            const int vtneighbor{nodeType(cx, cy)};
            if ((int)VertexType::Fluid == vtneighbor ||
                VertexType_Interface == vtneighbor) {
                fStage(xiIdx, 0, 0) = fStage(OPP[xiIdx], cx, cy);
                // fStage(xiIdx, 0, 0) =
                //     CalcBGKFeq(xiIdx, Rho(cx, cy), 0, 0, 1, 2);
                // for slip boundary condition at bottom
                //  TODO should check if it needs to be revised when body move
                // if (VG_JP == vg) {
                //     for (int xiIdx2 = lattIdx[0]; xiIdx2 <= lattIdx[1];
                //          xiIdx2++) {
                //         int cxx = (int)XI[xiIdx2 * LATTDIM];
                //         int cyy = (int)XI[xiIdx2 * LATTDIM + 1];
                //         if (cxx == cx && cyy == -cy) {
                //             fStage(xiIdx, 0, 0) = fStage(xiIdx2, 0, cy);
                //             break;
                //         }
                //     }
                // }
            }
        }
    }

#endif
}

void KerUnitInterBounceBack(ACC<Real>& f, const ACC<Real>& fStage, 
                        const ACC<int>& nodeType,
                        const int* lattIdx,
                        const ACC<Real>& Rho) {
#ifdef OPS_2D
    // VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
    const int vt{nodeType(0, 0)};
    Real u{0};
    Real v{0};
    Real rho{1};
    // Real Distance{0.7};
    if ((int)VertexType::Fluid == vt || (VertexType_Interface == vt)) {
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            int cx = (int)XI[xiIdx * LATTDIM];
            int cy = (int)XI[xiIdx * LATTDIM + 1];
            const int vtneighbor{nodeType(cx, cy)};
            if ((int)VertexType::Wall == vtneighbor) {
                if ((int)VertexType::Fluid != nodeType(-cx, -cy) &&
                    VertexType_Interface != nodeType(-cx, -cy)) {
                    f(OPP[xiIdx], 0, 0) =
                        ((fStage(xiIdx, 0, 0)+ DistanceforUIBBWall * fStage(OPP[xiIdx], 0, 0)) +
                         2 * WEIGHTS[xiIdx] * rho * (-cx * u - cy * v) * CS)/(DistanceforUIBBWall + 1);
                } else {
                    // if ((int)VertexType::Fluid != nodeType(-2 * cx, -2 * cy) &&
                    //     VertexType_Interface != nodeType(-2 * cx, -2 * cy)) {
                        f(OPP[xiIdx], 0, 0) =
                            (DistanceforUIBBWall * fStage(OPP[xiIdx], 0, 0) +
                             DistanceforUIBBWall * fStage(xiIdx, 0, 0) +
                             (1 - DistanceforUIBBWall) * fStage(xiIdx, -cx, -cy) +
                             2 * WEIGHTS[xiIdx] * rho * (-cx * u - cy * v) *
                                 CS) /
                            (DistanceforUIBBWall + 1);
                    // } else {
                    //     Real f_xw =
                    //         Distance * (1 + Distance) / 2 *
                    //             fStage(xiIdx, 0, 0) +
                    //         (1 - Distance) * (1 + Distance) *
                    //             fStage(xiIdx, -cx, -cy) -
                    //         Distance * (1 - Distance) / 2 *
                    //             fStage(xiIdx, -2 * cx, -2 * cy) +
                    //         2 * WEIGHTS[xiIdx] * rho * (-cx * u - cy * v) * CS;
                    //     f(OPP[xiIdx], 0, 0) =
                    //         2 / (1 + Distance) / (2 + Distance) * f_xw +
                    //         2 * Distance / (1 + Distance) *
                    //             fStage(OPP[xiIdx], 0, 0) -
                    //         Distance / (2 + Distance) *
                    //             fStage(OPP[xiIdx], -cx, -cy);
                    // }
                }
            }
        }
    }

#endif
}

void KeIntersectOfSolid(ACC<Real>& IntersectofSolid, 
                          const ACC<Real>& coordinates,
                          const ACC<int>& nodeType, const int* lattIdx,
                          const ACC<int>& SegmentNo) {
    const int vt{nodeType(0, 0)};
    const Real x{coordinates(0, 0, 0)};
    const Real y{coordinates(1, 0, 0)};
    if ((int)VertexType::Fluid == vt || VertexType_Interface == vt ||
        VertexType_Empty == vt) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            IntersectofSolid(xiIndex* (LATTDIM+1)+1,0,0) =-100;
            IntersectofSolid(xiIndex* (LATTDIM+1)+2,0,0) =-100;
            IntersectofSolid(xiIndex* (LATTDIM+1),0,0) =-100;
            const int vtneighbor = nodeType(cx, cy);
            if ((int)VertexType::ImmersedSolid == vtneighbor) {
                Real Intersect[]{-100, -100};
                Real Distance{-100};
#ifdef CIRCLE
                CalcIntersectIBB(cx, cy, x, y, RotationCenter, Radii, Intersect,
                                 Distance);
#endif  // CIRCLE
#ifdef SQUARE

                const Real P[]{x, y};
                CalcIntersectIBB(SegmentNo(0, 0), CS * TimeStep() * cx,
                                 CS * TimeStep() * cy, P, SolidVertex,
                                 SolidPlaneNormal, Polygon, NumofPol, Intersect,
                                 Distance);
                // ops_printf(
                // "sv = %d, x =%e, y =%e, Intersectx =%e,Intersecty =%e,
                // Distance =%e, xiIndex=%d \n", SegmentNo(0, 0), P[0], P[1],
                // Intersect[0], Intersect[1], Distance, xiIndex);
#endif  // SQUARE
                assert(Distance <= 1 + 1e-10 && Distance >= -0.0 - 1e-10);
                IntersectofSolid(xiIndex* (LATTDIM+1)+1,0,0) =Intersect[0];
                IntersectofSolid(xiIndex* (LATTDIM+1)+2,0,0) =Intersect[1];
                IntersectofSolid(xiIndex* (LATTDIM+1),0,0) =Distance;
            }
        }
    }
}

void KerUIBB_GIMEM(ACC<Real>& HydroForce, ACC<Real>& f, 
                          const ACC<Real>& fStage,
                          const ACC<Real>& coordinates,
                          const ACC<int>& nodeType, const int* lattIdx,
                          const ACC<Real>& VOF, const ACC<Real>& rho,
                          const ACC<Real>& Porosity,
                          const ACC<Real>& NormalDirections,
                          const ACC<Real>& U, 
                          const ACC<Real>& V,
                          const ACC<int>& SegmentNo,
                          const ACC<int>& Position,
                          const ACC<Real>& IntersectofSolid) {
    HydroForce(0, 0, 0) =0;
    HydroForce(1, 0, 0) =0;
    HydroForce(2, 0, 0) = 0;
    const int vt{nodeType(0, 0)};
    const int vp{Position(0, 0)};
    Real vof{VOF(0, 0)};
    const Real x{coordinates(0, 0, 0)};
    const Real y{coordinates(1, 0, 0)};
    const Real CentralX{RotationCenter[0]};
    const Real CentralY{RotationCenter[1]};
    Real MeshSize{CS * TimeStep()};
    const int DIM{SpaceDim()};
    const Real T{1.0};
    const int polyOrder{2};
    if ((int)VertexType::Fluid == vt || VertexType_Interface == vt ||
         (VertexType_Empty == vt &&
          (int)VertexPosition::In_Chamber == vp)) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor = nodeType(cx, cy);
            if ((int)VertexType::ImmersedSolid == vtneighbor) {
                Real Intersect[]{IntersectofSolid(xiIndex* (LATTDIM+1)+1,0,0), 
                                 IntersectofSolid(xiIndex* (LATTDIM+1)+2,0,0)};
                Real Distance{IntersectofSolid(xiIndex* (LATTDIM+1),0,0)};
                Real NormalDirDotLatticeVel{-100};
                Real IntersectFree[]{-100, -100};
                bool ImplementFSK{false};
                Real uSolid = VelocityofRoCenter[0] -
                              VelocityofRoCenter[2] * (Intersect[1] - CentralY);
                Real vSolid = VelocityofRoCenter[1] +
                              VelocityofRoCenter[2] * (Intersect[0] - CentralX);
                if (VertexType_Interface == vt) {
                    Real NormalDirection[SpaceDim()]{NormalDirections(0, 0, 0),
                                                     NormalDirections(1, 0, 0)};
                    NormalDirDotLatticeVel = NormalDirections(0, 0, 0) * cx +
                                             NormalDirections(1, 0, 0) * cy;
                    // bool Subvof = false;
                    // CalcIntersectFree(NormalDirection, vof, cx, cy, x, y,
                    //                   IntersectFree,Subvof);
                    // if (NormalDirDotLatticeVel > 0 + 1e-14 && vof > 0.5) {
                    //     const Real VectorInterF[SpaceDim()]{
                    //         IntersectFree[0] * MeshSize - x,
                    //         IntersectFree[1] * MeshSize - y,
                    //     };
                    //     const Real VectorInterS[SpaceDim()]{
                    //         Intersect[0] - x,
                    //         Intersect[1] - y,
                    //     };
                    //     const Real RelationofVector{
                    //         VectorInterF[0] * VectorInterS[0] +
                    //         VectorInterF[1] * VectorInterS[1]};
                    //     if (RelationofVector < 0 - 1e-15) {
                    //         ops_printf(
                    //             "Error, the two intersect vectors"
                    //             "are opposite to each other..."
                    //             "x=%e,y=%e,IntersectFreeX =%e IntersectFreeY =%e"
                    //             "IntersectX =%e IntersectY =%e,"
                    //             "NormalDirection[0]=%e, NormalDirection[1]=%e\n",
                    //             x,y,IntersectFree[0],IntersectFree[1],
                    //             Intersect[0],Intersect[1],
                    //             NormalDirection[0],NormalDirection[1]);
                    //         assert(RelationofVector > 0 - 1e-15);
                    //     } else {
                    //         if (abs(VectorInterF[0]) < abs(VectorInterS[0]) ||
                    //             abs(VectorInterF[1]) < abs(VectorInterS[1])) {
                    //             ImplementFSK = true;
                    //         }
                    //     }
                    // }
                    if (ImplementFSK||
                        (NormalDirDotLatticeVel > 0 + 1e-14)) {
                        ImplementFSK = true;
                        Real u = U(0, 0);
                        Real v = V(0, 0);
                // if(abs(IntersectFree[1]-y/MeshSize)>0.5){
                    IntersectFree[1] = (y/MeshSize+(vof - 0.5));
                // }
                        if ((int)VertexType::Fluid == nodeType(-cx, -cy) ||
                            (VertexType_Interface == nodeType(-cx, -cy))) {
                            Real DistanceFree = vof - 0.5;
                            u = U(0, 0) * (1 + DistanceFree) -
                                U(-cx, -cy) * DistanceFree;
                            v = V(0, 0) * (1 + DistanceFree) -
                                V(-cx, -cy) * DistanceFree;
                        }
                        Real rhoA = 1 + (IntersectFree[1]-(jMax+VOFofInletInterface-0.5) ) * CS *
                                                 TimeStep() * Gravity;
                        // the pressure inside OWC chamber
                        if ((int)VertexPosition::In_Chamber == vp) {
                            rhoA += DeltaRho;
                        }

                        Real feq =
                            CalcBGKFeq(xiIndex, rhoA, u, v, T, polyOrder);
                        Real feqopp =
                            CalcBGKFeq(OPP[xiIndex], rhoA, u, v, T, polyOrder);
                        f(OPP[xiIndex], 0, 0) =
                            feq + feqopp - fStage(xiIndex, 0, 0);
                    }
                }
                if ((int)VertexType::Fluid == vt || (VertexType_Interface == vt &&
                      (!ImplementFSK))) {
                    
                    if ((int)VertexType::Fluid != nodeType(-cx, -cy)
                            &&VertexType_Interface != nodeType(-cx, -cy)) {
                        f(OPP[xiIndex], 0, 0) =
                                (Distance * fStage(OPP[xiIndex], 0, 0) +
                                 Distance * fStage(xiIndex, 0, 0) +
                                 (1 - Distance) * fStage(xiIndex, 0, 0) +
                                 2 * WEIGHTS[xiIndex] * rho(0, 0) *
                                     (-cx * uSolid - cy * vSolid) * CS) /
                                (Distance + 1);
                    } else {
                        // if ((int)VertexType::Fluid !=
                        //         nodeType(-2 * cx, -2 * cy) &&
                        //     VertexType_Interface !=
                        //         nodeType(-2 * cx, -2 * cy)) {
                            f(OPP[xiIndex], 0, 0) =
                                (Distance * fStage(OPP[xiIndex], 0, 0) +
                                 Distance * fStage(xiIndex, 0, 0) +
                                 (1 - Distance) * fStage(xiIndex, -cx, -cy) +
                                 2 * WEIGHTS[xiIndex] * rho(0, 0) *
                                     (-cx * uSolid - cy * vSolid) * CS) /
                                (Distance + 1);
                        // } else {
                        //     Real f_xw =Distance*(1+Distance)/2*fStage(xiIndex, 0, 0)+
                        //                 (1-Distance)*(1+Distance)*fStage(xiIndex, -cx, -cy)-
                        //                 Distance*(1-Distance)/2*fStage(xiIndex, -2*cx, -2*cy)+
                        //                 2 * WEIGHTS[xiIndex] * rho(0, 0) * (-cx * uSolid - cy * vSolid) * CS;
                        //     f(OPP[xiIndex], 0, 0) = 2/(1+Distance)/(2+Distance)*f_xw+
                        //                             2*Distance/(1+Distance)*fStage(OPP[xiIndex], 0, 0)-
                        //                             Distance/(2+Distance)*fStage(OPP[xiIndex], -cx, -cy);
                                 
                        // }
                    }
                    const Real res{f(OPP[xiIndex], 0, 0)};
                    if (isnan(res) || res <= 0 || isinf(res)) {
                        ops_printf(
                            "Error! Distribution function = %e becomes invalid "
                            "at  "
                            "the lattice %i where rho=%e usolid=%e vsolid=%e "
                            "at "
                            "x=%e y=%e\n",
                            res, OPP[xiIndex], rho, uSolid, vSolid,
                            coordinates(0, 0, 0), coordinates(1, 0, 0));
                        assert(!(isnan(res) || res <= 0 || isinf(res)));
                    }
                    // MomentumExchangeIBB
// if ((int)VertexType::Fluid == vt) {
// if(y<10){
                    Real MassPostStream = (f(OPP[xiIndex], 0, 0)) * vof;
                    Real MassPreStream = (fStage(xiIndex, 0, 0)) * vof;
                    HydroForce(0, 0, 0) +=
                        -((MassPostStream) * (-cx * CS - uSolid) -
                          (MassPreStream) * (cx * CS - uSolid)) *
                        pow(MeshSize, DIM) / TimeStep();
                    HydroForce(1, 0, 0) +=
                        -((MassPostStream) * (-cy * CS - vSolid) -
                          (MassPreStream) * (cy * CS - vSolid)) *
                        pow(MeshSize, DIM) / TimeStep();
                    HydroForce(2, 0, 0) +=
                        -pow(MeshSize, DIM) / TimeStep() *
                        ((MassPostStream * (-cy * CS - vSolid) -
                          MassPreStream * (cy * CS - vSolid)) *
                             (Intersect[0] - CentralX) -
                         (MassPostStream * (-cx * CS - uSolid) -
                          MassPreStream * (cx * CS - uSolid)) *
                             (Intersect[1] - CentralY));
// }
// }
                }
                Real rhoA = 0;
                if ((int)VertexType::Fluid == vt ||
                    (VertexType_Interface == vt && (!ImplementFSK))) {
                    rhoA = 1;
                } else if (VertexType_Empty == vt &&
                           (int)VertexPosition::In_Chamber == vp) {
                    rhoA = -DeltaRho;
                    vof = 1;
                    //due to vof for air phase is zero
                    //negetive sign and given vof for generalization with air pressure
                }
                // buoyancy force due to atmosphere pressure and chamber pressure
                Real feq =
                    CalcBGKFeq(xiIndex, rhoA, uSolid, vSolid, T, polyOrder);
                Real feqopp = CalcBGKFeq(OPP[xiIndex], rhoA, uSolid, vSolid, T,
                                         polyOrder);
                Real BuoyPostStream = -feqopp * vof;  // - feqopp
                Real BuoyPreStream = -feq * vof;      // -feq
// if ((int)VertexType::Fluid == vt||VertexType_Empty == vt) {
// if(y<10){
                HydroForce(0, 0, 0) += -(BuoyPostStream * (-cx * CS - uSolid) -
                                         BuoyPreStream * (cx * CS - uSolid)) *
                                       pow(MeshSize, DIM) / TimeStep();
                HydroForce(1, 0, 0) += -(BuoyPostStream * (-cy * CS - vSolid) -
                                         BuoyPreStream * (cy * CS - vSolid)) *
                                       pow(MeshSize, DIM) / TimeStep();
                HydroForce(2, 0, 0) += -pow(MeshSize, DIM) / TimeStep() *
                                       ((BuoyPostStream * (-cy * CS - vSolid) -
                                         BuoyPreStream * (cy * CS - vSolid)) *
                                            (Intersect[0] - CentralX) -
                                        (BuoyPostStream * (-cx * CS - uSolid) -
                                         BuoyPreStream * (cx * CS - uSolid)) *
                                            (Intersect[1] - CentralY));
// }
// }
                      
                // ******************************//
            }
        }
    }
}

void KerUIBB_GIMEMBCs(ACC<Real>& f, 
                          const ACC<Real>& fStage,
                          const ACC<Real>& coordinates,
                          const ACC<int>& nodeType, const int* lattIdx,
                          const ACC<Real>& VOF, const ACC<Real>& rho,
                          const ACC<Real>& NormalDirections,
                          const ACC<Real>& U, 
                          const ACC<Real>& V,
                          const ACC<int>& Position,
                          const ACC<Real>& IntersectofSolid) {
    const int vt{nodeType(0, 0)};
    const int vp{Position(0, 0)};
    Real vof{VOF(0, 0)};
    const Real x{coordinates(0, 0, 0)};
    const Real y{coordinates(1, 0, 0)};
    const Real CentralX{RotationCenter[0]};
    const Real CentralY{RotationCenter[1]};
    Real MeshSize{CS * TimeStep()};
    const int DIM{SpaceDim()};
    const Real T{1.0};
    const int polyOrder{2};
    if ((int)VertexType::Fluid == vt || VertexType_Interface == vt ||
         (VertexType_Empty == vt &&
          (int)VertexPosition::In_Chamber == vp)) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor = nodeType(cx, cy);
            if ((int)VertexType::ImmersedSolid == vtneighbor) {
                Real Intersect[]{IntersectofSolid(xiIndex* (LATTDIM+1)+1,0,0), 
                                 IntersectofSolid(xiIndex* (LATTDIM+1)+2,0,0)};
                Real Distance{IntersectofSolid(xiIndex* (LATTDIM+1),0,0)};
                Real NormalDirDotLatticeVel{-100};
                Real IntersectFree[]{-100, -100};
                bool ImplementFSK{false};
                Real uSolid = VelocityofRoCenter[0] -
                              VelocityofRoCenter[2] * (Intersect[1] - CentralY);
                Real vSolid = VelocityofRoCenter[1] +
                              VelocityofRoCenter[2] * (Intersect[0] - CentralX);
                if (VertexType_Interface == vt) {
                    Real NormalDirection[SpaceDim()]{NormalDirections(0, 0, 0),
                                                     NormalDirections(1, 0, 0)};
                    NormalDirDotLatticeVel = NormalDirections(0, 0, 0) * cx +
                                             NormalDirections(1, 0, 0) * cy;
                    // bool Subvof = false;
                    // CalcIntersectFree(NormalDirection, vof, cx, cy, x, y,
                    //                   IntersectFree,Subvof);
                    // if (NormalDirDotLatticeVel > 0 + 1e-14 && vof > 0.5) {
                    //     const Real VectorInterF[SpaceDim()]{
                    //         IntersectFree[0] * MeshSize - x,
                    //         IntersectFree[1] * MeshSize - y,
                    //     };
                    //     const Real VectorInterS[SpaceDim()]{
                    //         Intersect[0] - x,
                    //         Intersect[1] - y,
                    //     };
                    //     const Real RelationofVector{
                    //         VectorInterF[0] * VectorInterS[0] +
                    //         VectorInterF[1] * VectorInterS[1]};
                    //     if (RelationofVector < 0 - 1e-15) {
                    //         ops_printf(
                    //             "Error, the two intersect vectors"
                    //             "are opposite to each other..."
                    //             "x=%e,y=%e,IntersectFreeX =%e IntersectFreeY =%e"
                    //             "IntersectX =%e IntersectY =%e,"
                    //             "NormalDirection[0]=%e, NormalDirection[1]=%e\n",
                    //             x,y,IntersectFree[0],IntersectFree[1],
                    //             Intersect[0],Intersect[1],
                    //             NormalDirection[0],NormalDirection[1]);
                    //         assert(RelationofVector > 0 - 1e-15);
                    //     } else {
                    //         if (abs(VectorInterF[0]) < abs(VectorInterS[0]) ||
                    //             abs(VectorInterF[1]) < abs(VectorInterS[1])) {
                    //             ImplementFSK = true;
                    //         }
                    //     }
                    // }
                    if (ImplementFSK||
                        (NormalDirDotLatticeVel > 0 + 1e-14)) {
                        ImplementFSK = true;
                        Real u = U(0, 0);
                        Real v = V(0, 0);
                // if(abs(IntersectFree[1]-y/MeshSize)>0.5){
                    IntersectFree[1] = (y/MeshSize+(vof - 0.5));
                // }
                        if ((int)VertexType::Fluid == nodeType(-cx, -cy) ||
                            (VertexType_Interface == nodeType(-cx, -cy))) {
                            Real DistanceFree = vof - 0.5;
                            u = U(0, 0) * (1 + DistanceFree) -
                                U(-cx, -cy) * DistanceFree;
                            v = V(0, 0) * (1 + DistanceFree) -
                                V(-cx, -cy) * DistanceFree;
                        }
                        Real rhoA = 1 + (IntersectFree[1]-(jMax+VOFofInletInterface-0.5) ) * CS *
                                                 TimeStep() * Gravity;
                        // the pressure inside OWC chamber
                        if ((int)VertexPosition::In_Chamber == vp) {
                            rhoA += DeltaRho;
                        }

                        Real feq =
                            CalcBGKFeq(xiIndex, rhoA, u, v, T, polyOrder);
                        Real feqopp =
                            CalcBGKFeq(OPP[xiIndex], rhoA, u, v, T, polyOrder);
                        f(OPP[xiIndex], 0, 0) =
                            feq + feqopp - fStage(xiIndex, 0, 0);
                    }
                }
                if ((int)VertexType::Fluid == vt || (VertexType_Interface == vt &&
                      (!ImplementFSK))) {
                    
                    if ((int)VertexType::Fluid != nodeType(-cx, -cy)
                            &&VertexType_Interface != nodeType(-cx, -cy)) {
                        f(OPP[xiIndex], 0, 0) =
                                (Distance * fStage(OPP[xiIndex], 0, 0) +
                                 Distance * fStage(xiIndex, 0, 0) +
                                 (1 - Distance) * fStage(xiIndex, 0, 0) +
                                 2 * WEIGHTS[xiIndex] * rho(0, 0) *
                                     (-cx * uSolid - cy * vSolid) * CS) /
                                (Distance + 1);
                    } else {
                        // if ((int)VertexType::Fluid !=
                        //         nodeType(-2 * cx, -2 * cy) &&
                        //     VertexType_Interface !=
                        //         nodeType(-2 * cx, -2 * cy)) {
                            f(OPP[xiIndex], 0, 0) =
                                (Distance * fStage(OPP[xiIndex], 0, 0) +
                                 Distance * fStage(xiIndex, 0, 0) +
                                 (1 - Distance) * fStage(xiIndex, -cx, -cy) +
                                 2 * WEIGHTS[xiIndex] * rho(0, 0) *
                                     (-cx * uSolid - cy * vSolid) * CS) /
                                (Distance + 1);
                        // } else {
                        //     Real f_xw =Distance*(1+Distance)/2*fStage(xiIndex, 0, 0)+
                        //                 (1-Distance)*(1+Distance)*fStage(xiIndex, -cx, -cy)-
                        //                 Distance*(1-Distance)/2*fStage(xiIndex, -2*cx, -2*cy)+
                        //                 2 * WEIGHTS[xiIndex] * rho(0, 0) * (-cx * uSolid - cy * vSolid) * CS;
                        //     f(OPP[xiIndex], 0, 0) = 2/(1+Distance)/(2+Distance)*f_xw+
                        //                             2*Distance/(1+Distance)*fStage(OPP[xiIndex], 0, 0)-
                        //                             Distance/(2+Distance)*fStage(OPP[xiIndex], -cx, -cy);
                                 
                        // }
                    }
                    const Real res{f(OPP[xiIndex], 0, 0)};
                    if (isnan(res) || res <= 0 || isinf(res)) {
                        ops_printf(
                            "Error! Distribution function = %e becomes invalid "
                            "at  "
                            "the lattice %i where rho=%e usolid=%e vsolid=%e "
                            "at "
                            "x=%e y=%e\n",
                            res, OPP[xiIndex], rho, uSolid, vSolid,
                            coordinates(0, 0, 0), coordinates(1, 0, 0));
                        assert(!(isnan(res) || res <= 0 || isinf(res)));
                    }
                }
            }
        }
    }
}

void KerUIBB_GIMEMFluidForce(ACC<Real>& HydroForce, const ACC<Real>& f, 
                          const ACC<Real>& fStage,
                          const ACC<int>& nodeType, const int* lattIdx,
                          const ACC<Real>& VOF,
                          const ACC<Real>& NormalDirections,
                          const ACC<int>& Position,
                          const ACC<Real>& IntersectofSolid) {
    HydroForce(0, 0, 0) =0;
    HydroForce(1, 0, 0) =0;
    HydroForce(2, 0, 0) = 0;
    const int vt{nodeType(0, 0)};
    const int vp{Position(0, 0)};
    Real vof{VOF(0, 0)};
    const Real CentralX{RotationCenter[0]};
    const Real CentralY{RotationCenter[1]};
    Real MeshSize{CS * TimeStep()};
    const int DIM{SpaceDim()};
    const Real T{1.0};
    const int polyOrder{2};
    if ((int)VertexType::Fluid == vt || VertexType_Interface == vt ||
         (VertexType_Empty == vt &&
          (int)VertexPosition::In_Chamber == vp)) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor = nodeType(cx, cy);
            if ((int)VertexType::ImmersedSolid == vtneighbor) {
                Real Intersect[]{IntersectofSolid(xiIndex* (LATTDIM+1)+1,0,0), 
                                 IntersectofSolid(xiIndex* (LATTDIM+1)+2,0,0)};
                Real Distance{IntersectofSolid(xiIndex* (LATTDIM+1),0,0)};
                Real NormalDirDotLatticeVel{-100};
                Real IntersectFree[]{-100, -100};
                bool ImplementFSK{false};
                Real uSolid = VelocityofRoCenter[0] -
                              VelocityofRoCenter[2] * (Intersect[1] - CentralY);
                Real vSolid = VelocityofRoCenter[1] +
                              VelocityofRoCenter[2] * (Intersect[0] - CentralX);
                if (VertexType_Interface == vt) {
                    NormalDirDotLatticeVel = NormalDirections(0, 0, 0) * cx +
                                             NormalDirections(1, 0, 0) * cy;
                    if (ImplementFSK||
                        (NormalDirDotLatticeVel > 0 + 1e-14)) {
                        ImplementFSK = true;
                    }
                }
                if ((int)VertexType::Fluid == vt || (VertexType_Interface == vt &&
                      (!ImplementFSK))) {
                    // MomentumExchangeIBB
// if ((int)VertexType::Fluid == vt) {
// if(y<10){
                    Real MassPostStream = (f(OPP[xiIndex], 0, 0)) * vof;
                    Real MassPreStream = (fStage(xiIndex, 0, 0)) * vof;
                    HydroForce(0, 0, 0) +=
                        -((MassPostStream) * (-cx * CS - uSolid) -
                          (MassPreStream) * (cx * CS - uSolid)) *
                        pow(MeshSize, DIM) / TimeStep();
                    HydroForce(1, 0, 0) +=
                        -((MassPostStream) * (-cy * CS - vSolid) -
                          (MassPreStream) * (cy * CS - vSolid)) *
                        pow(MeshSize, DIM) / TimeStep();
                    HydroForce(2, 0, 0) +=
                        -pow(MeshSize, DIM) / TimeStep() *
                        ((MassPostStream * (-cy * CS - vSolid) -
                          MassPreStream * (cy * CS - vSolid)) *
                             (Intersect[0] - CentralX) -
                         (MassPostStream * (-cx * CS - uSolid) -
                          MassPreStream * (cx * CS - uSolid)) *
                             (Intersect[1] - CentralY));
// }
// }
                }
                Real rhoA = 0;
                if ((int)VertexType::Fluid == vt ||
                    (VertexType_Interface == vt && (!ImplementFSK))) {
                    rhoA = 1;
                }
                // buoyancy force due to atmosphere pressure and chamber pressure
                Real feq =
                    CalcBGKFeq(xiIndex, rhoA, uSolid, vSolid, T, polyOrder);
                Real feqopp = CalcBGKFeq(OPP[xiIndex], rhoA, uSolid, vSolid, T,
                                         polyOrder);
                Real BuoyPostStream = -feqopp * vof;  // - feqopp
                Real BuoyPreStream = -feq * vof;      // -feq
// if ((int)VertexType::Fluid == vt||VertexType_Empty == vt) {
// if(y<10){
                HydroForce(0, 0, 0) += -(BuoyPostStream * (-cx * CS - uSolid) -
                                         BuoyPreStream * (cx * CS - uSolid)) *
                                       pow(MeshSize, DIM) / TimeStep();
                HydroForce(1, 0, 0) += -(BuoyPostStream * (-cy * CS - vSolid) -
                                         BuoyPreStream * (cy * CS - vSolid)) *
                                       pow(MeshSize, DIM) / TimeStep();
                HydroForce(2, 0, 0) += -pow(MeshSize, DIM) / TimeStep() *
                                       ((BuoyPostStream * (-cy * CS - vSolid) -
                                         BuoyPreStream * (cy * CS - vSolid)) *
                                            (Intersect[0] - CentralX) -
                                        (BuoyPostStream * (-cx * CS - uSolid) -
                                         BuoyPreStream * (cx * CS - uSolid)) *
                                            (Intersect[1] - CentralY));
// }
// }
                      
                // ******************************//
            }
        }
    }
}


void KerCalcBuoyForce(ACC<Real>& HydroForce, 
                          const ACC<Real>& coordinates,
                          const ACC<int>& nodeType, 
                          const int* lattIdx,
                          const ACC<Real>& Porosity,
                          const ACC<Real>& VOF,
                          const ACC<int>& SegmentNo,
                          const ACC<Real>& NormalDirections) {
    const int vt{nodeType(0, 0)};
    const int sn{SegmentNo(0, 0)};
    Real vof{VOF(0, 0)};
    const Real P[]{coordinates(0, 0, 0),coordinates(1, 0, 0)};
    const Real CentralX{RotationCenter[0]};
    const Real CentralY{RotationCenter[1]};
    Real MeshSize{CS * TimeStep()};
    const Real T{1.0};
    const int polyOrder{2};
    int DIM{SpaceDim()};
    bool BuoyancyRequired{false};
    Real NormalDirection[SpaceDim()]{NormalDirections(0, 0, 0),
                                     NormalDirections(1, 0, 0)};
    if ((int)VertexType::Fluid == vt || VertexType_Interface == vt) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor{nodeType(cx, cy)};
            if ((int)VertexType::ImmersedSolid == vtneighbor) {
                BuoyancyRequired = true;
            }
        }
        if (VertexType_Interface == vt&&BuoyancyRequired) {
            int numofempty{0};
            int numoffluid{0};
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                int cx = (int)XI[xiIndex * LATTDIM];
                int cy = (int)XI[xiIndex * LATTDIM + 1];
                const int vtneighbor{nodeType(cx, cy)};
                if (VertexType_Empty == vtneighbor) {
                    numofempty++;
                } else if ((int)VertexType::Fluid == vtneighbor) {
                    numoffluid++;
                }
            }
            if(numoffluid < numofempty){
                BuoyancyRequired=false;
            }
        }
    }
    if ((int)VertexType::ImmersedSolid == vt && Porosity(0,0)>0) {
        int numofempty{0};
        int numoffluid{0};
        int numofinterface{0};
        Real AverVof{0};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor{nodeType(cx, cy)};
            if (VertexType_Empty == vtneighbor) {
                numofempty++;
            } else if ((int)VertexType::Fluid == vtneighbor) {
                numoffluid++;
            } else if (VertexType_Interface == vtneighbor) {
                numofinterface++;
                AverVof += VOF(cx,cy);
            }
        }
        if (numofinterface==0&&numoffluid>0){
            vof =1;
        }else{
            vof = AverVof / numofinterface;
        }
        if(numoffluid >= numofempty){
            // = cause all solid nodes are counted
            BuoyancyRequired=true;
            int NumofNeighLiquid{0};
            // NormalDirection[0] = 0;
            // NormalDirection[1] = 0;
            // if(VertexType_Interface == nodeType(-1,0)||
            //     (int)VertexType::Fluid == nodeType(-1, 0)){
            //     NumofNeighLiquid++;
            //     NormalDirection[0] += NormalDirections(0, -1, 0);
            //     NormalDirection[1] += NormalDirections(1, -1, 0);
            // }else if(VertexType_Interface == nodeType(1, 0) ||
            //         (int)VertexType::Fluid == nodeType(1, 0)){
            //     NumofNeighLiquid++;
            //     NormalDirection[0] += NormalDirections(0, 1, 0);
            //     NormalDirection[1] += NormalDirections(1, 1, 0);
            // }
            // if (NumofNeighLiquid != 0) {
            //     NormalDirection[0]/=NumofNeighLiquid;
            //     NormalDirection[1]/=NumofNeighLiquid;
            // }
        }    
    }
    if (BuoyancyRequired) {
        // Real dfromcenter{sqrt((x - CentralX) * (x - CentralX) +
        //                       (y - CentralY) * (y - CentralY)) -
        //                  Radii};
        if (sn >=-1) {
            // ops_printf("starting buoyForce\n");
            // determine the cell where boundary locates
            int subgrids{1};
            if(-1==sn){
                subgrids =25;
            }else {
                subgrids =25;
            }
            const Real subMeshSize{MeshSize / subgrids};
            int subIdx[subgrids]{0};
            for (int num = 0; num < subgrids; num++) {
                subIdx[num] = num - (subgrids - (subgrids + 1) / 2);
            }
            for (int numx = 0; numx < subgrids; numx++) {
                for (int numy = 0; numy < subgrids; numy++) {
                    Real subP[]{P[0] + subMeshSize * subIdx[numx],
                                P[1] + subMeshSize * subIdx[numy]};
                    // ensure to count all subnodes outside solid
                    int segmentno = sn;
                    int inpolygon{2};
#ifdef CIRCLE
                    segmentno = sn,
                    inpolygon = InPolygon(segmentno, subP, RotationCenter, Radii);
#endif  // CIRCLE
#ifdef SQUARE
                    segmentno = sn;
                    Real d{0};
                    inpolygon = InPolygon(segmentno, subP, SolidVertex,
                                          SolidPlaneNormal, Polygon, NumofPol,d,
                                          subgrids);
#endif  // SQUARE
                    if (inpolygon>0) {
                        Real INTERSECT[3][9]{0};
                        int INPOLYGON[9]{0};
                        Real INTERSECTforRHOA{0};
                        int INPOLYGONforRHOA{0};
                        // ops_printf(" sn2 = %d, cx=%e, cy=%e\n",inpolygon,subP[0],subP[1]);
                        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1];
                             xiIndex++) {
                            int cx = (int)XI[xiIndex * LATTDIM];
                            int cy = (int)XI[xiIndex * LATTDIM + 1];
                            Real neighsubP[]{subP[0] + cx * subMeshSize,
                                             subP[1] + cy * subMeshSize};
                            Real Intersect[]{-100, -100};
                            Real Distance{-100};
                            inpolygon = 2;
#ifdef CIRCLE
                            segmentno = sn;
                            inpolygon =
                                InPolygon(segmentno, neighsubP, RotationCenter, Radii);
#endif  // CIRCLE
#ifdef SQUARE
                            segmentno = sn;
                            Real d{0};
                            inpolygon =
                                InPolygon(segmentno, neighsubP, SolidVertex,
                                          SolidPlaneNormal, Polygon, NumofPol,d,
                                          subgrids);
#endif  // SQUARE
                            if (inpolygon<0) {
                                // ops_printf(" sn3 = %d, npx=%e, npy=%e\n",sn,
                                    // neighsubP[0],neighsubP[1]);
                                // neighbor of subP is inside solid
#ifdef CIRCLE
                                CalcIntersectIBB(subMeshSize * cx,
                                                 subMeshSize * cy, subP[0],
                                                 subP[1], RotationCenter,
                                                 Radii, Intersect, Distance);
#endif  // CIRCLE
#ifdef SQUARE
                                segmentno = sn;
                                CalcIntersectIBB(
                                    segmentno, subMeshSize * cx, subMeshSize * cy,
                                    subP, SolidVertex, SolidPlaneNormal,
                                    Polygon, NumofPol, Intersect, Distance);
#endif  // SQUARE
        // ops_printf(
        //     "intersect! %e, %d, %d, %2.16f,%2.16f,%2.16f,%2.16f\n", Distance, cx, cy,
        //     P[0], P[1], Intersect[0],
        //     Intersect[1]);
                                assert(Distance <= 1+1e-10&&
                                       Distance >= -0.0-1e-10);
                                INPOLYGON[xiIndex] = inpolygon;
                                INTERSECT[0][xiIndex] = Intersect[0];
                                INTERSECT[1][xiIndex] = Intersect[1];
                                INTERSECT[2][xiIndex] = Distance;
                                INTERSECTforRHOA += Intersect[1];
                                INPOLYGONforRHOA ++;
                            }
                        }
                        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1];
                             xiIndex++) {
                                int cx = (int)XI[xiIndex * LATTDIM];
                                int cy = (int)XI[xiIndex * LATTDIM + 1];
                                if (INPOLYGON[xiIndex] < 0 && INPOLYGONforRHOA > 0) {
                                    Real Intersect[]{INTERSECT[0][xiIndex], INTERSECT[1][xiIndex]};
                                    Real Distance{INTERSECT[2][xiIndex]};
                                // if (vof >0 && vof<1){
                                //     Real IntersectFree[]{-100, -100};
                                //     bool Subvof = true;
                                //     CalcIntersectFree(NormalDirection, vof, cx, cy, 
                                //                       subP[0], subP[1],
                                //                       IntersectFree, Subvof);
                                //     if (Subvof){
                                //         vof = 1;
                                //     }else {
                                //         vof = 0;
                                //     }
                                // }
                                 // Real NormalDirection[SpaceDim()]{NormalDirections(0, 0, 0),
                    //                                  NormalDirections(1, 0, 0)};
                                // buoyancy force due to gravity
                                Real rhoA =
                                    (INTERSECTforRHOA/INPOLYGONforRHOA - 
                                        (jMax+VOFofInletInterface-0.5) * MeshSize) *
                                    Gravity;
                                Real feq = CalcBGKFeq(xiIndex, rhoA, 0, 0, T,
                                                      polyOrder);
                                Real feqopp = CalcBGKFeq(OPP[xiIndex], rhoA, 0,
                                                         0, T, polyOrder);
                                Real BuoyPostStream = -feqopp*vof;  // - feqopp
                                Real BuoyPreStream = -feq*vof;      // -feq
                                HydroForce(0, 0, 0) +=
                                    -(BuoyPostStream * (-cx * CS) -
                                      BuoyPreStream * (cx * CS)) *
                                    pow(subMeshSize, DIM) /
                                    (TimeStep() / subgrids);
                                HydroForce(1, 0, 0) +=
                                    -(BuoyPostStream * (-cy * CS) -
                                      BuoyPreStream * (cy * CS)) *
                                    pow(subMeshSize, DIM) /
                                    (TimeStep() / subgrids);
                                HydroForce(2, 0, 0) +=
                                    -pow(subMeshSize, DIM) /
                                    (TimeStep() / subgrids) *
                                    ((BuoyPostStream * (-cy * CS) -
                                      BuoyPreStream * (cy * CS)) *
                                         (Intersect[0] - CentralX) -
                                    (BuoyPostStream * (-cx * CS ) -
                                       BuoyPreStream * (cx * CS )) *
                                         (Intersect[1] - CentralY));
                            }
                        }
                    }
                }
            }
        }
    }
}

void KerTotalHydroForce(const ACC<int>& nodeType, 
                        const ACC<Real>& HydroForce,
                        const ACC<Real>& AirForce,
                        double* Fx, double* Fy, double* Mz,
                        double* FxAir, double* FyAir, double* MzAir) {                       
        *Fx = *Fx + HydroForce(0, 0, 0);
        *Fy = *Fy + HydroForce(1, 0, 0);
        *Mz = *Mz + HydroForce(2, 0, 0);
        *FxAir = *FxAir + AirForce(0, 0, 0);
        *FyAir = *FyAir + AirForce(1, 0, 0);
        *MzAir = *MzAir + AirForce(2, 0, 0);
}

void KerSetSolidVOF(ACC<Real>& Porosity, 
                          const ACC<Real>& coordinates,
                          ACC<int>& Position,
                          ACC<int>& SegmentNo) {
    const Real P[]{coordinates(0, 0, 0), coordinates(1, 0, 0)};
    Porosity(0, 0) = 1;//reinitialize
    Position(0, 0) = (int)VertexPosition::Out_Chamber;
    SegmentNo(0, 0) = -10;
    // const int dim{SpaceDim()};
#ifdef CIRCLE
    //  Initialise the Porosity here according to x,y
    // (x-0.5)^2+(y-0.75)^2=0.15^2
    // 64 96 for 128 nodes but it is not the (0.5,0.75)
    int vp{InPolygon(SegmentNo(0, 0), P, RotationCenter, Radii)};
    Position(0, 0) = vp;
    Porosity(0, 0) = 1 - CalcSolidVOF(RotationCenter, Radii, P);
#endif //CIRCLE
#ifdef SQUARE
    Real d{0};
    CalcSegmentNo(SegmentNo(0, 0), P, SolidVertex, SolidPlaneNormal, Polygon,
                      NumofPol);
    const int segmentno{SegmentNo(0, 0)};
    int vp{InPolygon(segmentno, P, SolidVertex, SolidPlaneNormal, Polygon,
                      NumofPol,d)};
    Position(0, 0) = vp;
    if (vp < 0 && -10 == SegmentNo(0, 0)) {
        Porosity(0, 0) = 0;
    }
    const int sn{SegmentNo(0, 0)};
    if (sn >= 0) {
        Real polygonVertex[]{SolidVertex[sn * SpaceDim()],
                             SolidVertex[sn * SpaceDim() + 1]};
        Real PlaneNormal[]{SolidPlaneNormal[sn * SpaceDim()],
                           SolidPlaneNormal[sn * SpaceDim() + 1]};
        if (sn>=Polygon){
            PlaneNormal[0] = -PlaneNormal[0];
            PlaneNormal[1] = -PlaneNormal[1];
        }
        // ops_printf("sn = %d, nx=%e, ny=%e\n",sn,PlaneNormal[0],PlaneNormal[1]);
        Porosity(0, 0) = 1 - CalcSolidVOF(polygonVertex, PlaneNormal, P);
    } else if (sn == -1) {
        const int subgrids{5};
        const Real subMeshSize{CS * TimeStep() / subgrids};
        int subIdx[subgrids]{0};
        for (int num = 0; num < subgrids; num++) {
            subIdx[num] = num - (subgrids - (subgrids + 1) / 2);
        }
        for (int numx = 0; numx < subgrids; numx++) {
            for (int numy = 0; numy < subgrids; numy++) {
                Real subP[]{P[0] + subMeshSize * subIdx[numx],
                            P[1] + subMeshSize * subIdx[numy]};
                Real d = 100;
                int vp{InPolygon(sn, subP, SolidVertex, SolidPlaneNormal,
                                 Polygon, NumofPol, d, subgrids)};
                if (vp < 0) {
                    // ops_printf("sn = %d, DisfromSeg1 =%e\n",sn,d);
                    Porosity(0, 0) =
                        Porosity(0, 0) -
                        (0.5 + 0.5 * d) * pow(1.0 / subgrids, SpaceDim());
                }
                if (vp > 0) {
                    // ops_printf("sn = %d, DisfromSeg2 =%e\n",sn,d);
                    Porosity(0, 0) =
                        Porosity(0, 0) -
                        (0.5 - 0.5 * d) * pow(1.0 / subgrids, SpaceDim());
                }  // double weight can further improve the accuracy of porosity
            }
        }
    }
#endif  // SQUARE
}


// define the node quality is No_Fluid or No_Empty 
// for node converted from solid
void KerUpdatevqSolid(ACC<int>& Quality,const ACC<int>& nodeType,
                      const int* lattIdx) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    Quality(0, 0) = (int)Standard;
    if (((int)VertexType::ImmersedSolid == vt ||
        VertexType_Solid2Interface == vt)) {
        int NumofFluid{0};
        int NumofEmpty{0};
        int NumofInterface{0};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor=nodeType(cx, cy);
            if ((int)VertexType::Fluid == vtneighbor) {
                NumofFluid++;
            } else if (VertexType_Empty == vtneighbor) {
                NumofEmpty++;
            } else if (VertexType_Interface == vtneighbor) {
                NumofInterface++;
            }
        }
        if ((0 == NumofFluid) && NumofEmpty > 0) {
            Quality(0, 0) = (int)No_Fluid;
        } else if (0 == NumofEmpty && NumofFluid > 0) {
            Quality(0, 0) = (int)No_Empty;
        } else if (0 == NumofFluid && 0 == NumofEmpty && 
                   0 == NumofInterface) {
            Quality(0, 0) = (int)Solid;
        } else {
            Quality(0, 0) = (int)Standard;
        }
    }
#endif  // OPS_2D    
}

void KerPreUpdateFlagIBB(ACC<int>& nodeType,
                            const ACC<int>& Quality,
                            const int* lattIdx,
                            const ACC<Real>& coordinates,
                            const ACC<int>& Position) {
    const int vt{nodeType(0,0)};
    const int vq{Quality(0,0)};
    const int vp{Position(0,0)};
    if ((int)VertexPosition::ImmerseSolid == vp) {        
         nodeType(0, 0) = (int)VertexType::ImmersedSolid;
    } else {
        if ((int)VertexType::ImmersedSolid == vt) {
            if ((int)No_Fluid == vq) {
                bool NeighborisNoempty{false};
                // for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1];
                //      xiIndex++) {
                //     int cx = (int)XI[xiIndex * LATTDIM];
                //     int cy = (int)XI[xiIndex * LATTDIM + 1];
                //     const int vqneighbor = Quality(cx, cy);
                //     if ((int)No_Empty == vqneighbor) {
                //         NeighborisNoempty = true;
                //         break;
                //     }
                // }
                if (NeighborisNoempty){
                    nodeType(0, 0) = VertexType_NewInterface;
                }else{
                    nodeType(0, 0) = VertexType_Empty;
                }
            } else if ((int)No_Empty == vq) {
                bool NeighborisNofluid{false};
                for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1];
                     xiIndex++) {
                    int cx = (int)XI[xiIndex * LATTDIM];
                    int cy = (int)XI[xiIndex * LATTDIM + 1];
                    const int vqneighbor = Quality(cx, cy);
                    if ((int)No_Fluid == vqneighbor) {
                        NeighborisNofluid = true;
                        break;
                    }
                }
                if (NeighborisNofluid){
                    nodeType(0, 0) = VertexType_NewInterface;
                }else{
                    nodeType(0,0) = VertexType_Solid2Fluid;
                }
            } else if ((int)Standard == vq) {
                nodeType(0, 0) = VertexType_NewInterface;
            } else if ((int)Solid == vq) {
                // need further address the cell
                // due to too large solid motion
                MultiLayerConvertedSolid = true;
                nodeType(0, 0) = VertexType_Solid2Interface;
            }
        }
    }  
}


void KerInitNewInterface( ACC<Real>& f,
                          ACC<Real>& fStage,
                          ACC<Real>& Rho,
                          ACC<Real>& U, ACC<Real>& V,
                          const ACC<int>& nodeType,
                          const int* lattIdx,
                          const Real* tauRef,
                          const ACC<Real>& coordinates,
                          const ACC<Real>& acceleration,
                          const ACC<Real>& vof,
                          const int* idx) {
#ifdef OPS_2D
    const int vt{nodeType(0, 0)};
    if (VertexType_Solid2Fluid == vt ||
        VertexType_NewInterface == vt) {
        Real TotalWeight{0.0};
        Real AverRho{0.0};
        Real AverU{0.0};
        Real AverV{0.0};
        const Real T{1};
        const int polyOrder{2};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor{nodeType(cx, cy)};
            if ((int)VertexType::Fluid == vtneighbor ||
                VertexType_Interface == vtneighbor ||
                VertexType_Interface2Fluid == vtneighbor) {
                Real rho{Rho(cx, cy)};
                Real u{U(cx, cy)};
                Real v{V(cx, cy)};
                // if use consistent initialize the following code can be
                // omitted
                AverRho += WEIGHTS[xiIndex] * rho;
                AverU += WEIGHTS[xiIndex] * u;
                AverV += WEIGHTS[xiIndex] * v;
                TotalWeight += WEIGHTS[xiIndex];
            }
        }
        AverRho /= TotalWeight;
        // AverRho =1 + (idx[1]+vof(0,0)-0.5 - 
        //                 (jMax+VOFofInletInterface-0.5) ) * CS * TimeStep() * Gravity;
        AverU /= TotalWeight;
        AverV /= TotalWeight;
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            fStage(xiIdx, 0, 0) =
                CalcBGKFeqFA(xiIdx, AverRho, AverU, AverV, T, polyOrder,
                             *tauRef, coordinates(0, 0, 0));
        }
        AverRho = 0;
        AverU = 0;
        AverV = 0;
        TotalWeight = 0;
        Real fneq[9] {};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor{nodeType(cx, cy)};
            if ((int)VertexType::Fluid == vtneighbor ||
                VertexType_Interface == vtneighbor ||
                VertexType_Interface2Fluid == vtneighbor) {
                Real rho{Rho(cx, cy)};
                Real u{U(cx, cy)};
                Real v{V(cx, cy)};
                //if use consistent initialize the following code can be omitted
                if (!ConsistentInitialize) {
                    rho = 0;
                    u = 0;
                    v = 0;
                    for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
                        Real cxx = CS * XI[xiIdx * LATTDIM];
                        Real cyy = CS * XI[xiIdx * LATTDIM + 1];
                        rho += f(xiIdx, cx, cy);
                        u += cxx * f(xiIdx, cx, cy);
                        v += cyy * f(xiIdx, cx, cy);
                    }
                    u = u / rho + (TimeStep() * acceleration(0, cx, cy) / 2);
                    v = v / rho + (TimeStep() * acceleration(1, cx, cy) / 2);
                    // for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
                    //     Real cxx = CS * XI[xiIdx * LATTDIM];
                    //     Real cyy = CS * XI[xiIdx * LATTDIM + 1];
                    //     fneq[xiIdx] += WEIGHTS[xiIndex] * (f(xiIdx, cx, cy)-
                    //             CalcBGKFeqFA(xiIdx, rho, u, v, T, polyOrder,
                    //          *tauRef, coordinates(0, 0, 0));
                    // } 

                }
                AverRho += WEIGHTS[xiIndex] * rho;
                // Real StaticPressure = -(j+cy-jMax)*CS*TimeStep()*Gravity;
                // AverRho += WEIGHTS[xiIndex]*(Rho(cx, cy)-StaticPressure);
                AverU += WEIGHTS[xiIndex] * u;
                AverV += WEIGHTS[xiIndex] * v;
                TotalWeight += WEIGHTS[xiIndex];
            }
        }
        AverRho /= TotalWeight;
        // AverRho =1 + (idx[1]+vof(0,0)-0.5 - 
        //                 (jMax+VOFofInletInterface-0.5) ) * CS * TimeStep() * Gravity;
        AverU /= TotalWeight;
        AverV /= TotalWeight;
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            f(xiIdx, 0, 0) =
                CalcBGKFeqFA(xiIdx, AverRho, AverU, AverV, T, polyOrder,
                             *tauRef, coordinates(0, 0, 0)) 
                    + fneq[xiIdx]/ TotalWeight;
        }
        Rho(0, 0) = AverRho;  // for KerCalcVolumeFraction
        U(0, 0) = AverU;
        V(0, 0) = AverV;
    }   // local
#endif  // OPS_2D
}

void KerSubCollideIBB(ACC<Real>& fStage, const ACC<Real>& f,
                      const ACC<Real>& coordinates, 
                      const ACC<int>& nodeType,
                      const ACC<Real>& Rho, const ACC<Real>& U,
                      const ACC<Real>& V, const Real* tauRef, 
                      const Real* dt,
                      const int* lattIdx, 
                      const ACC<Real>& acceleration) {
#ifdef OPS_2D
    const int vt{nodeType(0, 0)};
    bool CollisionRequired{
        (VertexType_Solid2Fluid == vt || VertexType_NewInterface == vt)};
    if (CollisionRequired) {
        Real rho{Rho(0, 0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        const Real Fx{acceleration(0, 0, 0)};
        const Real Fy{acceleration(1, 0, 0)};
        Real tau = (*tauRef);
        Real Mf[lattIdx[1] + 1]{0};
        Real MfStage[lattIdx[1] + 1]{0};
        //Turbulence model
        Real fstream[lattIdx[1]+1]{0};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            fstream[xiIndex] = f(xiIndex,0,0);
        }
        // moment for f
        for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                Mf[i] += TransformationM[i][xiIndex] * f(xiIndex, 0, 0);
            }
        }
        tau = Smagorinsky_LES(rho, u, v, tau, fstream, lattIdx);
        Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        // relaxation moment
        Real RelaxationRate[lattIdx[1]+1]{1, 
                                        0.3, 
                                        1, 
                                        1, 1.5, 
                                        1, 1.5, 
                                        dtOvertauPlusdt, 
                                        dtOvertauPlusdt};
        Real LinearAbsorb{0};
        const Real iCoordinate{coordinates(0,0,0)};
        if (iCoordinate <= BlockCoor[0]+ DampWidth|| 
            iCoordinate >= BlockCoor[1]-DampWidth) { 

            // Real NonLinearConst{10 * L0};       // unit 1/L
            Real StartPt{0};
            if (iCoordinate <= BlockCoor[0]+ DampWidth) {
                StartPt = BlockCoor[0]+ DampWidth;
            } else if (iCoordinate >= BlockCoor[1]- DampWidth) {
                StartPt = BlockCoor[1]- DampWidth;
            }
            Real DistanceScale = (iCoordinate - StartPt) / DampWidth;
            LinearAbsorb =TimeStep()*
                (- LinearConst * (exp(pow(DistanceScale, ns)) - 1) / (exp(1) - 1));
        }
        
         Real AbsorbWave[lattIdx[1]+1] {1,
                        1,
                        1,
                        1 + LinearAbsorb * 1/RelaxationRate[3],
                        1 + LinearAbsorb * 1/RelaxationRate[4],
                        1 + LinearAbsorb * 1/RelaxationRate[5],
                        1 + LinearAbsorb * 1/RelaxationRate[6],
                        1 + LinearAbsorb * 2/RelaxationRate[7],
                        1 + LinearAbsorb * 2/RelaxationRate[8]};
        // moment space
        
         Real EquilibriumMoment[lattIdx[1]+1]{rho, 
                                        rho * (-2 + (u * u + v * v) * 
                                        (1 + LinearAbsorb * 2/RelaxationRate[1])),
                                        rho * (1 - (u * u + v * v) * 
                                        (1 + LinearAbsorb * 2/RelaxationRate[2])),
                                        rho * u, -rho * u,
                                        rho * v, -rho * v,
                                        rho * (u * u - v * v), 
                                        rho * u * v};
         Real ForceMoment[lattIdx[1]+1]   {0,   
                                    rho * 2 * (Fx * u + Fy * v), 
                                    rho * (-2 * (Fx * u + Fy * v)), 
                                    rho * Fx, -rho * Fx, 
                                    rho * Fy, -rho * Fy, 
                                    rho * 2 * (Fx * u - Fy * v), 
                                    rho * (Fy * u + Fx * v)};
        //  Real EquilibriumMoment[lattIdx[1]+1]{rho, 
        //                                 rho * ((u * u + v * v))* 
        //                                 (1 + LinearAbsorb * 2/RelaxationRate[1]),
        //                                 0,
        //                                 rho * u, 0,
        //                                 rho * v, 0,
        //                                 rho * (u * u - v * v), 
        //                                 rho * u * v};
        //  Real ForceMoment[lattIdx[1]+1]   {0,   
        //                             rho * 2 * (Fx * u + Fy * v), 
        //                             0, 
        //                             rho * Fx, 0, 
        //                             rho * Fy, 0, 
        //                             rho * 2 * (Fx * u - Fy * v), 
        //                             rho * (Fy * u + Fx * v)};
        // MRT collision
        for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
             MfStage[i] =
                    Mf[i] - RelaxationRate[i] * 
                    (Mf[i] - AbsorbWave[i] * EquilibriumMoment[i]) +
                    (1 - RelaxationRate[i] / 2) * (*dt) * ForceMoment[i];
        }
        // inverse to population space
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            Real InversefStage{0};
            for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
                InversefStage += InverseM[xiIndex][i] * MfStage[i];
            }
            fStage(xiIndex, 0, 0) = InversefStage;
#ifdef CPU
            const Real res{fStage(xiIndex, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where rho=%e u=%e v=%e at "
                    "x=%e y=%e\n",
                    res, xiIndex, rho, u, v, coordinates(0, 0, 0),
                    coordinates(1, 0, 0));
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
#endif  // CPU
        }
    }
#endif  // OPS_2D
}

void KerSubStreamIBB(ACC<Real>& f, const ACC<Real>& fStage,
                     const ACC<Real>& coordinates,
                     const ACC<int>& nodeType, 
                     const ACC<Real>& Rho,
                     const ACC<Real>& U, 
                     const ACC<Real>& V,
                     const int* lattIdx,
                     const int* idx, 
                     const ACC<Real>& VolumeFraction,
                     const ACC<int>& SegmentNo,
                     const ACC<int>& Position,
                     const ACC<Real>& NormalDirections) {
#ifdef OPS_2D
    const int vt{nodeType(0, 0)};
    bool StreamRequired{
        (VertexType_Solid2Fluid == vt || VertexType_NewInterface == vt)};
    if ((int)VertexType::Fluid == vt || VertexType_Interface == vt) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor{nodeType(-cx, -cy)};
            if (VertexType_Solid2Fluid == vtneighbor ||
                VertexType_NewInterface == vtneighbor) {
                // StreamRequired = true;
                // break;
                f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
            }
        }
    }
    if (StreamRequired) {
        Real u{0};
        Real v{0};
        Real rho{Rho(0, 0)};
        Real vof{VolumeFraction(0,0)};
        const Real T{1};
        const int polyOrder{2};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor{nodeType(-cx, -cy)};
            // f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
            // // boundary condition for Free surface FSK
            // // and IBB for ImmerseSolid no-slip wall
            if (VertexType_Empty == vtneighbor) {
                // Real IntersectFree[]{-100, -100};
                // // Real NormalDirection[SpaceDim()]{NormalDirections(0, 0, 0),
                // //                                  NormalDirections(1, 0, 0)};
                // bool Subvof = false;
                // // CalcIntersectFree(NormalDirection, VolumeFraction(0, 0), 
                // //                 cx, cy, coordinates(0,0,0), coordinates(1,0,0), 
                // //                   IntersectFree, Subvof);
                // // if(abs(IntersectFree[1]-idx[1])>0.5){
                //     IntersectFree[1] = (idx[1]+(VolumeFraction(0, 0) - 0.5) );
                // // }
                // Real rhoA =
                //     1 + (IntersectFree[1] - (jMax+VOFofInletInterface-0.5) ) 
                //         * CS * TimeStep() * Gravity;
                // if ((int)VertexPosition::In_Chamber == Position(0, 0)) {
                //     rhoA += DeltaRho;
                // }
                // if ((int)VertexType::Fluid == nodeType(cx, cy) ||
                //     (VertexType_Interface == nodeType(cx, cy))) {
                //     Real Distance = vof - 0.5;
                //     u = U(0, 0) * (1 + Distance) - U(cx, cy) * Distance;
                //     v = V(0, 0) * (1 + Distance) - V(cx, cy) * Distance;
                // } else {
                //     u = U(0, 0);
                //     v = V(0, 0);
                // }
                // Real feq = CalcBGKFeq(xiIndex, rhoA, u, v, T, polyOrder);
                // Real feqopp =
                //     CalcBGKFeq(OPP[xiIndex], rhoA, u, v, T, polyOrder);
                // f(xiIndex, 0, 0) = feq + feqopp - fStage(OPP[xiIndex], 0, 0);
                f(xiIndex, 0, 0) = fStage(xiIndex, 0, 0);
            } else if ((int)VertexType::ImmersedSolid == vtneighbor ||
                        VertexType_Solid2Interface == vtneighbor) {
                // f(xiIndex, 0, 0) = fStage(xiIndex, 0, 0);
                Real Intersect[]{-100, -100};
                Real Distance{-100};
                const Real x{coordinates(0, 0, 0)};
                const Real y{coordinates(1, 0, 0)};
                const Real CentralX{RotationCenter[0]};
                const Real CentralY{RotationCenter[1]};
#ifdef CIRCLE
                CalcIntersectIBB(-cx, -cy, x, y, 
                              RotationCenter, Radii, Intersect, Distance);
#endif //CIRCLE
#ifdef SQUARE
                const Real P[]{x, y};
                CalcIntersectIBB(SegmentNo(0, 0), -CS * TimeStep() * cx,
                                 -CS * TimeStep() * cy, P, SolidVertex,
                                 SolidPlaneNormal, Polygon, NumofPol, Intersect,
                                 Distance);
#endif  // SQUARE
                Real uSolid = VelocityofRoCenter[0] -
                              VelocityofRoCenter[2] * (Intersect[1] - CentralY);
                Real vSolid = VelocityofRoCenter[1] +
                              VelocityofRoCenter[2] * (Intersect[0] - CentralX);
                if ((int)VertexType::Fluid != nodeType(cx, cy) &&
                    VertexType_Interface != nodeType(cx, cy)) {
                    f(xiIndex, 0, 0) = 
                            (Distance * fStage(xiIndex, 0, 0) +
                             Distance * fStage(OPP[xiIndex], 0, 0) +
                             (1 - Distance) * fStage(OPP[xiIndex], 0, 0) +
                             2 * WEIGHTS[xiIndex] * rho *
                                 (cx * uSolid + cy * vSolid) * CS) /
                            (Distance + 1);
                } else {
                    // if ((int)VertexType::Fluid != nodeType(2 * cx, 2 * cy) &&
                    //     VertexType_Interface != nodeType(2 * cx, 2 * cy)) {
                        f(xiIndex, 0, 0) =
                            (Distance * fStage(xiIndex, 0, 0) +
                             Distance * fStage(OPP[xiIndex], 0, 0) +
                             (1 - Distance) * fStage(OPP[xiIndex], cx, cy) +
                             2 * WEIGHTS[xiIndex] * rho *
                                 (cx * uSolid + cy * vSolid) * CS) /
                            (Distance + 1);
                    // } else {
                    //     Real f_xw = Distance * (1 + Distance) / 2 *
                    //                     fStage(OPP[xiIndex], 0, 0) +
                    //                 (1 - Distance) * (1 + Distance) *
                    //                     fStage(OPP[xiIndex], cx, cy) -
                    //                 Distance * (1 - Distance) / 2 *
                    //                     fStage(OPP[xiIndex], 2 * cx, 2 * cy) +
                    //                 2 * WEIGHTS[xiIndex] * rho *
                    //                     (cx * uSolid + cy * vSolid) * CS;
                    //     f(xiIndex, 0, 0) =
                    //         2 / (1 + Distance) / (2 + Distance) * f_xw +
                    //         2 * Distance / (1 + Distance) *
                    //             fStage(xiIndex, 0, 0) -
                    //         Distance / (2 + Distance) *
                    //             fStage(xiIndex, cx, cy);
                    // }
                }
            }else if ((int)VertexType::Wall == vtneighbor||
                        (int)VertexType::OutLet == vtneighbor){
                f(xiIndex, 0, 0) = fStage(xiIndex, 0, 0);
            } else{
                 f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
            }
            const Real res{f(xiIndex, 0, 0)};
            if (isnan(res) || res <= 0 || isinf(res)) {
                ops_printf(
                    "Error! Distribution function = %e becomes invalid at  "
                    "the lattice %i where rho=%e u=%e v=%e at "
                    "x=%e y=%e vtneighbor=%d vt=%d\n",
                    res, xiIndex, rho, u, v, coordinates(0, 0, 0),
                    coordinates(1, 0, 0),vtneighbor,vt);
                assert(!(isnan(res) || res <= 0 || isinf(res)));
            }
        }
    }
#endif  // OPS_2D
}

void KerFinalDFsIBB(ACC<Real>& f, const ACC<Real>& fStage,
                     const ACC<int>& nodeType, 
                     const int* lattIdx) {
#ifdef OPS_2D
    const int vt{nodeType(0, 0)};
    bool StreamRequired{
        (VertexType_Solid2Fluid == vt || VertexType_NewInterface == vt)};
    if (StreamRequired) {
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            f(xiIdx, 0, 0) = fStage(xiIdx, 0, 0);
        }
    }
#endif  // OPS_2D
}

void KerSubCalcMacroVarsIBB(ACC<Real>& Rho,
                        ACC<Real>& U, 
                        ACC<Real>& V,
                        ACC<Real>& Mass,
                        ACC<Real>& VolumeFraction,  
                        const ACC<int>& nodeType,
                        const ACC<Real>& Porosity,
                        const ACC<Real>& f, 
                        const int* lattIdx,
                        const ACC<Real>& acceleration,
                        ACC<Real>& ResidualError,
                        const bool* VOFCalcuRequired) {
#ifdef OPS_2D
    const int vt{nodeType(0, 0)};
    ResidualError(0, 0) =0;
    Real vof{0};
    bool StreamRequired{
        (VertexType_Solid2Fluid == vt || VertexType_NewInterface == vt)};
    if ((int)VertexType::Fluid == vt || VertexType_Interface == vt) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor{nodeType(cx, cy)};
            if (VertexType_Solid2Fluid == vtneighbor ||
                VertexType_NewInterface == vtneighbor) {
                StreamRequired = true;
                break;
            }
        }
    }
    if (StreamRequired) {
        Real u{0};
        Real v{0};
        Real rho{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            Real cx = CS * XI[xiIdx * LATTDIM];
            Real cy = CS * XI[xiIdx * LATTDIM + 1];
            rho += f(xiIdx, 0, 0);
            u += cx * f(xiIdx, 0, 0);
            v += cy * f(xiIdx, 0, 0);
        }
        if (VertexType_Solid2Fluid == vt || 
            VertexType_NewInterface == vt) {
            ResidualError(0, 0) = 
                        sqrt((rho - Rho(0, 0)) * (rho - Rho(0, 0))) /
                          sqrt((rho + Rho(0, 0)) * (rho + Rho(0, 0)));
            // ops_printf("Error= %e\n", ResidualError(0, 0));
        }
        Rho(0, 0) = rho;
        if (*VOFCalcuRequired){
            U(0, 0) = u / rho + (TimeStep() * acceleration(0, 0, 0) / 2);
            V(0, 0) = v / rho + (TimeStep() * acceleration(1, 0, 0) / 2);
        }
        
    }
    if (VertexType_NewInterface == vt&&(*VOFCalcuRequired)) {
        int NumofNeighLiquid{0};
        if (VertexType_Interface == nodeType(-1, 0) ||
            (int)VertexType::Fluid == nodeType(-1, 0)) {
            vof += VolumeFraction(-1, 0);
            NumofNeighLiquid++;
        }
        if (VertexType_Interface == nodeType(1, 0) ||
            (int)VertexType::Fluid == nodeType(1, 0)) {
            vof += VolumeFraction(1, 0);
            NumofNeighLiquid++;
        }
        if (NumofNeighLiquid != 0) {
            vof /= NumofNeighLiquid;
            VolumeFraction(0, 0) = vof;
        } else {
            VolumeFraction(0, 0) = 0;
        }
        Mass(0, 0) = VolumeFraction(0, 0) * 
                        Porosity(0, 0) * Rho(0, 0);
    }
#endif  // OPS_2D
}

void KerTotalDensityError(const ACC<int>& nodeType, 
                        const ACC<Real>& ResidualError,
                        double* Error) {
    const int vt{nodeType(0, 0)};
    if (VertexType_Solid2Fluid == vt || 
        VertexType_NewInterface == vt) {
        *Error += ResidualError(0, 0);
    }
}

void KerSetDFsofNewInterface(ACC<Real>& f, const ACC<Real>& fStage, 
                      const ACC<int>& nodeType,
                      const int* lattIdx) {
#ifdef OPS_2D
    const int vt{nodeType(0, 0)};
    if (VertexType_NewInterface == vt) {
         for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            f(xiIdx,0,0) = fStage(xiIdx,0,0);
        }
    }
#endif  // OPS_2D
}

void KerNonEQExtrapolation(ACC<Real>& f,
                           const ACC<int>& nodeType,
                           const ACC<int>& geometryProperty,
                           const Real *givenMacroVars,
                           const ACC<Real>& acceleration,
                           const int* lattIdx,
                           const ACC<Real>& coordinates,
                           const int* iter) {
#ifdef OPS_2D
    if ((int)VertexType::Inlet == nodeType(0, 0) ||
        (int)VertexType::OutLet == nodeType(0, 0)) {
        VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
        const Real T{1.0};
        const int polyOrder{2};
        int outgoingX{0};
        int outgoingY{0};
        switch (vg) {
            case VG_IP: {
                outgoingX = 1;
                outgoingY = 0;
            } break;
            case VG_IM: {
                outgoingX = -1;
                outgoingY = 0;
            } break;
            case VG_JP: {
                outgoingX = 0;
                outgoingY = 1;
            } break;
            case VG_JM: {
                outgoingX = 0;
                outgoingY = -1;
            } break;
            case VG_IPJP_I: {
                outgoingX = 1;
                outgoingY = 1;
            } break;
            case VG_IPJM_I: {
                outgoingX = 1;
                outgoingY = -1;
            } break;
            case VG_IMJM_I: {
                outgoingX = -1;
                outgoingY = -1;
            } break;
            case VG_IMJP_I: {
                outgoingX = -1;
                outgoingY = 1;
            } break;
            default:
                break;
        }
        Real rhoOutgoing{0};
        Real uOutgoing{0};
        Real vOutgoing{0};
        if (VertexType_Interface == nodeType(outgoingX, outgoingY) ||
            (int)VertexType::Fluid == nodeType(outgoingX, outgoingY)) {
            for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
                Real cx = CS * XI[xiIdx * LATTDIM];
                Real cy = CS * XI[xiIdx * LATTDIM + 1];
                rhoOutgoing += f(xiIdx, outgoingX, outgoingY);
                uOutgoing += cx * f(xiIdx, outgoingX, outgoingY);
                vOutgoing += cy * f(xiIdx, outgoingX, outgoingY);
                const Real res{f(xiIdx, outgoingX, outgoingY)};
                if (isnan(res) || res <= 0 || isinf(res)) {
                        ops_printf(
                            "Error! Distribution function = %e becomes invalid "
                            "at  "
                            "the lattice %i, vt=%d \n",
                            res, xiIdx, nodeType(outgoingX, outgoingY));
                        assert(!(isnan(res) || res <= 0 || isinf(res)));
                }
            }

            uOutgoing =
                uOutgoing / rhoOutgoing +
                (TimeStep() * acceleration(0, outgoingX, outgoingY) / 2);
            vOutgoing =
                vOutgoing / rhoOutgoing +
                (TimeStep() * acceleration(1, outgoingX, outgoingY) / 2);
            Real u{uOutgoing};
            Real v{vOutgoing};
            Real rho{rhoOutgoing};

            if ((int)VertexType::Inlet == nodeType(0, 0)) {
                // Real y{coordinates(1,0,0)};
                u = givenMacroVars[0];
                v = givenMacroVars[1];
            }
            if ((int)VertexType::OutLet == nodeType(0, 0)) {
                // v = 0;
                // rho = 1;
                // if (VG_JP==vg||VG_JM==vg){
                //     u = 0;
                //     v = 0;
                // }
            }
            for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
                Real feq = CalcBGKFeq(xiIdx, rho, u, v, T, polyOrder);
                Real feqOutgoing = CalcBGKFeq(xiIdx, rhoOutgoing, uOutgoing,
                                              vOutgoing, T, polyOrder);
                f(xiIdx, 0, 0) =
                    feq - feqOutgoing + f(xiIdx, outgoingX, outgoingY);
            }
        }
    }
#endif
}

void KerExtrapolation2Order(ACC<Real>& f,
                           const ACC<int>& nodeType,
                           const ACC<int>& geometryProperty,
                           const int* lattIdx,
                           const ACC<Real>& U, 
                           const ACC<Real>& V,
                           const ACC<Real>& fStage) {
#ifdef OPS_2D
    if (
        (int)VertexType::OutLet == nodeType(0, 0)) {
        VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
        int outgoingX{0};
        int outgoingY{0};
        switch (vg) {
            case VG_IP: {
                outgoingX = 1;
                outgoingY = 0;
            } break;
            case VG_IM: {
                outgoingX = -1;
                outgoingY = 0;
            } break;
            case VG_JP: {
                outgoingX = 0;
                outgoingY = 1;
            } break;
            case VG_JM: {
                outgoingX = 0;
                outgoingY = -1;
            } break;
            default:
                break;
        }
        if (VertexType_Interface == nodeType(outgoingX, outgoingY) ||
            (int)VertexType::Fluid == nodeType(outgoingX, outgoingY)) {
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            int cx = (int)XI[xiIdx * LATTDIM];
            int cy = (int)XI[xiIdx * LATTDIM + 1];
            // if (cx == outgoingX){//&&cy ==0
            // //     continue;
            // // }else if ((cx == outgoingX&&cy !=0)){
            //     f(xiIdx, 0, 0) =f(xiIdx, outgoingX, outgoingY);
            // // //     //2*f(xiIdx, outgoingX, outgoingY)-f(xiIdx, 2*outgoingX, 2*outgoingY);
            // }else{
                f(xiIdx, 0, 0) =f(xiIdx, outgoingX, outgoingY);
                // f(xiIdx, 0, 0) =fStage(xiIdx,-cx,-cy);
                // f(xiIdx, 0, 0) =2*f(xiIdx, outgoingX, outgoingY)-f(xiIdx, 2*outgoingX, 2*outgoingY);
            // }   
        }
        }
    }
#endif
}

void KerConvectiveBoundary(ACC<Real>& f,
                           const ACC<int>& nodeType,
                           const ACC<int>& geometryProperty,
                           const Real *givenMacroVars,
                           const int* lattIdx,
                           const ACC<Real>& coordinates) {
#ifdef OPS_2D
    if ((int)VertexType::OutLet == nodeType(0, 0)) {
        VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);

        int outgoingX{0};
        int outgoingY{0};
        switch (vg) {
            case VG_IP: {
                outgoingX = 1;
                outgoingY = 0;
            } break;
            case VG_IM: {
                outgoingX = -1;
                outgoingY = 0;
            } break;
            case VG_JP: {
                outgoingX = 0;
                outgoingY = 1;
            } break;
            case VG_JM: {
                outgoingX = 0;
                outgoingY = -1;
            } break;
            case VG_IPJP_I: {
                outgoingX = 1;
                outgoingY = 1;
            } break;
            case VG_IPJM_I: {
                outgoingX = 1;
                outgoingY = -1;
            } break;
            case VG_IMJM_I: {
                outgoingX = -1;
                outgoingY = -1;
            } break;
            case VG_IMJP_I: {
                outgoingX = -1;
                outgoingY = 1;
            } break;
            default:
                break;
        }
        if (VertexType_Interface == nodeType(outgoingX, outgoingY) ||
            (int)VertexType::Fluid == nodeType(outgoingX, outgoingY)) {
            for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
                // Real gama {givenMacroVars[0]/CS};
                //  in wave case
                Real uOutgoing{0};
                for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++){
                    Real cx = CS * XI[xiIdx * LATTDIM];
                    uOutgoing += cx * f(xiIdx, outgoingX, outgoingY);
                }
                Real gama {-outgoingX*uOutgoing};
                // f(xiIdx, 0, 0) =
                //     (f(xiIdx, 0, 0) + gama * f(xiIdx, outgoingX, outgoingY))/(1+gama);
                    // 2 order scheme
                f(xiIdx, 0, 0) =
                    (f(xiIdx, 0, 0) + 2*gama * f(xiIdx, outgoingX, outgoingY)-
                        0.5 * gama*f(xiIdx, 2*outgoingX, 2*outgoingY))/(1+1.5*gama);    
            }
        }
    }
#endif
}



#endif
