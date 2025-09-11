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

// The code is based on the algorithm proposed by N. Thürey, U. Rüde, C. Körner
// reference: Interactive Free Surface Fluids with the Lattice Boltzmann Method

#ifndef FREESURFACE2D_KERNEL_INC
#define FREESURFACE2D_KERNEL_INC
#include "ops_lib_core.h"

#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif
#include "model.h"
#include "type.h"
#include "flowfield_host_device.h"
#include "dam2d_host_device.h"
#define VertexType_Interface 2000
#define VertexType_Empty 2001
#define VertexType_NewInterface 2003
#define VertexType_Interface2Empty 2004
#define VertexType_Interface2Fluid 2005
#define VertexType_SlipWall 2006
const int BoundaryScheme_HalfwayBounceBack{100};
const int BoundaryScheme_HalfwayBounceForward{101};
#define FreeMpi_2D // an identification for parellel computing

// the dam size
const int iMax {102};
const int jMax {204};

enum VertexQuality { No_Fluid = 1, No_Empty = 2, Standard = 3 };

void KerSetInitialMacroVars(ACC<Real>& rho, ACC<Real>& u, ACC<Real>& v,
                            ACC<Real>& mass, ACC<Real>& vof, 
                            ACC<int>& Quality,
                            const ACC<int>& nodeType,
                            const ACC<Real>& coordinates, const int* idx,
                            const ACC<Real>& acceleration) {
    const int vt {nodeType(0, 0)};
    const int i{idx[0]};
    const int j{idx[1]};
    rho(0, 0) =
            1 + CS * TimeStep() * (idx[1] - jMax ) * acceleration(1, 0, 0);
    u(0, 0) = 0;
    v(0, 0) = 0;
    //  Initialise the mass and vof here according to x,y
    if ((int)VertexType::Fluid == vt) {
        vof(0, 0) = 1.0;
    } else if (VertexType_Interface == vt) {
        if (iMax == i && (j < jMax && j > 0)) {
            vof(0, 0) = 0.4;
        } else if (jMax == j && (i < iMax && i > 0)) {
            vof(0, 0) = 0.8;
        } else if (iMax == i && jMax == j) {
            vof(0, 0) = 0.16;
        }
    } else {
        vof(0, 0) = 0.0;
    }
    // dam breaking from right boundary
    // const int iMin = 513 - iMax;
    // if ((int)VertexType::Fluid == vt) {
    //     vof(0, 0) = 1.0;
    // } else if (VertexType_Interface == vt) {
    //     if (iMin == i && (j < jMax && j > 0)) {
    //         vof(0, 0) = 0.4;
    //     } else if (jMax == j && (i > iMin)) {
    //         vof(0, 0) = 0.8;
    //     } else if (iMin == i && jMax == j) {
    //         vof(0, 0) = 0.16;
    //     }
    // } else {
    //     vof(0, 0) = 0.0;
    // }
    mass(0, 0) = vof(0, 0) * rho(0, 0);
    Quality(0,0) = (int)Standard;
    //define the cell quality
}

// determine the initial interface location
void KerSetInitialNodeType(ACC<int>& nodeType, 
                           const ACC<Real>& coordinates,const int* idx) {
    // const Real x{coordinates(0, 0, 0)};
    // const Real y{coordinates(1, 0, 0)};
    const int i{idx[0]};
    const int j{idx[1]};
    nodeType(0, 0) = VertexType_Empty;
    // dam breaking from left boundary
    if ((i < iMax && i > 0) && (j < jMax && j > 0)) {
        nodeType(0, 0) = (int)VertexType::Fluid;
    } else if ((iMax == i && (j <= jMax && j > 0)) ||
               (jMax == j && (i <= iMax && i > 0))) {
        nodeType(0, 0) = VertexType_Interface;
    } else {
        nodeType(0, 0) = VertexType_Empty;
    }
    // dam breaking from right boundary
    // const int iMin = 513-iMax;
    // if ((i > iMin ) && (j < jMax && j > 0)) {
    //     nodeType(0, 0) = (int)VertexType::Fluid;
    // } else if ((iMin == i && (j <= jMax && j > 0)) ||
    //            (jMax == j && (i >= iMin))) {
    //     nodeType(0, 0) = VertexType_Interface;
    // } else {
    //     nodeType(0, 0) = VertexType_Empty;
    // }
}

// the kernel is just to specified which boundary is the slip boundary.
// if a new slip boundaryscheme is defined,the kernel can be omited.
void KerSetSlipBoundaryNodeType(ACC<int>& nodeType,
                                const ACC<Real>& coordinates,
                                const ACC<int>& geometry) {
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0);
    // Left wall
    if (VG_IP == vg) {
        nodeType(0, 0) = VertexType_SlipWall;
    }
    // Right wall
    if (VG_IM == vg) {
        nodeType(0, 0) = VertexType_SlipWall;
    }
    // Bottom wall
    if (VG_JP == vg) {
        nodeType(0, 0) = VertexType_SlipWall;
    }
    // Top wall
    if (VG_JM == vg) {
        nodeType(0, 0) = VertexType_SlipWall;
    }
    //the corner is still no slip wall.
}

void KerUpdateMacroBodyForce(ACC<Real>& acceleration,
                             const ACC<Real>& coordinates, const int* idx) {
    // here we put some constant for gravity to the following two lines
    acceleration(0, 0, 0) = 0;
    acceleration(1, 0, 0) = -4.2747e-4;
}

// no difference between original 'KerCalcDensity' and current 'KerCalcDensityfree'
void KerCalcDensityfree(ACC<Real>& Rho, const ACC<Real>& f,
                    const ACC<int>& nodeType, const int* lattIdx) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    if (vt != (int)VertexType::ImmersedSolid) {
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

// VertexType_Interface == vt is added compared with original KerCalcUForce
void KerCalcUForcefree(ACC<Real>& U, const ACC<Real>& f, const ACC<int>& nodeType,
                   const ACC<Real>& coordinates, const ACC<Real>& acceleration,
                   const ACC<Real>& Rho, const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    const Real x{coordinates(0, 0, 0)};
    const Real y{coordinates(1, 0, 0)};
    const int vt {nodeType(0, 0)};
    if (vt != (int) VertexType::ImmersedSolid) {
        Real u{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            u += CS * XI[xiIdx * LATTDIM] * f(xiIdx, 0, 0);
        }
        u /= Rho(0, 0);
        if ((int)VertexType::Fluid == vt || (int) VertexType::MDPeriodic == vt ||
            VertexType_Interface == vt) {
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

// VertexType_Interface == vt is added compared with original KerCalcVForce
void KerCalcVForcefree(ACC<Real>& V, const ACC<Real>& f, const ACC<int>& nodeType,
                   const ACC<Real>& coordinates, const ACC<Real>& acceleration,
                   const ACC<Real>& Rho, const Real* dt, const int* lattIdx) {
#ifdef OPS_2D
    const Real x{coordinates(0, 0, 0)};
    const Real y{coordinates(1, 0, 0)};
    const int vt {nodeType(0, 0)};
    if (vt != (int) VertexType::ImmersedSolid) {
        Real v{0};
        for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
            v += CS * XI[xiIdx * LATTDIM + 1] * f(xiIdx, 0, 0);
        }
        v /= Rho(0, 0);
        if ((int)VertexType::Fluid == vt || (int) VertexType::MDPeriodic == vt||
            VertexType_Interface == vt) {
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

// VertexType_Interface == vt
void KerCalcBodyForce1STfree(ACC<Real>& fStage, const ACC<Real>& acceleration,
                         const ACC<Real>& Rho, const ACC<int>& nodeType,
                         const int* lattIdx, const ACC<Real>& U,const ACC<Real>& V) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    if (vt == (int)VertexType::Fluid || vt == (int)VertexType::MDPeriodic||
        VertexType_Interface == vt) {
        Real rho{Rho(0, 0)};
        Real u{U(0, 0)};
        Real v{V(0, 0)};
        Real g[]{acceleration(0, 0, 0), acceleration(1, 0, 0)};
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

// VertexType_Interface == vt is added
void KerStreamfree(ACC<Real>& f, const ACC<Real>& fStage, const ACC<int>& nodeType,
               const ACC<int>& geometry, const int* lattIdx) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    VertexGeometryType vg = (VertexGeometryType)geometry(0, 0);
    for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
        int cx = (int)XI[xiIndex * LATTDIM];
        int cy = (int)XI[xiIndex * LATTDIM + 1];

        if (vt == (int)VertexType::Fluid || vt == (int)VertexType::MDPeriodic ||
            vt == (int)VertexType::VirtualBoundary||
            VertexType_Interface == vt) {
            f(xiIndex, 0, 0) = fStage(xiIndex, -cx, -cy);
            continue;
        }
        // Block boundary and immersed boundary
        if (vt != (int)VertexType::ImmersedSolid && vt != (int)VertexType::Fluid &&
            vt != (int)VertexType::VirtualBoundary && vt != (int)VertexType::MDPeriodic && 
            VertexType_Interface != vt) {
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

// turbulent model is used and VertexType_Interface== vt is added
void KerCollideBGKFreeSurface(ACC<Real>& fStage, const ACC<Real>& f,
                             const ACC<Real>& coordinates,
                             const ACC<int>& nodeType, const ACC<Real>& Rho,
                             const ACC<Real>& U, const ACC<Real>& V,
                             const Real* tauRef, const Real* dt,
                             const int* lattIdx,ACC<Real>&Taucopy) {
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
        const Real MeshSize = CS*(*dt);
        const Real C{0.18*0.18};
        Real NonEquilibriumMoment[4]{0,0,0,0};
        Real NonEquilibriumMoment2{0.0};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            Real cx = CS * XI[xiIndex * LATTDIM];
            Real cy = CS * XI[xiIndex * LATTDIM + 1];
            const Real feq{CalcBGKFeq(xiIndex, rho, u, v, T, polyOrder)};
            NonEquilibriumMoment[0] += (cx * cx * (f(xiIndex, 0, 0) - feq));
            NonEquilibriumMoment[1] += (cx * cy * (f(xiIndex, 0, 0) - feq));
            NonEquilibriumMoment[2] += (cy * cx * (f(xiIndex, 0, 0) - feq));
            NonEquilibriumMoment[3] += (cy * cy * (f(xiIndex, 0, 0) - feq));
        }
        for (int i=0;i<4;i++){
            NonEquilibriumMoment2 +=
                2 * NonEquilibriumMoment[i] * NonEquilibriumMoment[i];
        }
        tau = 0.5 *(tau - 0.5 * (*dt) +
               sqrt((tau + 0.5 * (*dt)) * (tau + 0.5 * (*dt)) +
                    2 * C * MeshSize * MeshSize * sqrt(NonEquilibriumMoment2)/rho));
        //Turbulence model
        
        Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            const Real feq{CalcBGKFeq(xiIndex, rho, u, v, T, polyOrder)};
            if (vt == (int)VertexType::Fluid || vt == (int)VertexType::MDPeriodic||
                VertexType_Interface== vt) {
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

// the current interface node has no fluid neighbor or no gas neighbor.
void UpdatNodeQuality(ACC<int>& Quality,const ACC<int>& nodeType,  
                      const int* lattIdx) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    Quality(0,0) = (int)Standard;
    if (VertexType_Interface == vt) {
        int NumofFluid{0};
        int NumofEmpty{0};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor=nodeType(cx, cy);
            if ((int)VertexType::Fluid == vtneighbor) {
                NumofFluid++;
            } else if (VertexType_Empty == vtneighbor) {
                NumofEmpty++;
            }
        }
        //determine the quality of the interface cells:
        if ((0==NumofFluid)&&NumofEmpty>0) {
            Quality(0,0) = (int)No_Fluid ;
        } else if ((0==NumofEmpty)&&NumofFluid>0) {
            Quality(0,0)=(int)No_Empty;
        } else {Quality(0,0)=(int)Standard;}
    }
#endif  // OPS_2D    
}

void KerMassExchange(ACC<Real>& Mass, const ACC<Real>& fStage,
                     const ACC<Real>& VolumeFraction, const ACC<int>& nodeType,
                     const ACC<int>& Quality, const ACC<int>& geometryProperty,
                     const ACC<Real>& coordinates, const int* lattIdx) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    VertexQuality vq = (VertexQuality)Quality(0, 0);
    if (VertexType_Interface == vt) {
        const Real vof {VolumeFraction(0, 0)} ;
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            Real deltam{0.0};
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            if (0 == cx && 0 == cy) {
                continue;
            }
            const int vtneighbor {nodeType(cx, cy)};
            switch (vtneighbor) {
                case VertexType_Interface: {
                    VertexQuality vqNeighbor = (VertexQuality)Quality(cx, cy);
                    //Force non-physical (or undesired) spatter and bubble
                    //to fill or empty: first step
                    if ((Standard == vq && Standard == vqNeighbor) ||
                        (No_Fluid == vq && No_Fluid == vqNeighbor) ||
                        (No_Empty == vq && No_Empty == vqNeighbor)) {
                        deltam = 0.5 * (vof + VolumeFraction(cx, cy)) *
                                 (fStage(OPP[xiIndex], cx, cy) -fStage(xiIndex, 0, 0));
                    } else if ((Standard == vq && No_Empty == vqNeighbor) ||
                               (No_Fluid == vq && Standard == vqNeighbor) ||
                               (No_Fluid == vq && No_Empty == vqNeighbor)) {
                        deltam = 0.5 * (vof + VolumeFraction(cx, cy)) *
                                 (-fStage(xiIndex, 0, 0));
                    } else if ((Standard == vq && No_Fluid == vqNeighbor) ||
                               (No_Empty == vq && Standard == vqNeighbor) ||
                               (No_Empty == vq && No_Fluid == vqNeighbor)) {
                        deltam = 0.5 * (vof + VolumeFraction(cx, cy)) *
                                 (fStage(OPP[xiIndex], cx, cy));
                    }
                    break;
                }
                case (int)VertexType::Fluid: {
                    deltam =
                        (fStage(OPP[xiIndex], cx, cy) - fStage(xiIndex, 0, 0));
                    break;
                }
                case VertexType_Empty: {
                    deltam = 0.0;
                    break;
                }
                case VertexType_SlipWall:{
                    const int vgneighbor = geometryProperty(cx, cy);
                    if (VG_IP == vgneighbor || VG_IM == vgneighbor) {
                        if (VertexType_Empty == nodeType(0, cy)) {
                            deltam = 0;
                        } else if (VertexType_Interface == nodeType(0, cy)) {
                            deltam = 0.5 * (vof + VolumeFraction(0, cy)) *
                                     (fStage(OPP[xiIndex], cx, cy) -
                                      fStage(xiIndex, 0, 0));
                        } else if ((int)VertexType::Fluid == nodeType(0, cy)) {
                            deltam = (fStage(OPP[xiIndex], cx, cy) -
                                      fStage(xiIndex, 0, 0));
                        }
                    } else if (VG_JP == vgneighbor || VG_JM == vgneighbor) {
                        if (VertexType_Empty == nodeType(cx, 0)) {
                            deltam = 0;
                        } else if (VertexType_Interface == nodeType(cx, 0)) {
                            deltam = 0.5 * (vof + VolumeFraction(cx, 0)) *
                                     (fStage(OPP[xiIndex], cx, cy) -
                                      fStage(xiIndex, 0, 0));
                        } else if ((int)VertexType::Fluid == nodeType(cx, 0)) {
                            deltam = (fStage(OPP[xiIndex], cx, cy) -
                                      fStage(xiIndex, 0, 0));
                        }
                    } else {
                        deltam = 0.0;
                    }
                    break;
                }
                default: {
                    deltam = 0.0;
                    // exchange with other boundary cells
                    break;
                }
            }
            Mass(0, 0) += deltam;
        }
    }
#endif  // OPS_2D
}

// FSK boundary condition implementation 
void KerReconstructDF(ACC<Real>& f, const ACC<Real>& fStage,
                     const ACC<Real>& VolumeFraction, const ACC<int>& nodeType,
                     const ACC<int>& geometryProperty,
                     const ACC<Real>& U, const ACC<Real>& V, 
                     const ACC<Real>& coordinates, const int* lattIdx,
                     const ACC<Real>& acceleration) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    if (VertexType_Interface == vt) {
        const Real rhoA {1.0} ;//the air pressure in free surface = 0
        const Real T{1.0};
        const int polyOrder{2};
        Real NormalDirection[2]{};
        int vgneighbors[3][3]{};
        Real vofneighbors[3][3]{};
        // bool NexttoBoundary{false};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            vgneighbors [cx + 1][cy + 1] = geometryProperty(cx, cy);
            vofneighbors[cx + 1][cy + 1] = VolumeFraction(cx, cy);
            // const int vtneighbor {nodeType(cx, cy)};
            // To determine if the interface cell is next to boundary wall:
            // if ((int)VertexType::Fluid != vtneighbor &&
            //     VertexType_Empty != vtneighbor &&
            //     VertexType_Interface != vtneighbor) {
            //     NexttoBoundary = true;
            // }
        }
        CalcNormalDirection(NormalDirection, vgneighbors, vofneighbors,
                            lattIdx);
        // reconstruct the missing DFs from empty cells and n*c<0
        // &&((int)VertexType::Fluid == vtneighbor ||
        // VertexType_Interface == vtneighbor)--> do not revise DF from wall
        // || NexttoBoundary-->reconsturct all DF interface next to wall
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor{nodeType(cx, cy)};
            if (cx == 0 && cy == 0) {
                continue;
            }
            // bool neighborisSlipboundary{false};
            // if (VertexType_SlipWall == vtneighbor) {
            //     VertexGeometryType vgneighbor =
            //         (VertexGeometryType)geometryProperty(cx, cy);
            //     if (VG_IP == vgneighbor || VG_IM == vgneighbor) {
            //         neighborisSlipboundary =
            //             (VertexType_Empty == nodeType(0, cy));
            //     } else if (VG_JP == vgneighbor || VG_JM == vgneighbor) {
            //         neighborisSlipboundary =
            //             (VertexType_Empty == nodeType(cx, 0));
            //     }
            // }
            Real NormalDirDotLatticeVel =
                NormalDirection[0] * cx + NormalDirection[1] * cy;
            if (//NormalDirDotLatticeVel > 1e-10 ||
                VertexType_Empty == vtneighbor) {
                Real u = U(0, 0);
                Real v = V(0, 0);  //- TimeStep() * acceleration(1, 0, 0) / 2;
                // Real Distance = VolumeFraction(0, 0) - 0.5;
                // if (VolumeFraction(0, 0) >= 0.5) {
                //     if ((int)VertexType::Fluid == nodeType(-cx, -cy) ||
                //         (VertexType_Interface == nodeType(-cx, -cy) &&
                //          VolumeFraction(-cx, -cy) >= 0.5)) {
                //         u = U(0, 0) * (1 + Distance) - U(-cx, -cy) * Distance;
                //         v = V(0, 0) * (1 + Distance) - V(-cx, -cy) * Distance;
                //     }
                // } else {
                //     if ((int)VertexType::Fluid == nodeType(-cx, -cy) ||
                //         (VertexType_Interface == nodeType(-cx, -cy) &&
                //          VolumeFraction(-cx, -cy) >= 0.5)) {
                //         // u = U(-cx, -cy);
                //         // v = V(-cx, -cy);
                //         if ((int)VertexType::Fluid ==
                //                 nodeType(-cx * 2, -cy * 2) ||
                //             (VertexType_Interface ==
                //                  nodeType(-cx * 2, -cy * 2) &&
                //              VolumeFraction(-cx * 2, -cy * 2) >= 0.5)) {
                //             u = U(-cx, -cy) * (2 + Distance) -
                //                 U(-cx * 2, -cy * 2) * (1 + Distance);
                //             v = V(-cx, -cy) * (2 + Distance) -
                //                 V(-cx * 2, -cy * 2) * (1 + Distance);
                //         }
                //     }
                // }
                if ((int)VertexType::Fluid == nodeType(-cx, -cy) ||
                    (VertexType_Interface == nodeType(-cx, -cy) )) {
                        // &&
                        // VolumeFraction(-cx, -cy) >= 0.5)
                        Real Distance = VolumeFraction(0, 0) - 0.5;
                        u = U(0, 0) * (1 + Distance) - U(-cx, -cy) * Distance;
                        v = V(0, 0) * (1 + Distance) - V(-cx, -cy) * Distance;
                }
                Real feq = CalcBGKFeq(xiIndex, rhoA, u, v, T, polyOrder);
                Real feqopp =
                    CalcBGKFeq(OPP[xiIndex], rhoA, u, v, T, polyOrder);
                f(OPP[xiIndex], 0, 0) = feq + feqopp - fStage(xiIndex, 0, 0);
            }
        }
    }
#endif  // OPS_2D
}

// change the flag of those interface cells that are filled or emptied
void KerPreUpdateFlag(ACC<int>& nodeType, const ACC<Real>& VolumeFraction,
                      const ACC<int>& Quality,const int* lattIdx) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    const int vq {Quality(0, 0)};
    if (VertexType_Interface == vt) {
        const Real vof {VolumeFraction(0, 0)};
        const Real OFFSET {0.001} ;  // should be a global variable
        // Force non-physical (or undesired) spatter and bubble
        // to fill or empty: second step
        if ((vof < 0.0 - OFFSET) || ((int)No_Fluid == vq && vof < 0.1)) {
            bool neighborisI2F{false};
            // for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            //     int cx = (int)XI[xiIndex * LATTDIM];
            //     int cy = (int)XI[xiIndex * LATTDIM + 1];
            //     const Real vofneighbor {VolumeFraction(cx, cy)};
            //     if (vofneighbor > 1.0 + OFFSET ||
            //         ((int)No_Empty == Quality(cx, cy) && vofneighbor > 0.9)) {
            //         neighborisI2F = true;
            //         break;
            //     }
            // }
            // do not change the flag of (interface next to the interface2fluid)
            // to interface2empty
            if (!neighborisI2F) {
                nodeType(0, 0) = VertexType_Interface2Empty;
            }
        } else if ((vof > 1.0 + OFFSET) || ((int)No_Empty == vq && vof > 0.9)) {
            bool neighborisI2E{false};
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                int cx = (int)XI[xiIndex * LATTDIM];
                int cy = (int)XI[xiIndex * LATTDIM + 1];
                const Real vofneighbor{VolumeFraction(cx, cy)};
                if ((vofneighbor < 0.0 - OFFSET) ||
                    ((int)No_Fluid == Quality(cx, cy) && vofneighbor < 0.1)) {
                    neighborisI2E = true;
                    break;
                }
            }
            if (!(neighborisI2E)) {
            nodeType(0, 0) = VertexType_Interface2Fluid;
            }
        }
    }
#endif  // OPS_2D
}

// those interface cells converted from empty need to be initialized by feq
void KerBuildNewInterface(ACC<Real>& f, ACC<int>& nodeType, ACC<Real>& Rho,
                          const ACC<Real>& U, const ACC<Real>& V,
                          const int* lattIdx,const ACC<Real>& acceleration) {
#ifdef OPS_2D
    const int vt{nodeType(0, 0)};
    if ((int)VertexType::Fluid == vt) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            if (VertexType_Interface2Empty == nodeType(cx, cy)) {
                nodeType(0, 0) = VertexType_Interface;
                break;
            }
        }
    } else if (VertexType_Empty == vt) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            if (VertexType_Interface2Fluid == nodeType(cx, cy)) {
                nodeType(0, 0) = VertexType_NewInterface;
                break;
            }
        }
        if (VertexType_NewInterface == nodeType(0, 0)) {
            int num {0};
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
                    AverRho += Rho(cx, cy);
                    AverU += U(cx, cy);
                    AverV += V(cx, cy);
                    num++;
                }
            }
            AverRho /= num;
            AverU /= num;
            AverV /= num;
            AverV += -TimeStep()*acceleration(1,0,0)/2;
            for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
                f(xiIdx, 0, 0) =
                    CalcBGKFeq(xiIdx, AverRho, AverU, AverV, T, polyOrder);
            }
            Rho(0, 0) = AverRho;  // for KerCalcVolumeFraction
        }                         // local
    }
#endif  // OPS_2D
}

// allocate excess mass by serial computing
// the local nodes are emptied and filled cells
void KerAllocateMass(ACC<Real>& Mass, const ACC<int>& nodeType,
                     const ACC<Real>& VolumeFraction, const ACC<Real>& Rho,
                     const ACC<int>& geometryProperty, const ACC<int>& Quality,
                     const int* lattIdx, const ACC<Real>& coordinates) {
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    if (VertexType_Interface2Fluid == vt || VertexType_Interface2Empty == vt) {
        Real NormalDirection[2]{};
        std::vector<Real> WeightforMass(lattIdx[1] + 1, 0.0);
        std::vector<Real> WeightforMassEven(lattIdx[1] + 1, 0.0);
        Real TotalWeight{0.0};
        Real ExcessMass{0.0};
        int vgneighbors[3][3]{};
        Real vofneighbors[3][3]{};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            vgneighbors[cx + 1][cy + 1] = geometryProperty(cx, cy);
            vofneighbors[cx + 1][cy + 1] = VolumeFraction(cx, cy);
        }
        CalcNormalDirection(NormalDirection, vgneighbors, vofneighbors,
                            lattIdx);
        // calculate normal direction along interface surface,
        if (VertexType_Interface2Fluid == vt) {
            ExcessMass = Mass(0, 0) - Rho(0, 0);
            Mass(0, 0) = Rho(0, 0);
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                int cx = (int)XI[xiIndex * LATTDIM];
                int cy = (int)XI[xiIndex * LATTDIM + 1];
                const int vtneighbor {nodeType(cx, cy)};
                Real NormalDirDotLatticeVel =
                    NormalDirection[0] * cx + NormalDirection[1] * cy;
                if (VertexType_Interface == vtneighbor ||
                    VertexType_NewInterface == vtneighbor) {
                    if (NormalDirDotLatticeVel>1e-8) {
                        WeightforMass[xiIndex] = NormalDirDotLatticeVel;
                    }
                    WeightforMassEven[xiIndex] = 1.0;
                }
                TotalWeight += WeightforMass[xiIndex];
            }
        }
        if (VertexType_Interface2Empty == vt) {
            ExcessMass = Mass(0, 0);
            Mass(0, 0) = 0.0;
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                int cx = (int)XI[xiIndex * LATTDIM];
                int cy = (int)XI[xiIndex * LATTDIM + 1];
                const int vtneighbor {nodeType(cx, cy)};
                Real NormalDirDotLatticeVel =
                    -(NormalDirection[0] * cx + NormalDirection[1] * cy);
                    //here n*c has been multiplyed by negative sign
                if (VertexType_Interface == vtneighbor ||
                    VertexType_NewInterface == vtneighbor) {
                    if (NormalDirDotLatticeVel>1e-8) {
                        WeightforMass[xiIndex] = NormalDirDotLatticeVel;
                    }
                    WeightforMassEven[xiIndex] = 1.0;
                }
                TotalWeight += WeightforMass[xiIndex];
            }
        }
        // EssentiallyEqual(&TotalWeight,&ZERO,EPS)
        if (TotalWeight < 1e-8) {
            WeightforMass.assign(WeightforMassEven.begin(),
                                 WeightforMassEven.end());
            TotalWeight = 0.0;
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                TotalWeight += WeightforMass[xiIndex];
            }
        }
        if (TotalWeight < 1e-8 && ((int)No_Fluid != Quality(0, 0) &&
                                   (int)No_Empty != Quality(0, 0))) {
            ops_printf(
                "Error! TotalWeight =%e becomes invalid! Maybe something "
                "wrong "
                "at x=%f y=%f\n",
                TotalWeight, coordinates(0, 0, 0), coordinates(1, 0, 0));
            assert(!(TotalWeight < 1e-8));
        }
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            if (WeightforMass[xiIndex] > 1e-8) {
                Mass(cx, cy) +=
                    ExcessMass * WeightforMass[xiIndex] / TotalWeight;
            }
        }
    }
#endif  // OPS_2D
}

// calculate the volume of fraction vof = m/rho
void KerCalcVolumeFraction(ACC<Real>& Mass, ACC<Real>& VolumeFraction,
                          ACC<int>& nodeType, const ACC<Real>& Rho ) {
    const int vt{nodeType(0, 0)};
    if (VertexType_Interface == vt || VertexType_NewInterface == vt) {
        nodeType(0, 0) = VertexType_Interface;
        VolumeFraction(0, 0) = Mass(0, 0) / Rho(0, 0);
    } else if ((int)VertexType::Fluid == vt ||
               VertexType_Interface2Fluid == vt) {
        nodeType(0, 0) = (int)VertexType::Fluid;
        Mass(0, 0) = Rho(0, 0);
        VolumeFraction(0, 0) = 1.0;
    } else if (VertexType_Empty == vt || VertexType_Interface2Empty == vt) {
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
                          const ACC<int>& geometryProperty,
                          const int* lattIdx,const ACC<Real>& Rho,
                          const ACC<Real>& U,const ACC<Real>& V) {
#ifdef OPS_2D
    const int vt{nodeType(0, 0)};
    VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
    for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
        int cx = (int)XI[xiIdx * LATTDIM];
        int cy = (int)XI[xiIdx * LATTDIM + 1];
        const int vtneighbor{nodeType(cx, cy)};
        if ((int)VertexType::Fluid == vtneighbor ||
            VertexType_Interface == vtneighbor) {
            if ((int)VertexType::Wall == vt) {
                fStage(xiIdx, 0, 0) = fStage(OPP[xiIdx], cx, cy);
            } else if (VertexType_SlipWall == vt) {
                if (VG_IP == vg || VG_IM == vg) {
                    for (int xiIdx2 = lattIdx[0]; xiIdx2 <= lattIdx[1];
                         xiIdx2++) {
                        int cxx = (int)XI[xiIdx2 * LATTDIM];
                        int cyy = (int)XI[xiIdx2 * LATTDIM + 1];
                        if (cxx == -cx && cyy == cy) {
                            fStage(xiIdx, 0, 0) = fStage(xiIdx2, cx, 0);
                            break;
                        }
                    }
                    if (VertexType_Empty == nodeType(cx, 0)) {
                        // Real u{U(cx, cy)};
                        Real v{V(cx, cy)};
                        Real rho{Rho(cx, cy)};
                        fStage(xiIdx, 0, 0) =
                            fStage(OPP[xiIdx], cx, cy) +
                            +(2 * WEIGHTS[OPP[xiIdx]]) * rho *
                                (0+ (Real)cy * v) * CS;
                    }
                } else if (VG_JM == vg || VG_JP == vg) {
                    for (int xiIdx2 = lattIdx[0]; xiIdx2 <= lattIdx[1];
                         xiIdx2++) {
                        int cxx = (int)XI[xiIdx2 * LATTDIM];
                        int cyy = (int)XI[xiIdx2 * LATTDIM + 1];
                        if (cxx == cx && cyy == -cy) {
                            fStage(xiIdx, 0, 0) = fStage(xiIdx2, 0, cy);
                            break;
                        }
                    }
                    if (VertexType_Empty == nodeType(0, cy)) {
                        Real u{U(cx, cy)};
                        // Real v{V(cx, cy)};
                        Real rho{Rho(cx, cy)};
                        fStage(xiIdx, 0, 0) =
                            fStage(OPP[xiIdx], cx, cy) +
                            +(2 * WEIGHTS[OPP[xiIdx]]) * rho *
                                ((Real)cx * u + 0) * CS;
                    }
                }
            }//Slip boundary condition is temporarily coded here.
        }
    }
#endif
}

//For slip boundary should new create a boundary scheme
void KerHalfwayBounceForward(ACC<Real>& fStage, const ACC<int>& nodeType,
                          const ACC<int>& geometryProperty,
                          const int* lattIdx,const ACC<Real>& Rho,
                          const ACC<Real>& U,const ACC<Real>& V) {
#ifdef OPS_2D
    const int vt{nodeType(0, 0)};
    VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
    for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
        int cx = (int)XI[xiIdx * LATTDIM];
        int cy = (int)XI[xiIdx * LATTDIM + 1];
        const int vtneighbor{nodeType(cx, cy)};
        if ((int)VertexType::Fluid == vtneighbor ||
            VertexType_Interface == vtneighbor) {
            if ((int)VertexType::Wall == vt) {
                fStage(xiIdx, 0, 0) = fStage(OPP[xiIdx], cx, cy);
            } else if (VertexType_SlipWall == vt) {
                if (VG_IP == vg || VG_IM == vg) {
                    for (int xiIdx2 = lattIdx[0]; xiIdx2 <= lattIdx[1];
                         xiIdx2++) {
                        int cxx = (int)XI[xiIdx2 * LATTDIM];
                        int cyy = (int)XI[xiIdx2 * LATTDIM + 1];
                        if (cxx == -cx && cyy == cy) {
                            fStage(xiIdx, 0, 0) = fStage(xiIdx2, cx, 0);
                            break;
                        }
                    }
                    if (VertexType_Empty == nodeType(cx, 0)) {
                        // Real u{U(cx, cy)};
                        Real v{V(cx, cy)};
                        Real rho{Rho(cx, cy)};
                        fStage(xiIdx, 0, 0) =
                            fStage(OPP[xiIdx], cx, cy) +
                            +(2 * WEIGHTS[OPP[xiIdx]]) * rho *
                                (0+ (Real)cy * v) * CS;
                    }
                } else if (VG_JM == vg || VG_JP == vg) {
                    for (int xiIdx2 = lattIdx[0]; xiIdx2 <= lattIdx[1];
                         xiIdx2++) {
                        int cxx = (int)XI[xiIdx2 * LATTDIM];
                        int cyy = (int)XI[xiIdx2 * LATTDIM + 1];
                        if (cxx == cx && cyy == -cy) {
                            fStage(xiIdx, 0, 0) = fStage(xiIdx2, 0, cy);
                            break;
                        }
                    }
                    if (VertexType_Empty == nodeType(0, cy)) {
                        Real u{U(cx, cy)};
                        // Real v{V(cx, cy)};
                        Real rho{Rho(cx, cy)};
                        fStage(xiIdx, 0, 0) =
                            fStage(OPP[xiIdx], cx, cy) +
                            +(2 * WEIGHTS[OPP[xiIdx]]) * rho *
                                ((Real)cx * u + 0) * CS;
                    }
                }
            }//Slip boundary condition is temporarily coded here.
        }
    }
#endif
}

// allocate excess mass for parellel computing
// a local node is always those interface cells are going to lose or obtain mass
// instead filled or emptied
#ifdef FreeMpi_2D
void KerCalcTotalWeight(ACC<Real>& TotalWeightsWrong, ACC<Real>& TotalWeightsRight,
                        ACC<Real>& NormalDirectionsX, ACC<Real>& NormalDirectionsY,
                        const ACC<Real>& VolumeFraction,
                        const ACC<int>& nodeType,
                        const ACC<int>& geometryProperty,
                        const ACC<int>& Quality, const int* lattIdx,
                        const ACC<Real>& coordinates) {
#ifdef OPS_2D
    const int vt{nodeType(0, 0)};
    NormalDirectionsX( 0, 0) = 0.0;
    NormalDirectionsY( 0, 0) = 0.0;
    TotalWeightsWrong( 0, 0) = 0.0;
    TotalWeightsRight( 0, 0) = 0.0;
    if (VertexType_Interface == vt ||
        VertexType_NewInterface == vt) {
        Real NormalDirection[2]{};
        int vgneighbors[3][3]{};
        Real vofneighbors[3][3]{};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            vgneighbors[cx + 1][cy + 1] = geometryProperty(cx, cy);
            vofneighbors[cx + 1][cy + 1] = VolumeFraction(cx, cy);
        }
        CalcNormalDirection(NormalDirection, vgneighbors, vofneighbors,
                            lattIdx);
        NormalDirectionsX(0, 0) = NormalDirection[0];
        NormalDirectionsY(0, 0) = NormalDirection[1];
    }
    if (VertexType_Interface2Fluid == vt || VertexType_Interface2Empty == vt) {
        Real NormalDirection[2]{};
        Real WeightforMass[lattIdx[1] + 1] {0.0};
        Real WeightforMassEven[lattIdx[1] + 1] {0.0};
        Real TotalWeight{0.0};
        int vgneighbors[3][3]{};
        Real vofneighbors[3][3]{};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            vgneighbors[cx + 1][cy + 1] = geometryProperty(cx, cy);
            vofneighbors[cx + 1][cy + 1] = VolumeFraction(cx, cy);
        }
        CalcNormalDirection(NormalDirection, vgneighbors, vofneighbors,
                            lattIdx);
        NormalDirectionsX( 0, 0) = NormalDirection[0];
        NormalDirectionsY( 0, 0) = NormalDirection[1];
        // calculate normal direction along interface surface,
        if (VertexType_Interface2Fluid == vt) {
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                int cx = (int)XI[xiIndex * LATTDIM];
                int cy = (int)XI[xiIndex * LATTDIM + 1];
                const int vtneighbor{nodeType(cx, cy)};
                Real NormalDirDotLatticeVel =
                    NormalDirection[0] * cx + NormalDirection[1] * cy;
                if (VertexType_Interface == vtneighbor ||
                    VertexType_NewInterface == vtneighbor) {
                    if (NormalDirDotLatticeVel > 1e-8) {
                        WeightforMass[xiIndex] = NormalDirDotLatticeVel;//1.0/((Real)cx*cx+(Real)cy*cy);
                    }
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
                // here n*c has been multiplyed by negative sign
                if (VertexType_Interface == vtneighbor ||
                    VertexType_NewInterface == vtneighbor) {
                    if (NormalDirDotLatticeVel < -1e-8) {
                        WeightforMass[xiIndex] = -NormalDirDotLatticeVel;//1.0/((Real)cx*cx+(Real)cy*cy);
                    }
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
                                   (int)No_Empty != Quality(0, 0))) {
            ops_printf(
                "Error! TotalWeight =%e becomes invalid! Maybe updating flag "
                "is wrong "
                "at x=%f y=%f\n",
                TotalWeight, coordinates(0, 0, 0), coordinates(1, 0, 0));
            // Indeed the 'totalweight=0' always reflects wrong flag updating
            assert(!(TotalWeight < 1e-8));
        }
        TotalWeightsRight(0, 0) = TotalWeight;
    }
#endif  // OPS_2D
}

void KerAllocateMassMpi(ACC<Real>& Mass, const ACC<int>& nodeType,
                        const ACC<Real>& Rho, 
                        const ACC<Real>& TotalWeightsWrong,
                        const ACC<Real>& TotalWeightsRight,
                        const ACC<Real>& NormalDirectionsX,
                        const ACC<Real>& NormalDirectionsY, const int* lattIdx) {
#ifdef OPS_2D
    const int vt = nodeType(0, 0);
    if (VertexType_Interface == vt || VertexType_NewInterface == vt) {
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            int cx = (int)XI[xiIndex * LATTDIM];
            int cy = (int)XI[xiIndex * LATTDIM + 1];
            const int vtneighbor = nodeType(cx, cy);
            if (VertexType_Interface2Fluid == vtneighbor ||
                VertexType_Interface2Empty == vtneighbor) {
                Real NormalDirDotLatticeVel =
                    -(NormalDirectionsX(cx, cy) * cx +
                      NormalDirectionsY(cx, cy) * cy);
                // Note:the negative sign from that we stand at interface
                Real Weight{0.0};
                Real TotalWeight = TotalWeightsRight(cx, cy);//TotalWeights( cx, cy)
                Real ExcessMass{0.0};
                if (VertexType_Interface2Fluid == vtneighbor) {
                    if (NormalDirDotLatticeVel > 1e-8) {
                        Weight = NormalDirDotLatticeVel;//1.0/((Real)cx*cx+(Real)cy*cy);
                    } else {
                       if (TotalWeightsWrong(cx, cy) < 1e-8) {
                           Weight = 1.0;
                       }
                    }
                    ExcessMass = Mass(cx, cy) - Rho(cx, cy);
                } else if (VertexType_Interface2Empty == vtneighbor ) {
                    if (NormalDirDotLatticeVel < -1e-8){
                        Weight = -NormalDirDotLatticeVel;//1.0/((Real)cx*cx+(Real)cy*cy);
                    }else{
                        if (TotalWeightsWrong(cx, cy) < 1e-8) {
                            Weight = 1.0;
                        }
                    }
                    ExcessMass = Mass(cx, cy);
                }
                Mass(0, 0) += ExcessMass * Weight/ TotalWeight;
            }
        }
    }
#endif  // OPS_2D
}

#endif //FreeMpi_2D


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

void KerCollideMRTFreeSurface(ACC<Real>& fStage, const ACC<Real>& f,
                              const ACC<Real>& coordinates,
                              const ACC<int>& nodeType, const ACC<Real>& Rho,
                              const ACC<Real>& U, const ACC<Real>& V,
                              const Real* tauRef, const Real* dt,
                              const int* lattIdx,
                              const ACC<Real>& acceleration,const int* idx ){
#ifdef OPS_2D
    const int vt {nodeType(0, 0)};
    // collisionRequired: means if collision is required at boundary
    // e.g., the ZouHe boundary condition explicitly requires collision
    bool collisionRequired = (vt != (int)VertexType::ImmersedSolid);
    if (collisionRequired) {
        const Real rho{Rho(0, 0)};
        const Real u{U(0, 0)};
        const Real v{V(0, 0)};
        const Real Fx{acceleration(0, 0, 0)};
        const Real Fy{acceleration(1, 0, 0)};
        const Real T{1};
        const int polyOrder{2};
        const Real MeshSize = CS * (*dt);
        const Real C{0.18 * 0.18};
        Real NonEquilibriumMoment[4]{0, 0, 0, 0};
        Real NonEquilibriumMoment2{0.0};
        Real tau = (*tauRef);
        Real Mf[lattIdx[1] + 1]{0};
        Real MfStage[lattIdx[1] + 1]{0};
        //Turbulence model
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            Real cx = CS * XI[xiIndex * LATTDIM];
            Real cy = CS * XI[xiIndex * LATTDIM + 1];
            const Real feq{CalcBGKFeq(xiIndex, rho, u, v, T, polyOrder)};
            NonEquilibriumMoment[0] += (cx * cx * (f(xiIndex, 0, 0) - feq));
            NonEquilibriumMoment[1] += (cx * cy * (f(xiIndex, 0, 0) - feq));
            NonEquilibriumMoment[2] += (cy * cx * (f(xiIndex, 0, 0) - feq));
            NonEquilibriumMoment[3] += (cy * cy * (f(xiIndex, 0, 0) - feq));
        }
        for (int i = 0; i < 4; i++) {
            NonEquilibriumMoment2 +=
                2 * NonEquilibriumMoment[i] * NonEquilibriumMoment[i];
        }
        tau = 0.5 *
              (tau - 0.5 * (*dt) +
               sqrt((tau + 0.5 * (*dt)) * (tau + 0.5 * (*dt)) +
                    2 * C * MeshSize * MeshSize * sqrt(NonEquilibriumMoment2)/rho));
        Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
        //moment space
         Real EquilibriumMoment[lattIdx[1]+1]{rho, 
                                        rho * (-2 + u * u + v * v),
                                        rho * (1 - u * u - v * v),
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
        Real RelaxationRate[lattIdx[1]+1]{0, 
                                        1, 
                                        1, 
                                        0, 1, 
                                        0, 1, 
                                        dtOvertauPlusdt, 
                                        dtOvertauPlusdt};
        //moment for f
        for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
            for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                Mf[i] += TransformationM[i][xiIndex] * f(xiIndex, 0, 0);
            }
        }
        // MRT collision
        for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
            if (vt == (int)VertexType::Fluid ||
                vt == (int)VertexType::MDPeriodic ||
                VertexType_Interface == vt || 
                (int)VertexType::Inlet == vt ||
                (int)VertexType::OutLet == vt) {
                MfStage[i] =
                    Mf[i] - RelaxationRate[i] * (Mf[i] - EquilibriumMoment[i]) +
                    (1 - RelaxationRate[i] / 2) * (*dt) * ForceMoment[i];
            } else {
                MfStage[i] =
                    Mf[i] - RelaxationRate[i] * (Mf[i] - EquilibriumMoment[i]);
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
#endif
