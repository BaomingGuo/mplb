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
 #include "dam2d_host_device.h"
 #define VertexType_Interface 2000
 #define VertexType_Empty 2001
 #define VertexType_NewInterface 2003
 #define VertexType_Interface2Empty 2004
 #define VertexType_Interface2Fluid 2005
 #define VertexType_SlipWall 2006
 const int iMax {73};
 const int jMax {iMax*2};
 const int Finaldepth{80};
 const int BoundaryScheme_HalfwayBounceBack{100};
 const int BoundaryScheme_HalfwayBounceForward{101};
 int jofInletInterface{jMax};
 Real VOFofInletInterface{0.99};
 const Real Gravity{0};//4.2748e-4
 bool ConsistentInitialize{false};
 #define FreeMpi_2D
 enum VertexQuality { No_Fluid = 1, No_Empty = 2, Standard = 3, Solid = -1 };
 
 
     void KerSetInitialMacroVars(ACC<Real>& rho, ACC<Real>& u, ACC<Real>& v,
                                 ACC<Real>& mass, ACC<Real>& vof, 
                                 const ACC<int>& nodeType,
                                 const int* idx,
                                 const ACC<Real>& acceleration,
                                 ACC<Real>& Porosity) {
         const int vt{nodeType(0, 0)};
         const int i{idx[0]};
         const int j{idx[1]};
         Porosity(0, 0) = 1;
         rho(0, 0) =
             1 + CS * TimeStep() * (idx[1] - jMax- (VOFofInletInterface - 0.5)) * acceleration(1, 0, 0);
         //+MeshSize*(j-Finaldepth)*acceleration(1,0,0);//+(jMax-0.5)*MeshSize*Gravity;
         u(0, 0) = 0;
         v(0, 0) = 0;
         // const int i{idx[0]};
         //  Initialise the mass and vof here according to x,y
         if ((int)VertexType::Fluid == vt) {
             vof(0, 0) = 1.0;
         } else if (VertexType_Interface == vt) {
             if (iMax == i && (j < jMax && j > 0)) {
                 vof(0, 0) = 0.99;
             } else if (jMax == j && (i < iMax && i > 0)) {
                 vof(0, 0) = 0.99;
             } else if (iMax == i && jMax == j) {
                 vof(0, 0) = 0.99;
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
         mass(0, 0) = vof(0, 0) * 1 * Porosity(0, 0);
         // porosity is key for sovling the vof of interface cell
     }
 
 void KerSetInitialNodeType(ACC<int>& nodeType, 
                            const ACC<Real>& coordinates,const int* idx,
                            ACC<int>& geometry) {
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
     if ((i > iMax*2 && i <iMax*2+int(0.024/0.146*iMax)+1) && (j > 0 && j < int(0.048/0.146*iMax)+1)) {
         nodeType(0, 0) = (int)VertexType::Wall;
     }
     if (i==iMax*2+1 && j > 1 && j <int(0.048/0.146*iMax)+1){
        geometry(0,0) = (int)VG_IM;
     }else if(j==int(0.048/0.146*iMax)&&i > iMax*2+1 && i <iMax*2+int(0.024/0.146*iMax)){
        geometry(0,0) = (int)VG_JP;
     }else if(i==iMax*2+1&&j==1){
        geometry(0,0) = (int)VG_IMJP_I;
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
     if ((coordinates(0, 0, 0) > iMax*2 * CS * TimeStep() &&
          coordinates(0, 0, 0) < (iMax*2+0.024/0.146*iMax+0.5 )* CS * TimeStep())) {
         nodeType(0, 0) = (int)VertexType::Wall;
     }
     //the corner is still no slip wall.
 }
 
 void KerInitialiseBGK2ndfree(ACC<Real>& f, ACC<Real>& fStage, const ACC<int>& nodeType,
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
             f(xiIdx, 0, 0) = CalcBGKFeqIncom(xiIdx, rho, u, v, T, polyOrder);
             fStage(xiIdx, 0, 0) = f(xiIdx, 0, 0);
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
     acceleration(1, 0, 0) = -4.2748e-4;//
     acceleration(0, 0, 0) = 0;
 }
 
 void KerUpdateMacroBodyForcetimeVOF(ACC<Real>& acceleration,
                              const ACC<Real>& coordinates, 
                              const int* idx,const ACC<Real>& vof) {
     // here we put some constant for gravity to the following two lines
     acceleration(1, 0, 0) = acceleration(1, 0, 0)*vof(0,0);
     acceleration(0, 0, 0) = acceleration(0, 0, 0)*vof(0,0);
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
         u /= 1;
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
         v /= 1;
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
         Real rho{1};
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
             const Real feq{CalcBGKFeqIncom(xiIndex, rho, u, v, T, polyOrder)};
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
 
 void KerCollideMRTFreeSurface(ACC<Real>& fStage, const ACC<Real>& f,
                               const ACC<Real>& coordinates,
                               const ACC<int>& nodeType, const ACC<Real>& Rho,
                               const ACC<Real>& U, const ACC<Real>& V,
                               const Real* tauRef, const Real* dt,
                               const int* lattIdx,
                               const ACC<Real>& acceleration,const int* idx,
                               const Real* endCoordinate, ACC<Real>&Tau){
 #ifdef OPS_2D
     // const Real iCoordinate{coordinates(1,0,0)};
     const int vt{nodeType(0, 0)};
     // collisionRequired: means if collision is required at boundary
     // e.g., the ZouHe boundary condition explicitly requires collision
     bool collisionRequired = (vt != (int)VertexType::ImmersedSolid&& 
                               VertexType_Empty != vt);
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
         tau = Smagorinsky_LES(rho, u, v, tau, fstream, lattIdx,Fx,Fy);
         Tau(0, 0) = tau;
         Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
         // printf("\nTotal Wall time l=%e\n", dtOvertauPlusdt);
         int NumofEmpty{0};
         for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
             int cx = (int)XI[xiIndex * LATTDIM];
             int cy = (int)XI[xiIndex * LATTDIM + 1];
             const int vtneighbor {nodeType(cx, cy)};
             if(VertexType_Empty==vtneighbor){
                 NumofEmpty++;
             }
         }
         // relaxation moment
         std::vector<Real> RelaxationRate;
         if (VertexType_Interface == vt&&NumofEmpty>=1){
             RelaxationRate={1, 
                                         1, 
                                         1, 
                                         1, 1, 
                                         1, 1, 
                                         dtOvertauPlusdt, 
                                         dtOvertauPlusdt};
         }else{
             RelaxationRate={1, 
                                         1, 
                                         1, 
                                         1, 1, 
                                         1, 1, 
                                         dtOvertauPlusdt, 
                                         dtOvertauPlusdt};
         }
 
         
         //moment space
         //moment for f
         for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
             for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                 Mf[i] += TransformationM[i][xiIndex] * f(xiIndex, 0, 0);
             }
         }
          Real EquilibriumMoment[lattIdx[1]+1]{rho, 
                                         (-2*rho + u * u + v * v),
                                         (1*rho - u * u - v * v),
                                         u, -u,
                                         v, -v,
                                         (u * u - v * v), 
                                         u * v};
          Real ForceMoment[lattIdx[1]+1]   {0,   
                                     2 * (Fx * u + Fy * v), 
                                     (-2 * (Fx * u + Fy * v)), 
                                     Fx, -Fx, 
                                     Fy, -Fy, 
                                     2 * (Fx * u - Fy * v), 
                                     (Fy * u + Fx * v)};
         // MRT collision
         for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
             if (vt == (int)VertexType::Fluid ||
                 vt == (int)VertexType::MDPeriodic ||
                 VertexType_Interface == vt || 
                 (int)VertexType::Inlet == vt ||
                 (int)VertexType::OutLet == vt) {
                 MfStage[i] =
                     Mf[i] - RelaxationRate[i] * 
                     (Mf[i] - EquilibriumMoment[i]) +
                     (1 - RelaxationRate[i] / 2) * (*dt) * ForceMoment[i];
             } else {
                 MfStage[i] =
                     Mf[i] - RelaxationRate[i] * 
                     (Mf[i] - EquilibriumMoment[i]);
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
                     res, vt, rho, u, v, coordinates(0, 0, 0),
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
     Quality(0,0) = (int)Standard;
     if (VertexType_Interface == vt) {
         int NumofFluid{0};
         int NumofEmpty{0};
         int NumofBoundary{0};
         for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
             int cx = (int)XI[xiIndex * LATTDIM];
             int cy = (int)XI[xiIndex * LATTDIM + 1];
             const int vtneighbor=nodeType(cx, cy);
             if ((int)VertexType::Fluid == vtneighbor) {
                 NumofFluid++;
             } else if (VertexType_Empty == vtneighbor) {
                 NumofEmpty++;
             } else if ((int)VertexType::Wall == vtneighbor ||
                        VertexType_SlipWall == vtneighbor) {
                 NumofBoundary++;
             }
         }
         //determine the quality of the interface cells:
         if ((0==NumofFluid)&&NumofEmpty>4) {
             Quality(0,0) = (int)No_Fluid ;
         } else if ((0==NumofEmpty)&&NumofFluid>4) {
             Quality(0,0)=(int)No_Empty;
         } else {Quality(0,0)=(int)Standard;}
         //  if (NumofBoundary>0){
         //      if ((0==NumofFluid)&&NumofEmpty>4) {
         //      Quality(0,0) = (int)No_Fluid ;
         //      } else if ((0==NumofEmpty)&&NumofFluid>4) {
         //      Quality(0,0)=(int)No_Empty;
         //     } else {Quality(0,0)=(int)Standard;}
         //  }
         
     }
 #endif  // OPS_2D    
 }
 
 void KerartificialSpatterBubble(ACC<int>& nodeType,  
                       const int* lattIdx) {
 #ifdef OPS_2D
     const int vt {nodeType(0, 0)};
     if (VertexType_Interface == vt) {
         int NumofFluid{0};
         int NumofEmpty{0};
         int NumofBoundary{0};
         for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
             int cx = (int)XI[xiIndex * LATTDIM];
             int cy = (int)XI[xiIndex * LATTDIM + 1];
             const int vtneighbor=nodeType(cx, cy);
             if ((int)VertexType::Fluid == vtneighbor) {
                 NumofFluid++;
             } else if (VertexType_Empty == vtneighbor) {
                 NumofEmpty++;
             } 
             else if ((int)VertexType::Wall == vtneighbor ||
                        VertexType_SlipWall == vtneighbor) {
                 NumofBoundary++;
             }
         }
         if (NumofBoundary+NumofEmpty==NUMXI-1) {
             nodeType(0, 0) = VertexType_Empty;
         } else if (NumofBoundary+NumofFluid == NUMXI-1) {
             nodeType(0, 0) =(int)VertexType::Fluid;
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
                      const ACC<Real>& NormalDirections,
                      const ACC<Real>& rho) {
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
             const Real vofneighbor{VolumeFraction(cx, cy)};
             Real Vofface = vof;
             Real Voffaceneighbor = vofneighbor; 
             // if  (vtneighbor == VertexType_Interface&&abs(cx)+abs(cy)==1){
             // const Real NormalDirection[2]
             //        {NormalDirections(0,0,0),NormalDirections(1,0,0)};
             // const Real NormalDirectionneighbor[2]
             //        {NormalDirections(0,cx,cy),NormalDirections(1,cx,cy)};
             // Vofface = CalcIntersectFree(NormalDirection, vof, cx, cy);
             // Voffaceneighbor = CalcIntersectFree(NormalDirectionneighbor,
             //                                    vofneighbor, -cx, -cy);
             // }
             const Real Area = 0.5 * (Vofface * porosity(0, 0) +
                                 Voffaceneighbor * porosity(cx, cy));
             switch (vtneighbor) {
                 case VertexType_Interface: {
                     VertexQuality vqNeighbor = (VertexQuality)Quality(cx, cy);
                     // Force non-physical (or undesired) spatter and bubble
                     // to fill or empty first step:
                     if ((Standard == vq && Standard == vqNeighbor) ||
                         (No_Fluid == vq && No_Fluid == vqNeighbor) ||
                         (No_Empty == vq && No_Empty == vqNeighbor)) {
                             
                         deltam = Area *
                                  (fStage(OPP[xiIndex], cx, cy)-fStage(xiIndex, 0, 0)
                                     - WEIGHTS[xiIndex]*(rho(cx,cy)-rho(0,0)));
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
                         (fStage(OPP[xiIndex], cx, cy) - fStage(xiIndex, 0, 0)
                         - WEIGHTS[xiIndex]*(rho(cx,cy)-rho(0,0)));
                     break;
                 }
                 case VertexType_Empty: {
                     deltam = 0.0;
                     break;
                 }
                 case (int)VertexType::ImmersedSolid: {
                     deltam = 0.0;//think no mass transfer
                     break;
                 }
                 case VertexType_SlipWall:{
                     const int vgneighbor = geometryProperty(cx, cy);
                     // if (VG_IP == vgneighbor || VG_IM == vgneighbor) {
                     //     if (VertexType_Empty == nodeType(0, cy)) {
                     //         deltam = 0;
                     //     } else if (VertexType_Interface == nodeType(0, cy)){
                     //         deltam = (vof +VolumeFraction(0,cy))*0.5 *porosity(0,0) * 
                     //                      (fStage(OPP[xiIndex], cx, cy) -
                     //                       fStage(xiIndex, 0, 0));
                     //     } else {
                     //         deltam = porosity(0,0) * 
                     //                      (fStage(OPP[xiIndex], cx, cy) -
                     //                       fStage(xiIndex, 0, 0));
                     //     }
                     // } else if (VG_JP == vgneighbor || VG_JM == vgneighbor) {
                     //     if (VertexType_Empty == nodeType(cx, 0)) {
                     //         deltam = 0;
                     //     } else if (VertexType_Interface == nodeType(cx, 0)){
                     //         deltam = (vof +VolumeFraction(cx,0))*0.5*porosity(0,0) * 
                     //                      (fStage(OPP[xiIndex], cx, cy) -
                     //                       fStage(xiIndex, 0, 0));
                     //         // In fact, interface and fluid need to be addressed
                     //         // separately
                     //     }else {
                     //         deltam = porosity(0,0) * 
                     //                      (fStage(OPP[xiIndex], cx, cy) -
                     //                       fStage(xiIndex, 0, 0));
                     //     }
                     // } else {
                         deltam = 0.0;
                     // }
                     break;
                 }
                 case (int)VertexType::Inlet: {
                     if (j + cy < *jofInletInterface) {
                         deltam = 0.5 * (porosity(0, 0) + porosity(cx, cy))*
                                  (fStage(OPP[xiIndex], cx, cy) -
                                   fStage(xiIndex, 0, 0));
                     } else if (j + cy == *jofInletInterface) {
                         deltam = Area *
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
                       const Real* endCoordinate,
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
         // reconstruct the missing DFs from empty cells
         // Warning:do not reconstruct DFs n*c>0 
         int NumofEmpty{0};
         for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
             int cx = (int)XI[xiIndex * LATTDIM];
             int cy = (int)XI[xiIndex * LATTDIM + 1];
             const int vtneighbor {nodeType(cx, cy)};
             if(VertexType_Empty==vtneighbor){
                 NumofEmpty++;
             }
         }
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
                  ((VertexType_Interface == vtneighbor
                     // VertexType_SlipWall == vtneighbor||
                     // (int)VertexType::Wall == vtneighbor
                     ) &&
                     NormalDirections(0,0,0)*cx+NormalDirections(1,0,0)*cy>
                     cos(0.0/180*PI)*sqrt(cx*cx+cy*cy))||
                 ((int)VertexType::OutLet == vtneighbor &&
                  VertexType_Empty == nodeType(0, cy))) {
                 Real u = U(0, 0);
                 Real v = V(0, 0);
                 Real IntersectFree[]{-100, -100};
                 // Real NormalDirection[SpaceDim()]{NormalDirections(0, 0, 0),
                 //                                  NormalDirections(1, 0, 0)};
                 // CalcIntersectFree(NormalDirection, VolumeFraction(0, 0), 
                 //                 cx, cy, coordinates(0,0,0), coordinates(1,0,0), IntersectFree);
                 // if(abs(IntersectFree[1]-j)>0.5){
                     Real ny = NormalDirections(1,0,0);
                     IntersectFree[1] = (j+(abs(ny) <= 0 ? 0 : ny/abs(ny))*
                                             (VolumeFraction(0, 0) - 0.5));
                 // }
                 // if ((int)VertexType::Fluid == nodeType(-cx, -cy) ) {
                     // if (VolumeFraction(0, 0) >= 0.5) {
                     //     &&
                     //  VolumeFraction(-cx, -cy) >= 0.5
                         Real Distance =  VolumeFraction(0, 0) - 0.5;
                         u = U(0, 0) ;
                         v = V(0, 0) ;
                         int CX{0};
                         int CY{0};
                         if (abs(NormalDirections(0,0,0))>=abs(NormalDirections(1,0,0))){
                             if (Distance>0){
                                 CX = 0;
                                 CY = cy;
                             }else{
                                 CX = 0;
                                 CY = -cy;
                             }
                         }else {
                             if (Distance>0){
                                 CX = cx;
                                 CY = 0;
                             }else{
                                 CX = -cx;
                                 CY = 0;
                             }
                         }
                         if (((int)VertexType::Fluid == nodeType(CX,CY) ||
                                         VertexType_Interface == nodeType(CX,CY))){
                             u = U(0, 0)*(1-abs(Distance))+U(CX,CY)*abs(Distance);
                             v = V(0, 0)*(1-abs(Distance))+V(CX,CY)*abs(Distance);
                         }
                         // if (Distance>0&&((int)VertexType::Fluid == nodeType(cx,0) ||
                         //                         VertexType_Interface == nodeType(cx,0))){
                         //     u = U(0, 0)*(1-Distance)+U(cx, 0)*Distance;
                         //     v = V(0, 0)*(1-Distance)+V(cx, 0)*Distance;
                         // }else if(Distance<0&&((int)VertexType::Fluid == nodeType(-cx,0) ||
                         //                         VertexType_Interface == nodeType(-cx,0))){
                         //     u = U(0, 0)*(1+Distance)-U(-cx, 0)*Distance;
                         //     v = V(0, 0)*(1+Distance)-V(-cx, 0)*Distance;
                         // }
                     // } //else {
                     //     u = U(-cx, -cy);
                     //     v = V(-cx, -cy);
                     // }
                 // }
                 rhoA = 1 + (IntersectFree[1] - (jMax+VOFofInletInterface-0.5) ) * CS *
                               TimeStep() * Gravity;
                 // Young-Laplace pressure
                 // Real Curvature{-3.0*4.0/27.0*(VolumeFraction(0, 0)+1)};
                 // Real Tension{0.0728/(0.1*0.146*(57.88*57.88))};
                 //  for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                 //     int cx = (int)XI[xiIndex * LATTDIM];
                 //     int cy = (int)XI[xiIndex * LATTDIM + 1];
                 //     Curvature += VolumeFraction(cx, cy)*4.0/27.0;
                 // }
                 // rhoA -=2* Tension *Curvature;
                 // rhoA =
                 // 1+(j-0.5+(VolumeFraction(0,0)-0.5))*CS*TimeStep()*Gravity;
                 Real feq = CalcBGKFeqIncom(xiIndex, rhoA, u, v, T, polyOrder);
                 Real feqopp =
                      CalcBGKFeqIncom(OPP[xiIndex], rhoA, u, v, T, polyOrder);
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
         if ((vof < 0.0 - OFFSET) || ((int)No_Fluid == vq && vof < 0.1)) {
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
                     ((int)No_Fluid == Quality(cx, cy) && vofneighbor < 0.1)) {
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
                         const int* lattIdx,const ACC<Real>& coordinates,
                         const ACC<Real>& U, const ACC<Real>& V) {
 #ifdef OPS_2D
     const int vt{nodeType(0, 0)};
     // NormalDirections(0, 0, 0) = 0.0;
     // NormalDirections(1, 0, 0) = 0.0;
     TotalWeightsWrong( 0, 0) = 0.0;
     TotalWeightsRight( 0, 0) = 0.0;
     if (VertexType_Interface2Fluid == vt || VertexType_Interface2Empty == vt) {
         const int j{idx[1]};
         Real NormalDirection[2]{};
         std::vector<Real> WeightforMass(lattIdx[1] + 1, 0.0);
         std::vector<Real> WeightforMassEven(lattIdx[1] + 1, 0.0);
         Real TotalWeight{0.0};
         int vgneighbors[3][3]{};
         Real vofneighbors[3][3]{};
         bool neighborisImmerseSolid{false};
         for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
             int cx = (int)XI[xiIndex * LATTDIM];
             int cy = (int)XI[xiIndex * LATTDIM + 1];
             vgneighbors[cx + 1][cy + 1] = geometryProperty(cx, cy);
             vofneighbors[cx + 1][cy + 1] = VolumeFraction(cx, cy);
             if ((int)VertexType::Wall == nodeType(cx, cy)) {
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
                         const ACC<Real>& f,
                         const ACC<Real>& U, const ACC<Real>& V) {
 #ifdef OPS_2D
     const int vt = nodeType(0, 0);
     if (VertexType_Interface == vt || VertexType_NewInterface == vt) {
         for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
             int cx = (int)XI[xiIndex * LATTDIM];
             int cy = (int)XI[xiIndex * LATTDIM + 1];
             const int vtneighbor = nodeType(cx, cy);
             const Real vofneighbor = VolumeFraction(cx, cy);
             if ((VertexType_Interface2Fluid == vtneighbor && vofneighbor > 1) ||
                 (VertexType_Interface2Empty == vtneighbor && vofneighbor < 0)) {
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
                     ExcessMass = Mass(cx, cy) - porosity(cx, cy) * 1;
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
         VolumeFraction(0, 0) = Mass(0, 0) / (porosity(0, 0) * 1);
     } else if ((int)VertexType::Fluid == vt ||
                VertexType_Interface2Fluid == vt) {
         nodeType(0, 0) = (int)VertexType::Fluid;
         Mass(0, 0) = porosity(0,0)* 1;
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
                             if (VertexType_Empty == nodeType(cx, 0)) {
                             // Real u{U(cx, cy)};
                                 Real v{V(cx, cy)};
                                 Real rho{1};
                                 Real weight{1};
                                 fStage(xiIdx, 0, 0) = 
                                     (1-weight)*fStage(xiIdx2, cx, cy)+
                                         weight*(fStage(OPP[xiIdx], cx, cy) +
                                         +(2 * WEIGHTS[OPP[xiIdx]]) * rho *
                                             (0+ (Real)cy * v) * CS);
                             }
                             break;
                         }
                     }
 
                 } else if (VG_JM == vg || VG_JP == vg) {
                     for (int xiIdx2 = lattIdx[0]; xiIdx2 <= lattIdx[1];
                          xiIdx2++) {
                         int cxx = (int)XI[xiIdx2 * LATTDIM];
                         int cyy = (int)XI[xiIdx2 * LATTDIM + 1];
                         if (cxx == cx && cyy == -cy) {
                             fStage(xiIdx, 0, 0) = fStage(xiIdx2, 0, cy);
                              if (VertexType_Empty == nodeType(0, cy)) {
                                 Real u{U(cx, cy)};
                                 // Real v{V(cx, cy)};
                                 Real rho{1};
                                 Real weight{1};
                                 fStage(xiIdx, 0, 0) =
                                     (1-weight)*fStage(xiIdx2, cx, cy)+
                                         weight*(fStage(OPP[xiIdx], cx, cy) +
                                         +(2 * WEIGHTS[OPP[xiIdx]]) * rho *
                                             ((Real)cx * u + 0) * CS);
                             }
                             break;
                         }
                     }
                    
                 }
             }//Slip boundary condition is temporarily coded here.
         }
     }
 #endif
 }
 
 void KerUnitInterBounceBack(ACC<Real>& f, const ACC<Real>& fStage,
                             const ACC<int>& nodeType, const int* lattIdx,
                             const ACC<Real>& Rho,
                             const ACC<Real>& U, const ACC<Real>& V,
                             const ACC<int>& geometryProperty) {
 #ifdef OPS_2D
     // VertexGeometryType vg = (VertexGeometryType)geometryProperty(0, 0);
     const int vt{nodeType(0, 0)};
     Real rho{Rho(0, 0)};
     Real Distance{0.5};
     if ((int)VertexType::Fluid == vt || (VertexType_Interface == vt)) {
         for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
             int cx = (int)XI[xiIdx * LATTDIM];
             int cy = (int)XI[xiIdx * LATTDIM + 1];
             const int vtneighbor{nodeType(cx, cy)};
             Real u{0};
             Real v{0};
             // if ((int)VertexType::Wall == vtneighbor){
             //     f(OPP[xiIdx], 0, 0) =
             //             fStage(xiIdx, 0, 0) ;
             // } else 
             if ((int)VertexType::Wall == vtneighbor|| 
                 VertexType_SlipWall == vtneighbor
                 ) {
                 if ((int)VertexType::Wall == vtneighbor) {
                     u = 0;
                     v = 0;
                 } else if (VertexType_SlipWall == vtneighbor) {
                     VertexGeometryType vg = (VertexGeometryType)geometryProperty(cx, cy);
                     if (VG_JM == vg || VG_JP == vg) {
                         u = U(cx, 0)*Distance+U(0, 0)*(1-Distance);
                         v = 0;
                         if (VertexType_Empty == nodeType(cx, 0)) {
                            u = U(0, 0);
                         } 
                     } else if (VG_IP == vg || VG_IM == vg) {
                         u = 0;
                         v = V(0, cy)*Distance+V(0, 0)*(1-Distance);
                         if (VertexType_Empty == nodeType(0, cy)) {
                            v = V(0, 0);
                         }
                     }
                 }
                 if ((int)VertexType::Fluid != nodeType(-cx, -cy) &&
                     VertexType_Interface != nodeType(-cx, -cy)) {
                     f(OPP[xiIdx], 0, 0) =
                         (Distance * fStage(OPP[xiIdx], 0, 0) +
                          Distance * fStage(xiIdx, 0, 0) +
                          (1 - Distance) * fStage(xiIdx, 0, 0) +
                          2 * WEIGHTS[xiIdx] * rho * (-cx * u - cy * v) * CS) /
                         (Distance + 1);
                 } else {
                     f(OPP[xiIdx], 0, 0) =
                         (Distance * fStage(OPP[xiIdx], 0, 0) +
                          Distance * fStage(xiIdx, 0, 0) +
                          (1 - Distance) * fStage(xiIdx, -cx, -cy) +
                          2 * WEIGHTS[xiIdx] * rho * (-cx * u - cy * v) * CS) /
                         (Distance + 1);
                 }
             }
         }
     }
 
 #endif
 }
 
 void KerInitNewInterface( ACC<Real>& f,
                           ACC<Real>& fStage,
                           ACC<Real>& Rho,
                           ACC<Real>& U, ACC<Real>& V,
                           const ACC<int>& nodeType,
                           const int* lattIdx,
                           const Real* endCoordinate,
                           const Real* tauRef,
                           const ACC<Real>& coordinates,
                           const ACC<Real>& acceleration) {
 #ifdef OPS_2D
     const int vt{nodeType(0, 0)};
     if (VertexType_NewInterface == vt) {
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
         AverU /= TotalWeight;
         AverV /= TotalWeight;
         for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
             fStage(xiIdx, 0, 0) =
                 CalcBGKFeqIncom(xiIdx, AverRho, AverU, AverV, T, polyOrder);
         }
         AverRho = 0;
         AverU = 0;
         AverV = 0;
         TotalWeight = 0;
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
                     u = u / 1 + (TimeStep() * acceleration(0, cx, cy) / 2);
                     v = v / 1 + (TimeStep() * acceleration(1, cx, cy) / 2);
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
         AverU /= TotalWeight;
         AverV /= TotalWeight;
         for (int xiIdx = lattIdx[0]; xiIdx <= lattIdx[1]; xiIdx++) {
             f(xiIdx, 0, 0) =
                 CalcBGKFeqIncom(xiIdx, AverRho, AverU, AverV, T, polyOrder);
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
         (VertexType_NewInterface == vt)};
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
         tau = Smagorinsky_LES(rho, u, v, tau, fstream, lattIdx,Fx,Fy);
         Real dtOvertauPlusdt = (*dt) / (tau + 0.5 * (*dt));
         // relaxation moment
         Real RelaxationRate[lattIdx[1]+1]{1, 
                                         1, 
                                         1, 
                                         1, 1, 
                                         1, 1, 
                                         dtOvertauPlusdt, 
                                         dtOvertauPlusdt};
 
         
         //moment space
         //moment for f
         for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
             for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
                 Mf[i] += TransformationM[i][xiIndex] * f(xiIndex, 0, 0);
             }
         }
         Real EquilibriumMoment[lattIdx[1]+1]{rho, 
                                         1 * (-2*rho + u * u + v * v),
                                         1 * (1*rho - u * u - v * v),
                                         1 * u, -1 * u,
                                         1 * v, -1 * v,
                                         1 * (u * u - v * v), 
                                         1 * u * v};
         Real ForceMoment[lattIdx[1]+1]   {0,   
                                     1 * 2 * (Fx * u + Fy * v), 
                                     1 * (-2 * (Fx * u + Fy * v)), 
                                     1 * Fx, -1 * Fx, 
                                     1 * Fy, -1 * Fy, 
                                     1 * 2 * (Fx * u - Fy * v), 
                                     1 * (Fy * u + Fx * v)};
         // MRT collision
         for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
             MfStage[i] = Mf[i] -
                          RelaxationRate[i] * (Mf[i] - EquilibriumMoment[i]) +
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
                      const ACC<Real>& NormalDirections,
                      const ACC<Real>& Mass) {
 #ifdef OPS_2D
     const int vt{nodeType(0, 0)};
     bool StreamRequired{
         (VertexType_NewInterface == vt)};
     if ((int)VertexType::Fluid == vt || VertexType_Interface == vt||
         VertexType_Interface2Fluid == vt) {
         for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
             int cx = (int)XI[xiIndex * LATTDIM];
             int cy = (int)XI[xiIndex * LATTDIM + 1];
             const int vtneighbor{nodeType(-cx, -cy)};
             if (VertexType_NewInterface == vtneighbor) {
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
         Real vof{Mass(0,0)/1};
         const Real T{1};
         const int polyOrder{2};
         for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
             int cx = (int)XI[xiIndex * LATTDIM];
             int cy = (int)XI[xiIndex * LATTDIM + 1];
             const int vtneighbor{nodeType(-cx, -cy)};
            
             // boundary condition for Free surface FSK
             // and IBB for ImmerseSolid no-slip wall
             if (VertexType_Empty == vtneighbor||
                 VertexType_Interface2Empty == vtneighbor) {
                 //     Real IntersectFree[]{-100, -100};
                 //     Real ny = NormalDirections(1,0,0);
                 //     IntersectFree[1] = (idx[1]+(abs(ny) <= 0 ? 0 : ny/abs(ny))*
                 //                             (vof - 0.5));
 
                 // Real rhoA =
                 //     1 + (IntersectFree[1] - (jMax+VOFofInletInterface-0.5) ) 
                 //         * CS * TimeStep() * Gravity;
                 // // if ((int)VertexPosition::In_Chamber == Position(0, 0)) {
                 // //     rhoA += DeltaRho;
                 // // }
                 // if ((int)VertexType::Fluid == nodeType(cx, cy) ||
                 //     (VertexType_Interface == nodeType(cx, cy))||
                 //     VertexType_Interface2Fluid == nodeType(cx, cy)) {
                 //     Real Distance = vof - 0.5;
                 //     u = U(0, 0) * (1 + Distance) - U(cx, cy) * Distance;
                 //     v = V(0, 0) * (1 + Distance) - V(cx, cy) * Distance;
                 // } else {
                 //     u = U(0, 0);
                 //     v = V(0, 0);
                 // }
                 // Real feq = CalcBGKFeqIncom(xiIndex, rhoA, u, v, T, polyOrder);
                 // Real feqopp =
                 //     CalcBGKFeqIncom(OPP[xiIndex], rhoA, u, v, T, polyOrder);
                 // f(xiIndex, 0, 0) = feq + feqopp - fStage(OPP[xiIndex], 0, 0);
                 f(xiIndex, 0, 0) = fStage(xiIndex, 0, 0);
             }else if (VertexType_SlipWall == vtneighbor){
                 // f(xiIndex, 0, 0) = fStage(OPP[xiIndex], 0, 0)+
                 //         2 * WEIGHTS[xiIndex] * 1 * (cx * U(0, 0) +cy * V(0, 0)) * CS;
                 f(xiIndex, 0, 0) = fStage(xiIndex, 0, 0);
             } else if ((int)VertexType::Wall == vtneighbor){
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
                         ACC<Real>& ResidualError) {
 #ifdef OPS_2D
     const int vt{nodeType(0, 0)};
     ResidualError(0, 0) =0;
     Real vof{0};
     bool StreamRequired{
         ( VertexType_NewInterface == vt)};
     if ((int)VertexType::Fluid == vt || VertexType_Interface == vt||
         VertexType_Interface2Fluid == vt) {
         for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
             int cx = (int)XI[xiIndex * LATTDIM];
             int cy = (int)XI[xiIndex * LATTDIM + 1];
             const int vtneighbor{nodeType(cx, cy)};
             if (
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
             // u += cx * f(xiIdx, 0, 0);
             // v += cy * f(xiIdx, 0, 0);
         }
         if (
             VertexType_NewInterface == vt) {
             ResidualError(0, 0) = 
                         sqrt((rho - Rho(0, 0)) * (rho - Rho(0, 0))) /
                           sqrt((rho + Rho(0, 0)) * (rho + Rho(0, 0)));
             // ops_printf("Error= %e\n", ResidualError(0, 0));
         }
         Rho(0, 0) = rho;
         // U(0, 0) = u / 1 + (TimeStep() * acceleration(0, 0, 0) / 2);
         // V(0, 0) = v / 1 + (TimeStep() * acceleration(1, 0, 0) / 2);
     }
     // if (VertexType_NewInterface == vt) {
     //     int NumofNeighLiquid{0};
     //     if (VertexType_Interface == nodeType(-1, 0) ||
     //         (int)VertexType::Fluid == nodeType(-1, 0)) {
     //         vof += VolumeFraction(-1, 0);
     //         NumofNeighLiquid++;
     //     }
     //     if (VertexType_Interface == nodeType(1, 0) ||
     //         (int)VertexType::Fluid == nodeType(1, 0)) {
     //         vof += VolumeFraction(1, 0);
     //         NumofNeighLiquid++;
     //     }
     //     if (NumofNeighLiquid != 0) {
     //         vof /= NumofNeighLiquid;
     //         VolumeFraction(0, 0) = vof;
     //     } else {
     //         VolumeFraction(0, 0) = 0;
     //     }
     //     Mass(0, 0) = VolumeFraction(0, 0) * 
     //                     Porosity(0, 0) * Rho(0, 0);
     // }
 #endif  // OPS_2D
 }
 
 void KerTotalDensityError(const ACC<int>& nodeType, 
                         const ACC<Real>& ResidualError,
                         double* Error) {
     // const int vt{nodeType(0, 0)};
     // if (VertexType_Solid2Fluid == vt || 
     //     VertexType_NewInterface == vt) {
         *Error += ResidualError(0, 0);
     // }
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
                 Real feq = CalcBGKFeqIncom(xiIdx, rho, u, v, T, polyOrder);
                 Real feqOutgoing = CalcBGKFeqIncom(xiIdx, rhoOutgoing, uOutgoing,
                                               vOutgoing, T, polyOrder);
                 f(xiIdx, 0, 0) =
                     feq - feqOutgoing + f(xiIdx, outgoingX, outgoingY);
             }
         }
     }
 #endif
 }
 
 
 #endif
 