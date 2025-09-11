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

/** @brief An example main source code of simulating 2D dambreaking by vof
 *  @author Baoming Guo Jianping Meng
 **/
#include <cmath>
#include <iostream>
#include <ostream>
#include <string>
#include "mplb.h"
#include "ops_seq_v2.h"
#include "dam2d_kernel.inc"
RealField VolumeFraction("vof");
RealField Mass("mass");
IntField Quality("Quality");
RealField Taucopy("tau");
RealField TotalWeightsWrong("TotalWeightsWrong");
RealField TotalWeightsRight("TotalWeightsRight");
RealField NormalDirectionsX("NormalDirectionsX");
RealField NormalDirectionsY("NormalDirectionsY");

// Provide macroscopic initial conditions
void SetInitialMacrosVars() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const int rhoId{compo.macroVars.at(Variable_Rho).id};
            ops_par_loop(
                KerSetInitialMacroVars, "KerSetInitialMacroVars", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(g_MacroVars().at(rhoId).at(blockIdx), 1,
                            LOCALSTENCIL, "Real", OPS_RW),
                ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIdx), 1,
                            LOCALSTENCIL, "Real", OPS_RW),
                ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIdx), 1,
                            LOCALSTENCIL, "Real", OPS_RW),
                ops_arg_dat(Mass.at(blockIdx), 1, LOCALSTENCIL, "Real", OPS_RW),
                ops_arg_dat(VolumeFraction.at(blockIdx), 1, LOCALSTENCIL,
                            "Real", OPS_RW),
                ops_arg_dat(Quality.at(blockIdx), 1, LOCALSTENCIL, "int",
                            OPS_RW),
                ops_arg_dat(g_NodeType().at(compo.id).at(block.ID()), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                            LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_idx(),
                ops_arg_dat(g_MacroBodyforce().at(compo.id).at(blockIdx), SpaceDim(),
                            LOCALSTENCIL, "Real", OPS_READ));
        }
    }
}

// The flag in free surface model is initialized
void SetInitialNodeType() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.BulkRange().begin(), block.BulkRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            ops_par_loop(KerSetInitialNodeType, "KerSetInitialNodeType",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_NodeType().at(compo.id).at(block.ID()),
                                     1, LOCALSTENCIL, "int", OPS_RW),
                         ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                         ops_arg_idx());
        }
    }
}

// due to no slip boundaryscheme, here I define a VertexType_SlipWall 
// to identify those will is implemented BounceForward(a slip scheme)
void SetSlipBoundaryNodeType() {
    for (const auto& boundary : BlockBoundaries()) {
        const Block& block{g_Block().at(boundary.blockIndex)};
        // const int surface{(int)boundary.boundarySurface};
        const int blockIndex{block.ID()};
        std::vector<int> range(2 * SpaceDim());
        range.assign(
            block.BoundarySurfaceRange().at(boundary.boundarySurface).begin(),
            block.BoundarySurfaceRange().at(boundary.boundarySurface).end());
        // for convenience:use existing bounceback to refer to halfwayBB
        ops_par_loop(
            KerSetSlipBoundaryNodeType, "KerSetSlipBoundaryNodeType",
            block.Get(), SpaceDim(), range.data(),
            ops_arg_dat(g_NodeType().at(boundary.componentID).at(blockIndex), 1,
                        ONEPTLATTICESTENCIL, "int", OPS_RW),
            ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(), LOCALSTENCIL,
                        "Real", OPS_READ),
            ops_arg_dat(g_GeometryProperty().at(blockIndex), 1, LOCALSTENCIL,
                        "int", OPS_READ));
    }
}

// if all flags like interface are included into enum VertexType
// it can fo back to original 'UpdateMacroVars'
void UpdateMacroVarsfree() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const Real* pdt{pTimeStep()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            for (auto& macroVar : compo.macroVars) {
                const int varId{macroVar.second.id};
                const VariableTypes varType{macroVar.first};
                switch (varType) {
                    case Variable_Rho:
                        ops_par_loop(
                            KerCalcDensityfree, "KerCalcDensityfree", block.Get(),
                            SpaceDim(), iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    case Variable_U_Force:
                        ops_par_loop(
                            KerCalcUForcefree, "KerCalcUForcefree", block.Get(),
                            SpaceDim(), iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(g_CoordinateXYZ()[blockIndex],
                                        SpaceDim(), LOCALSTENCIL, "double",
                                        OPS_READ),
                            ops_arg_dat(
                                g_MacroBodyforce().at(compo.id).at(blockIndex),
                                SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(pdt, 1, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    case Variable_V_Force:
                        ops_par_loop(
                            KerCalcVForcefree, "KerCalcVForcefree", block.Get(),
                            SpaceDim(), iterRng.data(),
                            ops_arg_dat(g_MacroVars().at(varId).at(blockIndex),
                                        1, LOCALSTENCIL, "double", OPS_RW),
                            ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                        "double", OPS_READ),
                            ops_arg_dat(
                                g_NodeType().at(compo.id).at(blockIndex), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                            ops_arg_dat(g_CoordinateXYZ()[blockIndex],
                                        SpaceDim(), LOCALSTENCIL, "double",
                                        OPS_READ),
                            ops_arg_dat(
                                g_MacroBodyforce().at(compo.id).at(blockIndex),
                                SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_dat(
                                g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, LOCALSTENCIL, "double", OPS_READ),
                            ops_arg_gbl(pdt, 1, "double", OPS_READ),
                            ops_arg_gbl(compo.index, 2, "int", OPS_READ));
                        break;
                    default:
                        break;
                }
            }
        }
    }
#endif  // OPS_2D
}

// Provide macroscopic body-force term
void UpdateMacroscopicBodyForce(const Real time=0.0) {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            // Specify general boundary type
            ops_par_loop(
                KerUpdateMacroBodyForce, "KerUpdateMacroBodyForce", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(g_MacroBodyforce().at(compo.id).at(blockIdx), SpaceDim(),
                            LOCALSTENCIL, "Real", OPS_RW),
                ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                            LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_idx());
        }
    }
}

// similarly, if all flags like interface are included into enum VertexType
void PreDefinedBodyForcefree() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const BodyForceType forceType{compo.bodyForceType};
            switch (forceType) {
                case BodyForce_1st:
                    ops_par_loop(
                        KerCalcBodyForce1STfree, "KerCalcBodyForce1STfree",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_Rho).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_RW),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ));
                    break;
                default:
                    ops_printf(
                        "The specified force type is not implemented!\n");
                    break;
            }
        }
    }
#endif  // OPS_2D
}

// With the function--PreDefinedCollision(),the only difference is
// the present one includes the turbulent model
void FreeSurfaceCollision() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const CollisionType collisionType{compo.collisionType};
            const Real tau{compo.tauRef};
            const Real* pdt{pTimeStep()};
            switch (collisionType) {
                case Collision_BGKIsothermal2nd:
                    ops_par_loop(
                        KerCollideBGKFreeSurface, "KerCollideBGKFreeSurface",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_READ),
                        ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_Rho).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                        ops_arg_dat(Taucopy[blockIndex], 1, LOCALSTENCIL,
                                        "double", OPS_RW));
                    break;
                default:
                    ops_printf(
                        "The specified collision type is not implemented!\n");
                    break;
            }
        }
    }
#endif  // OPS_2D
}

// MRT model
void FreeSurfaceCollisionMRT() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const CollisionType collisionType{compo.collisionType};
            const Real tau{compo.tauRef};
            const Real* pdt{pTimeStep()};
            switch (collisionType) {
                case Collision_BGKIsothermal2nd:
                    ops_par_loop(
                        KerCollideMRTFreeSurface, "KerCollideMRTFreeSurface",
                        block.Get(), SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_READ),
                        ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_Rho).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(&tau, 1, "double", OPS_READ),
                        ops_arg_gbl(pdt, 1, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_idx());
                    break;
                default:
                    ops_printf(
                        "The specified collision type is not implemented!\n");
                    break;
            }
        }
    }
#endif  // OPS_2D
}

// due to the lattice velocity in present LBM is different from the standard.
// the transformation matrix is revised
#ifdef OPS_2D
void MatrixforMRT2D(const int nc9) {
    for (int i = 3; i < nc9; i++) {
        for (int j = 0; j < nc9; j++) {
            if (i < 7) {
                TransformationM[i][j] *= CS;  // D2Q9
                InverseM[j][i] /= CS;
            } else {
                TransformationM[i][j] *= (CS * CS);
                InverseM[j][i] /= (CS * CS);
            }
        }
    }
}
#endif  // OPS_2D

// similarly, if all flags like interface are included into enum VertexType,
// the code can be easily merged into the original 'Stream'
void Streamfree() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};//(blocksizeX - 1)
        for (const auto& compo : g_Components()) {
            ops_par_loop(
                KerStreamfree, "KerStreamfree", block.Get(), SpaceDim(),
                iterRng.data(),
                ops_arg_dat(g_f().at(blockIndex), NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_fStage().at(blockIndex), NUMXI,
                            ONEPTLATTICESTENCIL, "double", OPS_READ),
                ops_arg_dat(g_NodeType().at(compo.first).at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty().at(blockIndex), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_gbl(compo.second.index, 2, "int", OPS_READ));
        }
    }
#endif  // OPS_2D
}

// Bounce back scheme is to revise the fStage here
void ImplementBoundaryFree() {
    for (const auto& boundary : BlockBoundaries()) {
        const Block& block{g_Block().at(boundary.blockIndex)};
        // const int surface{(int)boundary.boundarySurface};
        const int blockIndex{block.ID()};
        const int boundaryScheme{(int)boundary.boundaryScheme};
        Component compo{g_Components().at(boundary.componentID)};
        std::vector<int> range(2 * SpaceDim());
        range.assign(
            block.BoundarySurfaceRange().at(boundary.boundarySurface).begin(),
            block.BoundarySurfaceRange().at(boundary.boundarySurface).end());
        switch (boundaryScheme) {
            case (int)BoundaryScheme::BounceBack: {
                // for convenience:use existing bounceback to refer to halfwayBB
                ops_par_loop(
                    KerHalfwayBounceBack, "KerHalfwayBounceBack", block.Get(),
                    SpaceDim(), range.data(),
                    ops_arg_dat(g_fStage()[blockIndex], NUMXI,
                                ONEPTLATTICESTENCIL, "double", OPS_RW),
                    ops_arg_dat(
                        g_NodeType().at(boundary.componentID).at(blockIndex), 1,
                        ONEPTLATTICESTENCIL, "int", OPS_READ),
                    ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                                LOCALSTENCIL, "int", OPS_READ),
                    ops_arg_gbl(g_Components().at(boundary.componentID).index,
                                2, "int", OPS_READ),
                    ops_arg_dat(g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, ONEPTLATTICESTENCIL, "double", OPS_READ),
                    ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex), 1,
                                ONEPTLATTICESTENCIL, "double", OPS_READ),
                    ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex), 1,
                                ONEPTLATTICESTENCIL, "double", OPS_READ));
            } break;
            case BoundaryScheme_HalfwayBounceForward: {
                ops_par_loop(
                    KerHalfwayBounceForward, "KerHalfwayBounceForward",
                    block.Get(), SpaceDim(), range.data(),
                     ops_arg_dat(g_fStage()[blockIndex], NUMXI,
                                ONEPTLATTICESTENCIL, "double", OPS_RW),
                    ops_arg_dat(
                        g_NodeType().at(boundary.componentID).at(blockIndex), 1,
                        ONEPTLATTICESTENCIL, "int", OPS_READ),
                    ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                                LOCALSTENCIL, "int", OPS_READ),
                    ops_arg_gbl(g_Components().at(boundary.componentID).index,
                                2, "int", OPS_READ),
                    ops_arg_dat(g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIndex),
                                1, ONEPTLATTICESTENCIL, "double", OPS_READ),
                    ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex), 1,
                                ONEPTLATTICESTENCIL, "double", OPS_READ),
                    ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex), 1,
                                ONEPTLATTICESTENCIL, "double", OPS_READ));
            } break;
            default:
                break;
        }
    }
}

void MassExchangeFreeSurface() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            ops_par_loop(UpdatNodeQuality, "UpdatNodeQuality", block.Get(),
                         SpaceDim(), iterRng.data(),
                         ops_arg_dat(Quality.at(blockIndex), 1, LOCALSTENCIL,
                                     "int", OPS_RW),
                         ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex),
                                     1, ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_gbl(compo.index, 2, "int", OPS_READ));

            ops_par_loop(KerMassExchange, "KerMassExchange", block.Get(),
                         SpaceDim(), iterRng.data(),
                         ops_arg_dat(Mass.at(blockIndex), 1, LOCALSTENCIL,
                                     "Real", OPS_RW),
                         ops_arg_dat(g_fStage()[blockIndex], NUMXI,
                                     ONEPTLATTICESTENCIL, "double", OPS_READ),
                         ops_arg_dat(VolumeFraction.at(blockIndex), 1,
                                     ONEPTLATTICESTENCIL, "Real", OPS_READ),
                         ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex),
                                     1, ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_dat(Quality.at(blockIndex), 1,
                                     ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                                     ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                     LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_gbl(compo.index, 2, "int", OPS_READ));

            ops_par_loop(KerReconstructDF, "KerReconstructDF", block.Get(),
                         SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                     "double", OPS_RW),
                         ops_arg_dat(g_fStage()[blockIndex], NUMXI,
                                     LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_dat(VolumeFraction.at(blockIndex), 1,
                                     ONEPTLATTICESTENCIL, "Real", OPS_READ),
                         ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex),
                                     1, ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                                     ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
                                     1, LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
                                     1, LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                     LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                         ops_arg_dat(g_MacroBodyforce().at(compo.id).at(blockIndex), SpaceDim(),
                            LOCALSTENCIL, "Real", OPS_READ));

        }
    }
#endif
}

void UpdateFlagAllocateMass() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            ops_par_loop(KerPreUpdateFlag, "KerPreUpdateFlag", block.Get(),
                         SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex),
                                     1, LOCALSTENCIL, "int", OPS_RW),
                         ops_arg_dat(VolumeFraction.at(blockIndex), 1,
                                     ONEPTLATTICESTENCIL, "Real", OPS_READ),
                         ops_arg_dat(Quality.at(blockIndex), 1,
                                     ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_gbl(compo.index, 2, "int", OPS_READ));                       

            ops_par_loop(
                KerBuildNewInterface, "KerBuildNewInterface", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "int", OPS_RW),
                ops_arg_dat(g_MacroVars()
                                .at(compo.macroVars.at(Variable_Rho).id)
                                .at(blockIndex),
                            1, ONEPTLATTICESTENCIL, "double", OPS_RW),
                ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "double", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "double", OPS_READ),
                ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                ops_arg_dat(g_MacroBodyforce().at(compo.id).at(blockIndex), SpaceDim(),
                            LOCALSTENCIL, "Real", OPS_READ));
#ifndef FreeMpi_2D
            ops_par_loop(
                KerAllocateMass, "KerAllocateMass", block.Get(), SpaceDim(),
                iterRng.data(),
                ops_arg_dat(Mass.at(blockIndex), 1, ONEPTLATTICESTENCIL, "Real",
                            OPS_RW),
                ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "int", OPS_READ),
                ops_arg_dat(VolumeFraction.at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_MacroVars()
                                .at(compo.macroVars.at(Variable_Rho).id)
                                .at(blockIndex),
                            1, LOCALSTENCIL, "double", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                            ONEPTLATTICESTENCIL, "int", OPS_READ),
                ops_arg_dat(Quality.at(blockIndex), 1, LOCALSTENCIL,
                            "int", OPS_READ),
                ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                            LOCALSTENCIL, "double", OPS_READ));
#endif
#ifdef FreeMpi_2D
            ops_par_loop(
                KerCalcTotalWeight, "KerCalcTotalWeight", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(TotalWeightsWrong.at(blockIndex), 1, LOCALSTENCIL,
                            "Real", OPS_RW),
                ops_arg_dat(TotalWeightsRight.at(blockIndex), 1, LOCALSTENCIL,
                            "Real", OPS_RW),
                ops_arg_dat(NormalDirectionsX.at(blockIndex), 1, LOCALSTENCIL,
                            "Real", OPS_RW),
                ops_arg_dat(NormalDirectionsY.at(blockIndex), 1, LOCALSTENCIL,
                            "Real", OPS_RW),
                ops_arg_dat(VolumeFraction.at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                            ONEPTLATTICESTENCIL, "int", OPS_READ),
                ops_arg_dat(Quality.at(blockIndex), 1, LOCALSTENCIL, "int",
                            OPS_READ),
                ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                            LOCALSTENCIL, "double", OPS_READ));

            ops_par_loop(
                KerAllocateMassMpi, "KerAllocateMassMpi", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(Mass.at(blockIndex), 1, ONEPTLATTICESTENCIL, 
                            "Real", OPS_RW),
                ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "int", OPS_READ),
                ops_arg_dat(g_MacroVars()
                                .at(compo.macroVars.at(Variable_Rho).id)
                                .at(blockIndex),
                            1, ONEPTLATTICESTENCIL, "double", OPS_READ),
                ops_arg_dat(TotalWeightsWrong.at(blockIndex), 1, 
                            ONEPTLATTICESTENCIL, "Real", OPS_READ),
                ops_arg_dat(TotalWeightsRight.at(blockIndex), 1, 
                            ONEPTLATTICESTENCIL, "Real", OPS_READ),
                ops_arg_dat(NormalDirectionsX.at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "Real", OPS_READ),
                ops_arg_dat(NormalDirectionsY.at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "Real", OPS_READ),
                ops_arg_gbl(compo.index, 2, "int", OPS_READ));
#endif
        }
    }
#endif
}

void CalcVolumeFraction() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            ops_par_loop(
                KerCalcVolumeFraction, "KerCalcVolumeFraction", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(Mass.at(blockIdx), 1, LOCALSTENCIL, 
                            "Real", OPS_RW),
                ops_arg_dat(VolumeFraction.at(blockIdx), 1, LOCALSTENCIL,
                            "Real", OPS_RW),
                ops_arg_dat(g_NodeType().at(compo.id).at(blockIdx), 1,
                            LOCALSTENCIL, "int", OPS_RW),
                ops_arg_dat(g_MacroVars()
                                .at(compo.macroVars.at(Variable_Rho).id)
                                .at(blockIdx),
                            1, LOCALSTENCIL, "double", OPS_READ));
        }
    }
}

void StreamCollisionFreeSurface(const Real time) {
#if DebugLevel >= 1
    ops_printf("Calculating the macroscopic variables...\n");
#endif
#ifdef OPS_2D
    UpdateMacroVarsfree();
#endif
    // CopyBlockEnvelopDistribution(g_fStage(), g_f());
#if DebugLevel >= 1
    ops_printf("Calculating the mesoscopic body force term...\n");
#endif
    UpdateMacroscopicBodyForce(time);
#ifdef OPS_2D
    PreDefinedBodyForcefree();
#endif
#if DebugLevel >= 1
    ops_printf("Calculating the collision term...\n");
#endif
#ifdef OPS_2D
    // PreDefinedCollision();
    // FreeSurfaceCollision();
    FreeSurfaceCollisionMRT();
#endif

#if DebugLevel >= 1
    ops_printf("Updating the halos...\n");
#endif
    TransferHalos();  // have not read,may need to revise too

#if DebugLevel >= 1
    ops_printf("Streaming...\n");
#endif

#ifdef OPS_2D
    ImplementBoundaryFree();
    // ImplementBoundary();
#endif

#ifdef OPS_2D
    Streamfree();
#endif

#ifdef OPS_2D
    MassExchangeFreeSurface();
#endif
    UpdateMacroVarsfree();
    CalcVolumeFraction();

#ifdef OPS_2D
    UpdateFlagAllocateMass();
    CalcVolumeFraction();
#endif
}


template <typename T>
void Iteratefree(void (*cycle)(T), const SizeType steps,
             const SizeType checkPointPeriod, const SizeType start = 0) {
    ops_printf("Starting the iteration...\n");
    for (SizeType iter = start; iter < start + steps; iter++) {
        const Real time{iter * TimeStep()};
        cycle(time);
        // if (iter>2000&&iter<3000){
        if (((iter + 1) % checkPointPeriod) == 0) {
            ops_printf("%d iterations!\n", iter + 1);
#ifdef OPS_3D
            UpdateMacroVars3D();
#endif
#ifdef OPS_2D
            UpdateMacroVarsfree();
#endif
            WriteFlowfieldToHdf5((iter + 1));
            // WriteDistributionsToHdf5((iter + 1));
            WriteNodePropertyToHdf5((iter + 1));
            VolumeFraction.WriteToHDF5("DamBreaking", iter + 1);
            // Mass.WriteToHDF5("DamBreaking", iter + 1);
            // Quality.WriteToHDF5("DamBreaking", iter + 1);
            // Taucopy.WriteToHDF5("DamBreaking", iter + 1);
            // TotalWeightsWrong.WriteToHDF5("DamBreaking", iter + 1);
            // TotalWeightsRight.WriteToHDF5("DamBreaking", iter + 1);
            NormalDirectionsX.WriteToHDF5("DamBreaking", iter + 1);
            NormalDirectionsY.WriteToHDF5("DamBreaking", iter + 1);
        // }
        }
    }
    ops_printf("Simulation finished! Exiting...\n");
    DestroyModel();
}

void simulate(const Configuration& config) {
    DefineCase(config.caseName, config.spaceDim, config.transient);
    DefineBlocks(config.blockIds, config.blockNames, config.blockSize,
                 config.meshSize, config.startPos);
    VolumeFraction.CreateFieldFromScratch(g_Block());
    Mass.CreateFieldFromScratch(g_Block());
    Quality.CreateFieldFromScratch(g_Block());
    Taucopy.CreateFieldFromScratch(g_Block());
    TotalWeightsWrong.CreateFieldFromScratch(g_Block());
    TotalWeightsRight.CreateFieldFromScratch(g_Block());
    // TotalWeights.SetDataDim(2);
    NormalDirectionsX.CreateFieldFromScratch(g_Block());
    NormalDirectionsY.CreateFieldFromScratch(g_Block());
    // NormalDirections.SetDataDim(SpaceDim());
    DefineComponents(config.compoNames, config.compoIds, config.lattNames,
                     config.tauRef, config.currentTimeStep);
    DefineMacroVars(config.macroVarTypes, config.macroVarNames,
                    config.macroVarIds, config.macroCompoIds,
                    config.currentTimeStep);
    DefineCollision(config.CollisionTypes, config.CollisionCompoIds);
    DefineBodyForce(config.bodyForceTypes, config.bodyForceCompoIds);
    DefineScheme(config.schemeType);
    DefineInitialCondition(config.initialTypes, config.initialConditionCompoId);
    for (auto& bcConfig : config.blockBoundaryConfig) {
        DefineBlockBoundary(bcConfig.blockIndex, bcConfig.componentID,
                            bcConfig.boundarySurface, bcConfig.boundaryScheme,
                            bcConfig.macroVarTypesatBoundary,
                            bcConfig.givenVars, bcConfig.boundaryType);
    }
    Partition();
    ops_diagnostic_output();
    SetTimeStep(config.meshSize / SoundSpeed());
    if (config.currentTimeStep == 0) {
        UpdateMacroscopicBodyForce();
        SetInitialNodeType();
        SetInitialMacrosVars();
        PreDefinedInitialCondition();
        SetSlipBoundaryNodeType();
        MatrixforMRT2D(9);
    };
    // VolumeFraction.WriteToHDF5("DamBreaking", 1);

    if (config.transient) {
        Iteratefree(StreamCollisionFreeSurface, config.timeStepsToRun,
                config.checkPeriod, config.currentTimeStep);
    } else {
        Iteratefree(StreamCollisionFreeSurface, config.convergenceCriteria,
                config.checkPeriod, config.currentTimeStep);
    }
}

int main(int argc, const char** argv) {
    // OPS initialisation where a few arguments can be passed to set
    // the simulation
    ops_init(argc, argv, 4);
    bool configFileFound{false};
    std::string configFileName;
    GetConfigFileFromCmd(configFileFound, configFileName, argc, argv);
    double ct0, ct1, et0, et1;
    ops_timers(&ct0, &et0);
    // start a new simulaton from a configuration file
    if (configFileFound) {
        ReadConfiguration(configFileName);
        simulate(Config());
    }
    ops_timers(&ct1, &et1);
    ops_printf("\nTotal Wall time %lf\n", et1 - et0);
    // Print OPS performance details to output stream
    ops_timing_output(std::cout);
    ops_exit();
}