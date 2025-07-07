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
 *  @author Jianping Meng
 **/
#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include "mplb.h"
#include "ops_seq_v2.h"
#include "dam2d_kernel.inc"
RealField VolumeFraction("vof");
RealField Mass("mass");
RealField Tau("tau");
RealField TotalWeightsWrong("TotalWeightsWrong");
RealField TotalWeightsRight("TotalWeightsRight");
RealField NormalDirections("NormalDirections");
RealField HydrodynamicForce("HydrodynamicForce");
RealField Porosity("Porosity");
RealField ResErrorFresh("ResErrorFresh");
IntField Quality("Quality");
IntField Position("Position");
IntField SegmentNo("SegmentNo");
// RealField VofOld("VofOld");
// RealField MassOld("MassOld");
// RealField fOld("fOld");
// RealField fStageOld("fStageOld");
// IntField NodetypeOld("NodetypeOld");
// IntField GeometryOld("GeometryOld");
// RealField MacroBodyforceOld("MacroBodyforceOld");
// Provide macroscopic initial conditions
//*********************Initialization******************//
ops_reduction TotalDensityError;
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
                ops_arg_dat(g_NodeType().at(compo.id).at(block.ID()), 1,
                            LOCALSTENCIL, "int", OPS_READ),
                ops_arg_idx(),
                ops_arg_dat(g_MacroBodyforce().at(compo.id).at(blockIdx),
                            SpaceDim(), LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_dat(Porosity.at(blockIdx), 1, LOCALSTENCIL,
                                     "Real", OPS_RW));
        }
    }
}

void SetInitialNodeType() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.BulkRange().begin(), block.BulkRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            // Specify general boundary type
            ops_par_loop(KerSetInitialNodeType, "KerSetInitialNodeType",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_NodeType().at(compo.id).at(block.ID()),
                                     1, LOCALSTENCIL, "int", OPS_RW),
                         ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                     LOCALSTENCIL, "Real", OPS_READ),
                         ops_arg_idx(),
                         ops_arg_dat(g_GeometryProperty().at(blockIdx), 1, LOCALSTENCIL,
                        "int", OPS_RW));
        }
    }
}

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

void PreDefinedInitialConditionfree() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const int compoId{compo.id};
            const InitialType initialType{compo.initialType};
            switch (initialType) {
                case Initial_BGKFeq2nd: {
                    ops_par_loop(
                        KerInitialiseBGK2ndfree, "KerInitialiseBGK2ndfree", block.Get(),
                        SpaceDim(), iterRng.data(),
                        ops_arg_dat(g_f()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_fStage()[blockIndex], NUMXI, LOCALSTENCIL,
                                    "double", OPS_WRITE),
                        ops_arg_dat(g_NodeType().at(compoId).at(blockIndex), 1,
                                    LOCALSTENCIL, "int", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                        .at(compo.macroVars.at(Variable_Rho).id)
                                        .at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
                                    1, LOCALSTENCIL, "double", OPS_READ),
                        ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                        ops_arg_dat(
                            g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), LOCALSTENCIL, "double", OPS_READ));
                } break;
                default:
                    ops_printf(
                        "The specified initial type is not implemented!\n");
                    break;
            }
        }
    }
    // TODO this may be better arranged.
    if (!IsTransient()) {
        CopyCurrentMacroVar();
    }
#endif  // OPS_2D
}
//*********************Initialization******************//

//********************LBM evolution*******************//
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
                            KerCalcDensityfree, "KerCalcDensity", block.Get(),
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

void UpdateMacroscopicBodyForcetimeVOF() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            // Specify general boundary type
            ops_par_loop(
                KerUpdateMacroBodyForcetimeVOF, "KerUpdateMacroBodyForceVOF", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(g_MacroBodyforce().at(compo.id).at(blockIdx), SpaceDim(),
                            LOCALSTENCIL, "Real", OPS_RW),
                ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                            LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_idx(),
                ops_arg_dat(VolumeFraction.at(blockIdx), 1, LOCALSTENCIL,
                            "Real", OPS_READ));
        }
    }
}


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
// the present one includes the turbulent model.(note for merge)
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
            const int blocksizeX {block.Size().at(0)};
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
                        ops_arg_dat(Tau[blockIndex], 1, LOCALSTENCIL,
                                    "double", OPS_RW),
                        ops_arg_idx(),
                        ops_arg_gbl(&blocksizeX, 1, "int", OPS_READ));
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
            const Real endblock{(block.Size().at(1)-1)*CS*TimeStep()};
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
                                    ONEPTLATTICESTENCIL, "int", OPS_READ),
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
                        ops_arg_idx(),
                        ops_arg_gbl(&endblock, 1, "double", OPS_READ),
                        ops_arg_dat(Tau.at(blockIndex), 1, LOCALSTENCIL,
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

void Streamfree() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.BulkRange().begin(), block.BulkRange().end());
        const int blockIndex{block.ID()};
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

void ImplementBoundaryFree() {
    for (const auto& boundary : BlockBoundaries()) {
        const Block& block{g_Block().at(boundary.blockIndex)};
        // VertexType BoundaryType{boundary.boundaryType};
        // const int blockSizeY {block.Size().at(1)-1};
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
            default:
                break;
        }
    }
}

void ImplementBoundaryUIBBFree() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.BulkRange().begin(), block.BulkRange().end());
        const int blockIndex{block.ID()};
        for (const auto& compo : g_Components()) {
            ops_par_loop(
                KerUnitInterBounceBack, "KerUnitInterBounceBack", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(g_f().at(blockIndex), NUMXI, LOCALSTENCIL, "double",
                            OPS_RW),
                ops_arg_dat(g_fStage().at(blockIndex), NUMXI,
                            ONEPTLATTICESTENCIL, "double", OPS_READ),
                ops_arg_dat(g_NodeType().at(compo.first).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "int", OPS_READ),
                ops_arg_gbl(compo.second.index, 2, "int", OPS_READ),
                ops_arg_dat(g_MacroVars()
                                .at(compo.second.macroVars.at(Variable_Rho).id)
                                .at(blockIndex),
                            1, LOCALSTENCIL, "double", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.second.uId).at(blockIndex), 1,
                                ONEPTLATTICESTENCIL, "double", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.second.vId).at(blockIndex), 1,
                                ONEPTLATTICESTENCIL, "double", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                                ONEPTLATTICESTENCIL, "int", OPS_READ));
        }
    }
#endif  // OPS_2D
}

void ImplementOpenBoundaryFree(const int iter) {
    for (const auto& boundary : BlockBoundaries()) {
        const Block& block{g_Block().at(boundary.blockIndex)};
        // VertexType BoundaryType{boundary.boundaryType};
        const int blockIndex{block.ID()};
        const int boundaryScheme{(int)boundary.boundaryScheme};
        Component compo{g_Components().at(boundary.componentID)};
        std::vector<int> range(2 * SpaceDim());
        range.assign(
            block.BoundarySurfaceRange().at(boundary.boundarySurface).begin(),
            block.BoundarySurfaceRange().at(boundary.boundarySurface).end());
        switch (boundaryScheme) {
            case (int)BoundaryScheme::BounceBack: {
                // should add another boundaryscheme for openboundary
                ops_par_loop(
                    KerNonEQExtrapolation, "KerNonEQExtrapolation", block.Get(),
                    SpaceDim(), range.data(),
                    ops_arg_dat(g_f()[blockIndex], NUMXI, ONEPTLATTICESTENCIL,
                                "double", OPS_RW),
                    ops_arg_dat(
                        g_NodeType().at(boundary.componentID).at(blockIndex), 1,
                        ONEPTLATTICESTENCIL, "int", OPS_READ),
                    ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                                ONEPTLATTICESTENCIL, "int", OPS_READ),
                    ops_arg_gbl(boundary.givenVars.data(), 2, "double",
                                OPS_READ),
                    ops_arg_dat(
                        g_MacroBodyforce().at(compo.id).at(blockIndex),
                        SpaceDim(), ONEPTLATTICESTENCIL, "double", OPS_READ),
                    ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                    ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                                    LOCALSTENCIL, "double", OPS_READ),
                    ops_arg_gbl(&iter, 1, "int", OPS_READ));
            } break;
            default:
                break;
        }
    }
}
//********************LBM evolution*******************//

//***************Free surface algorithem**************//
void MassExchangeFreeSurface() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        const Real endblock {(block.Size().at(1)-1)*CS*TimeStep()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
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
                                     1, ONEPTLATTICESTENCIL, "double", OPS_READ),
                         ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex),
                                     1, ONEPTLATTICESTENCIL, "double", OPS_READ),
                         ops_arg_idx(),
                         ops_arg_gbl(&jofInletInterface, 1, "int", OPS_READ),
                         ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                         ops_arg_dat(Position.at(blockIndex), 1, LOCALSTENCIL, "int",
                            OPS_READ),
                         ops_arg_gbl(&endblock, 1, "double", OPS_READ),
                         ops_arg_dat(Tau[blockIndex], 1, LOCALSTENCIL,
                                    "double", OPS_READ),
                         ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                            LOCALSTENCIL, "double", OPS_READ),
                         ops_arg_dat(NormalDirections.at(blockIndex), SpaceDim(),
                            LOCALSTENCIL, "Real", OPS_READ));
            ops_par_loop(
                KerUpdatNodeQuality, "KerUpdatNodeQuality", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(Quality.at(blockIndex), 1, LOCALSTENCIL, "int",
                            OPS_RW),
                ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "int", OPS_READ),
                ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                ops_arg_dat(Porosity.at(blockIndex), 1, LOCALSTENCIL, "Real",
                            OPS_READ),
                ops_arg_dat(NormalDirections.at(blockIndex), SpaceDim(),
                            LOCALSTENCIL, "Real", OPS_READ),
                ops_arg_dat(VolumeFraction.at(blockIndex), 1, LOCALSTENCIL,
                            "Real", OPS_READ),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                            LOCALSTENCIL, "double", OPS_READ),
                ops_arg_dat(SegmentNo.at(blockIndex), 1, LOCALSTENCIL, "int",
                            OPS_READ));

            ops_par_loop(KerMassExchange, "KerMassExchange", block.Get(),
                         SpaceDim(), iterRng.data(),
                         ops_arg_dat(Mass.at(blockIndex), 1, LOCALSTENCIL,
                                     "Real", OPS_RW),
                         ops_arg_dat(g_fStage()[blockIndex], NUMXI,
                                     ONEPTLATTICESTENCIL, "double", OPS_READ),
                         ops_arg_dat(g_f()[blockIndex], NUMXI,
                                     ONEPTLATTICESTENCIL, "double", OPS_READ),
                         ops_arg_dat(VolumeFraction.at(blockIndex), 1,
                                     ONEPTLATTICESTENCIL, "Real", OPS_READ),
                         ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex),
                                     1, ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_dat(Quality.at(blockIndex), 1,
                                     ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                                     ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_idx(),
                         ops_arg_gbl(&jofInletInterface, 1, "int", OPS_READ),
                         ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                         ops_arg_dat(Porosity.at(blockIndex), 1, 
                                     ONEPTLATTICESTENCIL, "Real", OPS_READ),
                         ops_arg_dat(NormalDirections.at(blockIndex), 2,
                                     ONEPTLATTICESTENCIL, "Real", OPS_READ),
                        ops_arg_dat(g_MacroVars()
                                .at(compo.macroVars.at(Variable_Rho).id)
                                .at(blockIndex),
                            1, ONEPTLATTICESTENCIL, "double", OPS_READ));


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
        const Real endblock {(block.Size().at(1)-1)*CS*TimeStep()};//block mesh numbers
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const Real tau{compo.tauRef};
            ops_par_loop(KerPreUpdateFlag, "KerPreUpdateFlag", block.Get(),
                         SpaceDim(), iterRng.data(),
                         ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex),
                                     1, LOCALSTENCIL, "int", OPS_RW),
                         ops_arg_dat(VolumeFraction.at(blockIndex), 1,
                                     ONEPTLATTICESTENCIL, "Real", OPS_READ),
                         ops_arg_dat(Quality.at(blockIndex), 1,
                                     ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_idx(),
                         ops_arg_gbl(&jofInletInterface, 1, "int", OPS_READ),
                         ops_arg_gbl(compo.index, 2, "int", OPS_READ));

            ops_par_loop(
                KerBuildNewInterface, "KerBuildNewInterface", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "int", OPS_RW),
                ops_arg_idx(),
                ops_arg_gbl(&jofInletInterface, 1, "int", OPS_READ),
                ops_arg_gbl(compo.index, 2, "int", OPS_READ));
            ops_par_loop(
                KerInitNewInterface, "KerInitNewInterface", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(g_f()[blockIndex], NUMXI, ONEPTLATTICESTENCIL,
                            "double", OPS_RW),
                ops_arg_dat(g_fStage()[blockIndex], NUMXI, ONEPTLATTICESTENCIL,
                            "double", OPS_RW),            
                ops_arg_dat(g_MacroVars()
                                .at(compo.macroVars.at(Variable_Rho).id)
                                .at(blockIndex),
                            1, ONEPTLATTICESTENCIL, "double", OPS_RW),
                ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "double", OPS_RW),
                ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "double", OPS_RW),
                ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "int", OPS_READ),
                ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                ops_arg_gbl(&endblock, 1, "Real", OPS_READ),
                ops_arg_gbl(&tau, 1, "double", OPS_READ),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                            LOCALSTENCIL, "double", OPS_READ),
                ops_arg_dat(g_MacroBodyforce().at(compo.id).at(blockIndex),
                            SpaceDim(), ONEPTLATTICESTENCIL, "double",
                            OPS_READ));
#ifdef FreeMpi_2D
            ops_par_loop(
                KerCalcTotalWeight, "KerCalcTotalWeight", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(TotalWeightsWrong.at(blockIndex), 1, LOCALSTENCIL,
                            "Real", OPS_RW),
                ops_arg_dat(TotalWeightsRight.at(blockIndex), 1, LOCALSTENCIL,
                            "Real", OPS_RW),
                ops_arg_dat(NormalDirections.at(blockIndex), 2, LOCALSTENCIL,
                            "Real", OPS_RW),
                ops_arg_dat(VolumeFraction.at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "Real", OPS_READ),
                ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "int", OPS_READ),
                ops_arg_dat(g_GeometryProperty()[blockIndex], 1,
                            ONEPTLATTICESTENCIL, "int", OPS_READ),
                ops_arg_dat(Quality.at(blockIndex), 1, LOCALSTENCIL, "int",
                            OPS_READ),
                ops_arg_idx(),
                ops_arg_gbl(&jofInletInterface, 1, "int", OPS_READ),
                ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                ops_arg_dat(g_CoordinateXYZ()[blockIndex], SpaceDim(),
                            LOCALSTENCIL, "double", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "double", OPS_RW),
                ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "double", OPS_RW));

            ops_par_loop(
                KerAllocateMassMpi, "KerAllocateMassMpi", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(Mass.at(blockIndex), 1, ONEPTLATTICESTENCIL, "Real",
                            OPS_RW),
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
                ops_arg_dat(NormalDirections.at(blockIndex), 2,
                            ONEPTLATTICESTENCIL, "Real", OPS_READ),
                ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                ops_arg_dat(VolumeFraction.at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "Real", OPS_READ),
                ops_arg_dat(Porosity.at(blockIndex), 1, ONEPTLATTICESTENCIL,
                            "Real", OPS_READ),
                ops_arg_dat(g_f()[blockIndex], NUMXI, ONEPTLATTICESTENCIL, 
                            "double", OPS_READ),
                ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "double", OPS_RW),
                ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "double", OPS_RW));
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
                ops_arg_dat(Mass.at(blockIdx), 1, LOCALSTENCIL, "Real", OPS_RW),
                ops_arg_dat(VolumeFraction.at(blockIdx), 1, LOCALSTENCIL,
                            "Real", OPS_RW),
                ops_arg_dat(g_NodeType().at(compo.id).at(blockIdx), 1,
                            LOCALSTENCIL, "int", OPS_RW),
                ops_arg_dat(g_MacroVars()
                                .at(compo.macroVars.at(Variable_Rho).id)
                                .at(blockIdx),
                            1, LOCALSTENCIL, "double", OPS_READ),
                ops_arg_dat(Porosity.at(blockIdx), 1, LOCALSTENCIL, "Real",
                            OPS_READ),
                ops_arg_dat(g_f()[blockIdx], NUMXI, LOCALSTENCIL, "double",
                            OPS_READ),
                ops_arg_gbl(compo.index, 2, "int", OPS_READ));
        }
    }
}

void CalcNormalDirectionFree() {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            ops_par_loop(KerCalcNormalDirection, "KerCalcNormalDirection",
                         block.Get(), SpaceDim(), iterRng.data(),
                         ops_arg_dat(NormalDirections.at(blockIdx), 2,
                                     LOCALSTENCIL, "Real", OPS_RW),
                         ops_arg_dat(VolumeFraction.at(blockIdx), 1,
                                     ONEPTLATTICESTENCIL, "Real", OPS_READ),
                         ops_arg_dat(g_NodeType().at(compo.id).at(blockIdx),
                                     1, ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_dat(g_GeometryProperty()[blockIdx], 1,
                                     ONEPTLATTICESTENCIL, "int", OPS_READ),
                         ops_arg_gbl(compo.index, 2, "int", OPS_READ));
        }
    }
}

void ArtificialreduceSpaBub() {
#ifdef OPS_2D
    for (const auto& idBlock : g_Block()) {
        const Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIndex{block.ID()};
        for (const auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            ops_par_loop(
                KerartificialSpatterBubble, "KerartificialSpatterBubble", block.Get(),
                SpaceDim(), iterRng.data(),
                ops_arg_dat(g_NodeType().at(compo.id).at(blockIndex), 1,
                            ONEPTLATTICESTENCIL, "int", OPS_RW),
                ops_arg_gbl(compo.index, 2, "int", OPS_READ));
        }
    }
#endif
}

void SubIterationIBB(const Real convergenceCriteria,
                     const int IterSteps) {
    for (auto idBlock : g_Block()) {
        Block& block{idBlock.second};
        std::vector<int> iterRng;
        iterRng.assign(block.WholeRange().begin(), block.WholeRange().end());
        const int blockIdx{block.ID()};
        // const Real endblock {(block.Size().at(0)-1)*CS*TimeStep()};//block mesh numbers
        for (auto& idCompo : g_Components()) {
            const Component& compo{idCompo.second};
            const Real tau{compo.tauRef};
            const Real* pdt{pTimeStep()};
            for (int SubIter = 0; SubIter <= IterSteps; SubIter++) {
                Real Error{0};
                ops_par_loop(
                    KerSubCollideIBB, "KerSubCollideIBB", block.Get(),
                    SpaceDim(), iterRng.data(),
                    ops_arg_dat(g_fStage()[blockIdx], NUMXI, LOCALSTENCIL,
                                "double", OPS_RW),
                    ops_arg_dat(g_f()[blockIdx], NUMXI, LOCALSTENCIL, "double",
                                OPS_READ),
                    ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                LOCALSTENCIL, "double", OPS_READ),
                    ops_arg_dat(g_NodeType().at(compo.id).at(blockIdx), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                    ops_arg_dat(g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIdx),
                                1, LOCALSTENCIL, "double", OPS_READ),
                    ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIdx), 1,
                                LOCALSTENCIL, "double", OPS_READ),
                    ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIdx), 1,
                                LOCALSTENCIL, "double", OPS_READ),
                    ops_arg_gbl(&tau, 1, "double", OPS_READ),
                    ops_arg_gbl(pdt, 1, "double", OPS_READ),
                    ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                    ops_arg_dat(g_MacroBodyforce().at(compo.id).at(blockIdx),
                                SpaceDim(), LOCALSTENCIL, "double", OPS_READ));
                ops_par_loop(
                    KerSubStreamIBB, "KerSubStreamIBB", block.Get(), SpaceDim(),
                    iterRng.data(),
                    ops_arg_dat(g_f()[blockIdx], NUMXI, LOCALSTENCIL, "double",
                                OPS_RW),
                    ops_arg_dat(g_fStage()[blockIdx], NUMXI,
                                ONEPTLATTICESTENCIL, "double", OPS_READ),
                    ops_arg_dat(g_CoordinateXYZ()[blockIdx], SpaceDim(),
                                LOCALSTENCIL, "double", OPS_READ),
                    ops_arg_dat(g_NodeType().at(compo.id).at(blockIdx), 1,
                                ONEPTLATTICESTENCIL, "int", OPS_READ),
                    ops_arg_dat(g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIdx),
                                1, ONEPTLATTICESTENCIL, "double", OPS_READ),
                    ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIdx), 1,
                                ONEPTLATTICESTENCIL, "double", OPS_READ),
                    ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIdx), 1,
                                ONEPTLATTICESTENCIL, "double", OPS_READ),
                    ops_arg_gbl(compo.index, 2, "int", OPS_READ), ops_arg_idx(),
                    ops_arg_dat(VolumeFraction.at(blockIdx), 1,
                                ONEPTLATTICESTENCIL, "Real", OPS_READ),
                    ops_arg_dat(SegmentNo.at(blockIdx), 1, LOCALSTENCIL,
                                "int", OPS_READ),
                    ops_arg_dat(Position.at(blockIdx), 1, LOCALSTENCIL, "int",
                            OPS_READ),
                    ops_arg_dat(NormalDirections.at(blockIdx), SpaceDim(),
                            LOCALSTENCIL, "Real", OPS_READ),
                    ops_arg_dat(Mass.at(blockIdx), 1, LOCALSTENCIL,
                                "double", OPS_READ));
                ops_par_loop(
                    KerSubCalcMacroVarsIBB, "KerSubCalcMacroVarsIBB",
                    block.Get(), SpaceDim(), iterRng.data(),
                    ops_arg_dat(g_MacroVars()
                                    .at(compo.macroVars.at(Variable_Rho).id)
                                    .at(blockIdx),
                                1, LOCALSTENCIL, "double", OPS_RW),
                    ops_arg_dat(g_MacroVars().at(compo.uId).at(blockIdx), 1,
                                LOCALSTENCIL, "double", OPS_RW),
                    ops_arg_dat(g_MacroVars().at(compo.vId).at(blockIdx), 1,
                                LOCALSTENCIL, "double", OPS_RW),
                    ops_arg_dat(Mass.at(blockIdx), 1, LOCALSTENCIL,
                                "double", OPS_RW),
                    ops_arg_dat(VolumeFraction.at(blockIdx), 1,
                                ONEPTLATTICESTENCIL, "Real", OPS_RW),
                    ops_arg_dat(g_NodeType().at(compo.id).at(blockIdx), 1,
                                ONEPTLATTICESTENCIL, "int", OPS_READ),
                    ops_arg_dat(Porosity.at(blockIdx), 1, LOCALSTENCIL,
                                "double", OPS_READ),
                    ops_arg_dat(g_f()[blockIdx], NUMXI, LOCALSTENCIL, "double",
                                OPS_READ),
                    ops_arg_gbl(compo.index, 2, "int", OPS_READ),
                    ops_arg_dat(g_MacroBodyforce().at(compo.id).at(blockIdx),
                                SpaceDim(), LOCALSTENCIL, "double", OPS_READ),
                    ops_arg_dat(ResErrorFresh.at(blockIdx), 1,
                                LOCALSTENCIL, "double", OPS_RW));
                ops_par_loop(
                    KerTotalDensityError, "KerTotalDensityError", block.Get(),
                    SpaceDim(), iterRng.data(),
                    ops_arg_dat(g_NodeType().at(compo.id).at(blockIdx), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                    ops_arg_dat(ResErrorFresh.at(blockIdx), 1, LOCALSTENCIL,
                                "Real", OPS_READ),
                    ops_arg_reduce(TotalDensityError, 1, "double", OPS_INC));
                    ops_reduction_result(TotalDensityError, &Error);
                
                if (Error < convergenceCriteria||SubIter >= IterSteps) {
                    break;
                }
            }
            ops_par_loop(
                    KerSetDFsofNewInterface, "KerSetDFsofNewInterface", block.Get(),
                    SpaceDim(), iterRng.data(),
                    ops_arg_dat(g_f()[blockIdx], NUMXI, LOCALSTENCIL,
                                "double", OPS_RW),
                    ops_arg_dat(g_fStage()[blockIdx], NUMXI, LOCALSTENCIL, "double",
                                OPS_READ),
                    ops_arg_dat(g_NodeType().at(compo.id).at(blockIdx), 1,
                                LOCALSTENCIL, "int", OPS_READ),
                    ops_arg_gbl(compo.index, 2, "int", OPS_READ));
        }
    }
}

    // Address MRT transmission Matrix//
    typedef std::vector<Real> vec;
    typedef std::vector<vec> mat;
    enum Moments {
    rho = 0,
    e = 1,
    E = 2,
    jx = 3,
    qx = 4,
    jy = 5,
    qy = 6,
    pxx = 7,
    pxy = 8
    };

    void TransformationM2D(const int nc) {
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                Real cx = XI[j * LATTDIM];
                Real cy = XI[j * LATTDIM + 1];
                switch ((Moments)i) {
                    case rho: {
                        TransformationM[i][j] = 1;
                    } break;
                    case e: {
                        TransformationM[i][j] =
                            3*(cx * cx + cy * cy) - 4;
                    } break;
                    case E: {
                        TransformationM[i][j] =
                            (9 * pow((cx * cx + cy * cy), 2) -
                             21 * (cx * cx + cy * cy) + 8)/2.0;
                    } break;
                    case jx: {
                        TransformationM[i][j] = cx;
                    } break;
                    case qx: {
                        TransformationM[i][j] =
                            cx * (3 * (cx * cx + cy * cy) - 5);
                    } break;
                    case jy: {
                        TransformationM[i][j] = cy;
                    } break;
                    case qy: {
                        TransformationM[i][j] =
                            cy * (3 * (cx * cx + cy * cy) - 5);
                    } break;
                    case pxx: {
                        TransformationM[i][j] =
                            cx * cx - cy * cy;
                    } break;
                    case pxy: {
                        TransformationM[i][j] = cx * cy;
                    } break;
                    default: {
                    } break;
                }
            }
        }
    }

    // void TransformationM2D(const int nc) {
    //     for (int i = 0; i < nc; i++) {
    //         for (int j = 0; j < nc; j++) {
    //             Real cx = (int)XI[j * LATTDIM];
    //             Real cy = (int)XI[j * LATTDIM + 1];
    //             switch ((Moments)i) {
    //                 case rho: {
    //                     TransformationM[i][j] = 1;
    //                 } break;
    //                 case e: {
    //                     TransformationM[i][j] =
    //                         3*(cx * cx + cy * cy) - 2;
    //                 } break;
    //                 case E: {
    //                     TransformationM[i][j] =
    //                         9 * pow((cx * cx + cy * cy), 2) -
    //                          15 * (cx * cx + cy * cy) + 2;
    //                 } break;
    //                 case jx: {
    //                     TransformationM[i][j] = cx;
    //                 } break;
    //                 case qx: {
    //                     TransformationM[i][j] =
    //                         cx * (3 * (cx * cx + cy * cy) - 4);
    //                 } break;
    //                 case jy: {
    //                     TransformationM[i][j] = cy;
    //                 } break;
    //                 case qy: {
    //                     TransformationM[i][j] =
    //                         cy * (3 * (cx * cx + cy * cy) - 4);
    //                 } break;
    //                 case pxx: {
    //                     TransformationM[i][j] =
    //                         cx * cx - cy * cy;
    //                 } break;
    //                 case pxy: {
    //                     TransformationM[i][j] = cx * cy;
    //                 } break;
    //                 default: {
    //                 } break;
    //             }
    //         }
    //     }
    // }

    Real det(mat cofactorm, const int n){
        Real value{1};
        for (int i = 0; i < n-1; i++) {
            vec temp(n);
            int k {0};
            for (k = i; k < n; k++){
                if (abs(cofactorm[k][i])>1e-10){
                    break;
                }
            }
            if (k<n){
                if (k!=i){
                    for (int j = i; j < n; j++) {
                        temp[j] = cofactorm[i][j];
                        cofactorm[i][j] = cofactorm[k][j];
                        cofactorm[k][j] = temp[j];
                    }
                    value = -value;
                }
                if (abs(cofactorm[i][i])<1e-10){
                    return 0;
                }
                for (int j = i+1; j < n; j++){
                    Real a = -cofactorm[j][i]/cofactorm[i][i];
                    for (int k = i; k < n; k++){
                        cofactorm[j][k] += a*cofactorm[i][k];
                    }
                }
            }
            value = value *cofactorm[i][i];
        }
        return value*cofactorm[n-1][n-1];
    }//O(n^3)

    void CalcInverseM(const int nc) {
        mat Cofactor(nc, vec(nc));
        mat CofactorM(nc-1, vec(nc-1));
        mat matTransformationM(nc, vec(nc));
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                matTransformationM[i][j] = TransformationM[i][j];
            }
        }
        Real detofTransformationM{det(matTransformationM, nc)};
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                for (int k = 0; k < nc - 1; k++) {
                    for (int t = 0; t < nc - 1; t++) {
                        CofactorM[k][t] = 
                            TransformationM[k >= i ? k + 1 : k]
                                            [t >= j ? t + 1 : t];
                    }
                }
                Cofactor[j][i] = pow(-1,i+j)*det(CofactorM,nc - 1);
                InverseM[j][i] = Cofactor[j][i]/detofTransformationM;
                if (abs(InverseM[j][i])<1e-10){
                    InverseM[j][i] = 0;
                }
            }
        }
        for (int i = 0; i < nc; i++) {
            Real b{0};
            for (int j = 0; j < nc; j++) {
                b += InverseM[i][j]*TransformationM[j][i];
            }
            if (abs(b - 1) > 1e-10) {
                ops_printf("Wrong inverseM b=%e,i=%d\n", b,i);
                assert(abs(b - 1) < 1e-10);
            }
        }
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                Real b{0};
                if (i != j) {
                    for (int k = 0; k < nc; k++) {
                        b += InverseM[i][k] * TransformationM[k][j];
                    }
                    if (abs(b) > 1e-10) {
                        ops_printf("Wrong inverseM b=%e,i=%d\n", b, i);
                        assert(abs(b) < 1e-10);
                    }
                }
            }
        }
        // for (int i = 0; i < nc; i++) {
        //     for (int j = 0; j < nc; j++) {
        //         // if (InverseM1[i][j]!=InverseM[i][j]){
        //              printf("%2.16f",InverseM1[i][j]-InverseM[i][j]);
        //         // }
        //     }
        //     printf("Wrong inverseM \n");
        // }
    }

    void MatrixforMRT2D(const int nc9) {
    for (int i = 3; i < nc9; i++) {
        for (int j = 0; j < nc9; j++) {
            if (i < 7) {
                TransformationM[i][j] *= CS;  // D2Q9
                InverseM[j][i] /= CS;
            } else {
                TransformationM[i][j] *= CS*CS;
                InverseM[j][i] /= CS*CS;
            }
        }
    }
}
    // address MRT transmission Matrix//

    void StreamCollisionDam(const Real time, const int iter,
                            const int RK) {
#if DebugLevel >= 1
    ops_printf("Calculating the macroscopic variables...\n");
#endif
    UpdateMacroscopicBodyForce(time);
    UpdateMacroscopicBodyForcetimeVOF();
#ifdef OPS_2D
    UpdateMacroVarsfree();
#endif
    // CopyBlockEnvelopDistribution(g_fStage(), g_f());
#if DebugLevel >= 1
    ops_printf("Calculating the mesoscopic body force term...\n");
#endif



#ifdef OPS_2D
    PreDefinedBodyForcefree();
#endif
#if DebugLevel >= 1
    ops_printf("Calculating the collision term...\n");
#endif
#ifdef OPS_2D
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
    // ImplementBoundaryFree();
#endif

#ifdef OPS_2D
    Streamfree();
    ImplementBoundaryUIBBFree();//less stable than halfway BB

    CalcNormalDirectionFree();
#endif    
#ifdef OPS_2D
    MassExchangeFreeSurface();    
#endif
    CalcVolumeFraction();

#ifdef OPS_2D
    UpdateFlagAllocateMass();
    // SubIterationIBB(1e-12,0);
    CalcVolumeFraction();
    // ArtificialreduceSpaBub();
    // ImplementOpenBoundaryFree(iter);
#endif


}


template <typename T>
void Iteratefree(void (*cycle)(T,const int iter, const int rk), const SizeType steps,
             const SizeType checkPointPeriod, const SizeType start = 0) {
    ops_printf("Starting the iteration...\n"); 
    for (SizeType iter = start; iter < start + steps; iter++) {
        const Real time{iter * TimeStep()};
        int RK {0};
        cycle(time, iter, RK);
        // if (iter>1000){
        if (((iter + 1) % checkPointPeriod) == 0) {
            ops_printf("%d iterations!\n", iter + 1);

#ifdef OPS_3D
            UpdateMacroVars3D();
#endif
#ifdef OPS_2D
            // UpdateMacroVarsfree();
#endif
            WriteFlowfieldToHdf5((iter + 1));
            // WriteDistributionsToHdf5((iter + 1));
            WriteNodePropertyToHdf5((iter + 1));
            VolumeFraction.WriteToHDF5("Dam", iter + 1);
            NormalDirections.WriteToHDF5("Dam", iter + 1);
            // Mass.WriteToHDF5("Dam", iter + 1);
            // Quality.WriteToHDF5("Dam", iter + 1);
        }
        // }
    }
    ops_printf("Simulation finished! Exiting...\n");
    DestroyModel();
}

template <typename T>
void Iteratefree(void (*cycle)(T,const int iter,const int rk), const Real convergenceCriteria,
             const SizeType checkPointPeriod, const SizeType start = 0) {
    SizeType iter{start};
    Real residualError{1};
    do {
        const Real time{iter * TimeStep()};
        int RK{0};
        cycle(time,iter,RK);
        iter = iter + 1;
        if ((iter % checkPointPeriod) == 0) {
#ifdef OPS_3D
            UpdateMacroVars3D();
#endif
#ifdef OPS_2D
            // UpdateMacroVarsfree();
#endif
            CalcResidualError();
            residualError = GetMaximumResidual(checkPointPeriod);
            DispResidualError(iter+1, checkPointPeriod);
            WriteFlowfieldToHdf5(iter+1);
            // WriteDistributionsToHdf5(iter+1);
            WriteNodePropertyToHdf5(iter+1);
            VolumeFraction.WriteToHDF5("Dam", iter+1);
            Mass.WriteToHDF5("Dam", iter+1);
            Porosity.WriteToHDF5("Dam", iter+1);
        }
    } while (residualError >= convergenceCriteria);

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
    Position.CreateFieldFromScratch(g_Block());
    SegmentNo.CreateFieldFromScratch(g_Block());
    Tau.CreateFieldFromScratch(g_Block());
    TotalWeightsWrong.CreateFieldFromScratch(g_Block());
    TotalWeightsRight.CreateFieldFromScratch(g_Block());
    // TotalWeights.SetDataDim(2);
    NormalDirections.SetDataDim(SpaceDim());
    NormalDirections.CreateFieldFromScratch(g_Block());
    HydrodynamicForce.SetDataDim(3);
    HydrodynamicForce.CreateFieldFromScratch(g_Block());
    Porosity.CreateFieldFromScratch(g_Block());
    ResErrorFresh.CreateFieldFromScratch(g_Block());
    TotalDensityError = ops_decl_reduction_handle(
                    sizeof(double), "double", "DError");
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
    SetTimeStep(config.meshSize / SoundSpeed());
    ops_diagnostic_output();
    if (config.currentTimeStep == 0) {
        UpdateMacroscopicBodyForce();
        SetSlipBoundaryNodeType();
        SetInitialNodeType();
        SetInitialMacrosVars();
        PreDefinedInitialConditionfree();
        
        TransformationM2D(9);
        CalcInverseM(9);
        MatrixforMRT2D(9);
    };
    if (config.transient) {
        Iteratefree(StreamCollisionDam, config.timeStepsToRun,
                config.checkPeriod, config.currentTimeStep);
    } else {
        Iteratefree(StreamCollisionDam, config.convergenceCriteria,
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
