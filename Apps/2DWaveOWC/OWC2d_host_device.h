
#ifndef FREESURFACE_HOST_DEVICE_H
#define FREESURFACE_HOST_DEVICE_H
#ifndef OPS_FUN_PREFIX
#define OPS_FUN_PREFIX
#endif

//*****wave parameters*****//
// std::string CaseName;
Real WaterDepth;
Real WaveAmplitude;
Real WavenumberKH;
Real WaveLength;
Real WavePeriod;
Real V0;
Real SoundofWater;
Real Ma;
Real L0;//water depth
// wave absorption
Real LinearConst;  // unit 1/T
Real ns;  // empirical coefficient
Real DampWidth;
Real AbsorbLength;
Real AbsorbConst;
Real Gravity;//1.0791e-4
Real Absorption[3]{};
Real BlockCoor[2]{0};
//*****wave parameters*****//

//*****aerodynamic parameters*****//
const Real AirCoe{1.4};
const Real PhysicalRhoOfFluid{1000};
Real ReferencePOfFluid;
Real RefRhoOfAir;
Real RefPOfAir;
Real ChamberCoef;
//*****aerodynamic parameters*****//

static inline OPS_FUN_PREFIX Real CalcBGKFeqFA(
    const int l, const Real rho, const Real u, const Real v, const Real T,
    const int polyOrder, const Real tau, const Real iCoor) {
    Real cu{(CS * XI[l * LATTDIM] * u + CS * XI[l * LATTDIM + 1] * v)};
    Real c2{(CS * XI[l * LATTDIM] * CS * XI[l * LATTDIM] +
             CS * XI[l * LATTDIM + 1] * CS * XI[l * LATTDIM + 1])};
    Real cu2{cu * cu};
    Real u2{u * u + v * v};
    Real LinearAbsorb{0};
    if (iCoor <= BlockCoor[0]+DampWidth|| 
        iCoor >= BlockCoor[1]-DampWidth ) {
        //|| 
        // iCoor >= blockcoor[1]-DampWidth
        // Real NonLinearConst{10 * L0};       // unit 1/L
        Real StartPt{0};
        if (iCoor <= BlockCoor[0]+DampWidth) {
            StartPt = BlockCoor[0]+DampWidth;
        } else if (iCoor >= BlockCoor[1]-DampWidth) {
            StartPt = BlockCoor[1]-DampWidth;
        }
        Real DistanceScale = (iCoor - StartPt) / DampWidth;
        LinearAbsorb =TimeStep()*
            (-LinearConst * (exp(pow(DistanceScale, ns)) - 1) / (exp(1) - 1));
            // ops_printf("Updating the halos...%e, %e\n",  DampWidth,iCoor);   
    }
    Real dtOvertauPlusdt = TimeStep() / (tau + 0.5 * TimeStep());
    Real res = 1.0 + cu * (1 + LinearAbsorb * 1/dtOvertauPlusdt) +
               0.5 * (cu2 - u2 + (T - 1.0) * (c2 - LATTDIM)) *
               (1 + LinearAbsorb * 2/dtOvertauPlusdt);
            
    if ((polyOrder) >= 3) {
        res = res +
              cu * (cu2 - 3.0 * u2 + 3.0 * (T - 1.0) * (c2 - LATTDIM - 2.0)) /
                  6.0;
    }
    if ((polyOrder) >= 4) {
        res =
            res + (cu2 * cu2 - 6.0 * cu2 * u2 + 3.0 * u2 * u2) / 24.0 +
            (T - 1.0) * ((c2 - (LATTDIM + 2)) * (cu2 - u2) - 2.0 * cu2) / 4.0 +
            (T - 1.0) * (T - 1.0) *
                (c2 * c2 - 2.0 * (LATTDIM + 2) * c2 + LATTDIM * (LATTDIM + 2)) /
                8.0;
    }
    return WEIGHTS[l] * rho * res;
}

static inline OPS_FUN_PREFIX void CalcNormalDirection(Real* NormalDirection,
                const int (*vgs)[3],const Real (*vofs)[3], const int* lattIdx) {
    Real weight;
    for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
        int cx = (int)XI[xiIndex * LATTDIM];
        int cy = (int)XI[xiIndex * LATTDIM + 1];
        Real vof{vofs[cx+1][cy+1]};
        VertexGeometryType vg = (VertexGeometryType) vgs[cx+1][cy+1];
        if (0 == cx||0 == cy){
            weight = 2.0;
        }else{
            weight = 1.0;
        }
        switch(vg){
            case VG_IP: {
                    vof = vofs[cx+1-cx][cy+1];
            } break;
            case VG_IM: {
                    vof = vofs[cx+1-cx][cy+1];
            } break;
            case VG_JP: {
                    vof = vofs[cx+1][cy+1-cy];
            } break;
            case VG_JM: {
                    vof = vofs[cx+1][cy+1-cy];
            } break;
            default:{
            } break;
        }
        if (VG_IPJP_I == vg || VG_IPJM_I == vg || 
            VG_IMJP_I == vg || VG_IMJM_I == vg || 
            VG_IPJP_O == vg || VG_IPJM_O == vg || 
            VG_IMJP_O == vg || VG_IMJM_O == vg) {
                    NormalDirection[0] =
                        0.5 * (vofs[1][1] - vofs[1 - cx][1]) * (-cx);
                    NormalDirection[1] =
                        0.5 * (vofs[1][1] - vofs[1][1 - cy]) * (-cy);
                    break;
        } 
        else if (VG_ImmersedSolid == vg){
                    NormalDirection[0] =
                        0.5 * (vofs[0][1] - vofs[2][1]);
                    NormalDirection[1] =
                        0.5 * (vofs[1][0] - vofs[1][2]);
                    break;
        }
         else {
            NormalDirection[0] += -weight * cx * vof / 8.0;
            NormalDirection[1] += -weight * cy * vof / 8.0;
        }
    }
}

static inline OPS_FUN_PREFIX Real CalcBodyForcefree(const int xiIndex, const Real rho,
                   const Real* acceleration,const Real& u, const Real& v) {
    Real cf{0};
    for (int i = 0; i < LATTDIM; i++) {
        cf += CS * XI[xiIndex * LATTDIM + i] * acceleration[i];
    }
    Real cu{
        (CS * XI[xiIndex * LATTDIM] * u + CS * XI[xiIndex * LATTDIM + 1] * v)};
    Real uf{u * acceleration[0] + v * acceleration[1]};
    return WEIGHTS[xiIndex] * rho * (cf + cu * cf - uf);
}

static inline OPS_FUN_PREFIX Real Smagorinsky_LES(
                const Real rho, const Real u, const Real v, 
                const Real tau,const Real* f, const int* lattIdx) {
        Real tautotal;
        const Real MeshSize = CS*TimeStep();
        const Real C{0.1*0.1};
        const Real U[LATTDIM]{u,v};
        // const Real T{1};
        // const int polyOrder{2};
        Real NonEquilibriumMoment[LATTDIM*LATTDIM]{0};
        Real NonEquilibriumMoment2{0.0};
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            Real cx = CS * XI[xiIndex * LATTDIM];
            Real cy = CS * XI[xiIndex * LATTDIM + 1];
            Real Lattice_Direction[LATTDIM]{cx,cy};
            // const Real feq{CalcBGKFeq(xiIndex, rho, u, v,w, T, polyOrder)};
            for (int i =0; i< LATTDIM; i++){
                for (int j =0; j< LATTDIM; j++){
                    NonEquilibriumMoment[i*LATTDIM+j] += 
                        Lattice_Direction[i]*Lattice_Direction[j]*
                            f[xiIndex] ;
                }
            }
        }
        for (int l=0;l<LATTDIM*LATTDIM;l++){
            int j = l % LATTDIM;
            int i = l / LATTDIM;
            NonEquilibriumMoment[i*LATTDIM+j] -= (
                        rho*((i==j? 1.0 : 0.0)+U[i]*U[j]));
            NonEquilibriumMoment2 +=
                2 * NonEquilibriumMoment[l] * NonEquilibriumMoment[l];
        }
        tautotal = 0.5 * (tau - 0.5 * TimeStep() +
                     sqrt((tau + 0.5 * TimeStep()) * (tau + 0.5 * TimeStep()) +
                          2 * C * MeshSize * MeshSize *
                              sqrt(NonEquilibriumMoment2) / rho));
    return tautotal;
}

static inline OPS_FUN_PREFIX Real Smagorinsky_LES(
                const Real rho, const Real u, const Real v, 
                const Real tau,const Real* f, const int* lattIdx,
                const Real *Mf,const Real InitialTau, 
                const Real LinearAbsorb, const Real IterSteps=5) {
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
        Real tautotal{InitialTau};
        const Real MeshSize = CS*TimeStep();
        const Real C{0.1*0.1};
        const Real U[LATTDIM]{u,v};
        // const Real T{1};
        // const int polyOrder{2};
        Real EquilibriumMoment[lattIdx[1]+1]{rho, 
                                        rho * (-2 + (u * u + v * v) ),
                                        rho * (1 - (u * u + v * v) ),
                                        rho * u, -rho * u,
                                        rho * v, -rho * v,
                                        rho * (u * u - v * v), 
                                        rho * u * v};
    Real TAU{InitialTau};
    for (int SubIter = 0; SubIter <= IterSteps; SubIter++) {
        Real NonEquilibriumMoment[LATTDIM*LATTDIM]{0};
        Real NonEquilibriumMoment2{0.0};
        Real MfStage[lattIdx[1] + 1]{0};
        Real fStage[lattIdx[1] + 1]{0};
        Real dtOvertauPlusdt = (TimeStep()) / (TAU + 0.5 * (TimeStep()));
        // relaxation moment
        Real RelaxationRate[lattIdx[1]+1]{1, 
                                        0.3, 
                                        1, 
                                        1, 1, 
                                        1, 1, 
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
        for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
                MfStage[i] =
                    -RelaxationRate[i] * 
                    (Mf[i] -  AbsorbWave[i]*EquilibriumMoment[i]);
        }
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            Real InversefStage{0};
            for (int i = lattIdx[0]; i <= lattIdx[1]; i++) {
                InversefStage += InverseM[xiIndex][i] * MfStage[i];
            }
            fStage[xiIndex] = InversefStage;
        }
        for (int xiIndex = lattIdx[0]; xiIndex <= lattIdx[1]; xiIndex++) {
            Real cx = CS * XI[xiIndex * LATTDIM];
            Real cy = CS * XI[xiIndex * LATTDIM + 1];
            Real Lattice_Direction[LATTDIM]{cx,cy};
            // const Real feq{CalcBGKFeq(xiIndex, rho, u, v,w, T, polyOrder)};
            for (int i =0; i< LATTDIM; i++){
                for (int j =0; j< LATTDIM; j++){
                    NonEquilibriumMoment[i*LATTDIM+j] += 
                        Lattice_Direction[i]*Lattice_Direction[j]*
                            fStage[xiIndex];
                }
            }
        }
        for (int l=0;l<LATTDIM*LATTDIM;l++){
            NonEquilibriumMoment2 +=
                2 * NonEquilibriumMoment[l] * NonEquilibriumMoment[l];
        }
        tautotal = TAU;
        TAU = tau + C * MeshSize * MeshSize * 
                        sqrt(NonEquilibriumMoment2) / rho / (2*TimeStep());
        if (abs(TAU - tautotal) < tau||SubIter >= IterSteps) {
                    break;
        }
    }
    return tautotal;
}

static inline OPS_FUN_PREFIX int InPolygon(int& seriesnum, const Real* P,
                                            const Real* RotationCenter,
                                            const Real& Radii) {
    Real distance{0};
    for (int i = 0; i < SpaceDim(); i++) {
        distance += pow((P[i] - RotationCenter[i]), 2);
    }
    if (abs(sqrt(distance)-Radii)<=2*CS*TimeStep()){
        seriesnum = -1;
    }
    if (distance - Radii * Radii < 0.0-1e-13) {
        return -1;
    }else {
        return 2;
    }
}

static inline OPS_FUN_PREFIX void CalcIntersectIBB(const int& cx, const int& cy,
                                                   const Real& x, const Real& y,
                                                   const Real* RotationCenter,
                                                   const Real& Radii,
                                                   Real* Root, Real& Distance) {
    const Real centerX = RotationCenter[0];
    const Real centerY = RotationCenter[1];
    Real MeshSize{CS * TimeStep()};
    if (abs((x + cx * MeshSize - centerX) * (x + cx * MeshSize - centerX) +
            (y + cy * MeshSize - centerY) * (y + cy * MeshSize - centerY) -
            Radii * Radii) <= 1e-13) {
        Root[0] = x + cx * MeshSize;
        Root[1] = y + cy * MeshSize;
        Distance = 1;
    } else if (abs((x - centerX) * (x - centerX) +
                   (y - centerY) * (y - centerY) - Radii * Radii) <= 1e-13) {
        Root[0] = x;
        Root[1] = y;
        Distance = 0;
    } else {
        if (cx != 0) {
            Real slope = (Real)cy / cx;
            Real intercept = y - slope * x;
            Real a = 1 + slope * slope;
            Real b = 2 * slope * (intercept - centerY) - 2 * centerX;
            Real c = centerX * centerX - Radii * Radii +
                     (intercept - centerY) * (intercept - centerY);
            Real result = b * b - 4 * a * c;
            if (result < 0 - 1e-12) {
                    ops_printf("No intersect!\n");
                    assert(result >= 0);
            }
            Root[0] = (-b + sqrt(result)) / (2 * a);
            Root[1] = slope * Root[0] + intercept;
            Distance = sqrt((Root[0] - x) * (Root[0] - x) +
                            (Root[1] - y) * (Root[1] - y)) /
                       (sqrt((Real)cx * cx + (Real)cy * cy) * MeshSize) *
                       (Root[0] - x) * cx / (abs((Root[0] - x) * cx) + 1e-16);
            // ops_printf("No intersect1! %2.16f\n", Distance);
            if (Distance >= 1 + 1e-14 || Distance < 0 - 1e-14) {
                    Real X = (-b - sqrt(result)) / (2 * a);
                    Real Y = slope * X + intercept;
                    Real D = sqrt((X - x) * (X - x) + (Y - y) * (Y - y)) /
                             (sqrt((Real)cx * cx + (Real)cy * cy) * MeshSize) *
                             (X - x) * cx / (abs((X - x) * cx) + 1e-16);
                    if (abs(D) < abs(Distance)) {
                        Root[0] = X;
                        Root[1] = Y;
                        Distance = D;
                    }
            }
        } else if (cx == 0) {
            Root[0] = x;
            Root[1] =
                sqrt(Radii * Radii - (x - centerX) * (x - centerX)) + centerY;
            Distance = (Root[1] - y) / (cy * (MeshSize));
            // ops_printf("No intersect1! %2.16f\n", Distance);
            if (Distance >= 1 + 1e-14 || Distance < 0 - 1e-14) {
                    Real Y =
                        -sqrt(Radii * Radii - (x - centerX) * (x - centerX)) +
                        centerY;
                    Real D = (Y - y) / (cy * (MeshSize));
                    if (abs(D) < abs(Distance)) {
                        Root[1] = Y;
                        Distance = D;
                    }
                    // ops_printf("No intersect2! %e\n", (Distance - 1));
            }
            // ops_printf("No intersect! %d,%d,%e,%e,%e\n",cx,cy,x,y,MeshSize);
        }
    }
}

// static inline OPS_FUN_PREFIX void CalcIntersectIBB(const Real& cx, const Real& cy,
//                                                 const Real& x, const Real& y,
//                                                 const Real* RotationCenter,
//                                                 const Real& Radii, Real& RootX,
//                                                 Real& RootY, Real& Distance) {
//     const Real centerX = RotationCenter[0];
//     const Real centerY = RotationCenter[1];
//     if (abs((x + cx - centerX) * (x + cx - centerX) +
//             (y + cy - centerY) * (y + cy - centerY) - Radii * Radii) <= 1e-14) {
//        Root[0] = x + cx;
//         RootY = y + cy;
//         Distance = 1;
//         // ops_printf("No intersect2! %e, %e, %e, %e, %e\n", (Distance - 1) *
//         // 1e10,
//         //            cx, cy, x, y);
//     } else if (abs((x - centerX) * (x - centerX) +
//                    (y - centerY) * (y - centerY) - Radii * Radii) <= 1e-14) {
//         RootX = x;
//         RootY = y;
//         Distance = 0;
//     } else {
//         if (cx > 1e-14 || cx < -1e-14) {
//             Real slope = (Real)cy / cx;
//             Real intercept = y - slope * x;
//             Real a = 1 + slope * slope;
//             Real b = 2 * slope * (intercept - centerY) - 2 * centerX;
//             Real c = centerX * centerX - Radii * Radii +
//                      (intercept - centerY) * (intercept - centerY);
//             Real result = b * b - 4 * a * c;
//             if (result < 0 - 1e-14) {
//                     ops_printf("No intersect! %e, %e, %e, %e, %e\n", result, cx,
//                                cy, x, y);
//                     assert(result >= 0);
//             }
//             RootX = (-b + sqrt(result)) / (2 * a);
//             RootY = slope * RootX + intercept;
//             Distance =
//                 sqrt((RootX - x) * (RootX - x) + (RootY - y) * (RootY - y)) /
//                 (sqrt((Real)cx * cx + (Real)cy * cy));
//             assert(result >= 0);
//             if (Distance >= 1 + 1e-14) {
//                     RootX = (-b - sqrt(result)) / (2 * a);
//                     RootY = slope * RootX + intercept;
//                     Distance = sqrt((RootX - x) * (RootX - x) +
//                                     (RootY - y) * (RootY - y)) /
//                                (sqrt((Real)cx * cx + (Real)cy * cy));
//             }
//         } else {
//             RootX = x;
//             RootY =
//                 sqrt(Radii * Radii - (x - centerX) * (x - centerX)) + centerY;
//             Distance = sqrt((RootY - y) * (RootY - y)) / (sqrt((Real)cy * cy));
//             // ops_printf("term %e, %e, %e,%e ,%e\n",Distance,x,y,RootX,
//             //                   RootY);
//             if (Distance >= 1 + 1e-14) {
//                     RootY =
//                         -sqrt(Radii * Radii - (x - centerX) * (x - centerX)) +
//                         centerY;
//                     Distance =
//                         sqrt((RootY - y) * (RootY - y)) / (sqrt((Real)cy * cy));
//             }
//         }
//     }
// }

static inline OPS_FUN_PREFIX void CalcIntersectIBB(const Real& cx, const Real& cy,
                                                const Real& x, const Real& y,
                                                const Real* RotationCenter,
                                                const Real& Radii,
                                                Real* Root, Real& Distance) {
    const Real centerX = RotationCenter[0];
    const Real centerY = RotationCenter[1];
    if (abs((x + cx - centerX) * (x + cx - centerX) +
            (y + cy - centerY) * (y + cy - centerY) - Radii * Radii) <= 1e-13) {
        Root[0] = x + cx;
        Root[1] = y + cy;
        Distance = 1;
        // ops_printf("No intersect2! %e, %e, %e, %e, %e\n", (Distance - 1) * 1e10,
        //            cx, cy, x, y);
    } else if (abs((x  - centerX) * (x - centerX) +
            (y - centerY) * (y - centerY) - Radii * Radii) <= 1e-13) {
        Root[0] = x;
        Root[1] = y;
        Distance = 0;
    } else {
        if (abs(cx)>1e-10) {
            Real slope = (Real)cy / cx;
            Real intercept = y - slope * x;
            Real a = 1 + slope * slope;
            Real b = 2 * slope * (intercept - centerY) - 2 * centerX;
            Real c = centerX * centerX - Radii * Radii +
                     (intercept - centerY) * (intercept - centerY);
            Real result = b * b - 4 * a * c;
            if (result < 0 - 1e-14) {
                    ops_printf("No intersect! %e, %e, %e, %e, %e\n", result, cx,
                               cy, x, y);
                    assert(result >= 0);
            }
            Root[0] = (-b + sqrt(result)) / (2 * a);
            Root[1] = slope * Root[0] + intercept;
            Distance =
                sqrt((Root[0] - x) * (Root[0] - x) + (Root[1] - y) * (Root[1] - y)) /
                (sqrt((Real)cx * cx + (Real)cy * cy))*
                (Root[0] - x)*cx/(abs((Root[0] - x)*cx)+1e-16);
            assert(result >= 0);
            if (Distance >= 1 + 1e-14||Distance <0 - 1e-14){
                    Real X = (-b - sqrt(result)) / (2 * a);
                    Real Y = slope * X + intercept;
                    Real D = sqrt((X - x) * (X - x) + (Y - y) * (Y - y)) /
                             (sqrt((Real)cx * cx + (Real)cy * cy))*
                             (X - x)*cx/(abs((X - x)*cx)+1e-16);
                    if (abs(D) < abs(Distance)) {
                        Root[0] = X;
                        Root[1] = Y;
                        Distance = D;
                    }
            }
        } else {
            Root[0] = x;
            Root[1] =
                sqrt(Radii * Radii - (x - centerX) * (x - centerX)) + centerY;
            Distance = (Root[1] - y) / cy;
            // ops_printf("term %e, %e, %e,%e ,%e\n",Distance,x,y,Root[0],
            //                   Root[1]);
            if (Distance >= 1 + 1e-14||Distance <0 - 1e-14){
                    Real Y =
                        -sqrt(Radii * Radii - (x - centerX) * (x - centerX)) +
                        centerY;
                    Real D = (Y - y) / cy;

                    if (abs(D) < abs(Distance)) {
                        Root[1] = Y;
                        Distance = D;
                    }
            }
        }
    }
}

static inline OPS_FUN_PREFIX void CalcIntersectFree(const Real* NormalDirection,
                                                const Real& VOF,
                                                const int& cx, const int& cy,
                                                const Real& x, const Real& y,
                                                Real* Root,
                                                bool& Subvof) {
    Real PlaneConst{0};
    Real MeshSize{CS * TimeStep()};
    Real i = (x ) / MeshSize;
    //BlockCoor[0] denotes the end of x coordinate, not y
    Real j = (y ) / MeshSize;
    // determine the free surface plane constant
    Real nx{abs(NormalDirection[0])};
    Real ny{abs(NormalDirection[1])};
    const Real normalization{abs(nx)+abs(ny)};
    Real m{0};
    Real vof{0};
    // bring to the standard case 0<m1<m2 in Reference
    // doi:10.1006/jcph.2000.6567
    if (nx <= 0 + 1e-14 && ny <= 0+1e-14) {
        return;
    }
    if (nx < ny) {
        m = nx / normalization;
    } else {
        m = ny / normalization;
    }
    if (m >= 0 + 1e-14) {
        if( VOF >= 1 || VOF <= 0) {
            vof = 0;
        } else if (VOF > 0.5){
            vof = 1 - VOF;
        } else {
            vof = VOF;
        }
        const Real Criticalvof{m / (2 * (1 - m))};
        if (vof <= Criticalvof) {
            if (vof <= 0 + 1e-12) {
                    PlaneConst = 0;
            } else {
                    PlaneConst = sqrt(2 * m * (1 - m) * vof);
            }

        } else if (vof > Criticalvof && vof <= 0.5) {
            PlaneConst = vof * (1 - m) + 0.5 * m;
        }
        if (VOF > 0.5) {
            PlaneConst = 1 - PlaneConst;
        }
    } else if (m < 0 + 1e-14) {
        PlaneConst = VOF;
    }  // for nx=0 or ny=0
    if (isnan(PlaneConst)|| isinf(PlaneConst)) {
        ops_printf("Error,  PlaneConst =%e \n", PlaneConst);
        assert(!(isnan(PlaneConst) || isinf(PlaneConst)));
    }
    // if (vof>0) {
    //     printf("Error,  PlaneConst =%e \n", PlaneConst);
    // }
    //************************//
    // intersection between free surface plane and D2Q9 lattice direction
    Real Nx{NormalDirection[0]};
    Real Ny{NormalDirection[1]};
    if (abs(cx * Nx + cy * Ny) > 0 + 1e-14) {
        Real C = PlaneConst * normalization + Ny * (j - 0.5) +
                 Nx * (i - 0.5);  // free surface parameter
        if (Nx < 0) {
            C += Nx;
        }
        if (Ny < 0) {
            C += Ny;
        }
        if (cx != 0) {
            Real slope = (Real)cy / cx;
            Real intercept = j - slope * i;  // lattice direction D2Q9
            Root[0] = (C - Ny * intercept) / (slope * Ny - Nx * (-1));
            Root[1] = slope * Root[0] + intercept;
        } else {
            Root[0] = (Real)i;
            Root[1] = (C - Nx * i) / Ny;
        }
        if (vof <= 0.01) {
            Root[0] = (Real)i + cx * (VOF - 0.5);
            Root[1] = (Real)j + cy * (VOF - 0.5);
        }// to reduce the error caused by wrong normal directions
    }
    if (Subvof){
        //can further choose another ortho direction
         if ((Root[0]-i)*Nx+(Root[1]-j)*Ny>=0){
            Subvof = true;
         }else{
            Subvof = false;
         }
    }
   
}

static inline OPS_FUN_PREFIX Real VectorProduct(const Real* p1, const Real* p2,
                                                const Real* p) {
    Real pp1[]{p1[0] - p[0], p1[1] - p[1]};
    Real pp2[]{p2[0] - p[0], p2[1] - p[1]};
    return pp1[0] * pp2[1] - pp1[1] * pp2[0];
}

static inline OPS_FUN_PREFIX Real ScalarProduct(const Real* p1, const Real* p2,
                                                const Real* p) {
    Real pp1[]{p1[0] - p[0], p1[1] - p[1]};
    Real pp2[]{p2[0] - p[0], p2[1] - p[1]};
    return pp1[0] * pp2[0] + pp1[1] * pp2[1];
}

static inline OPS_FUN_PREFIX Real ScalarProduct(const Real* p1, const Real* p2,
                                                const Real* p,const Real* n) {
    Real pp1[]{p1[0] - p[0], p1[1] - p[1]};
    Real pp2[]{p2[0] - p[0], p2[1] - p[1]};
    return ((pp1[0] * n[0] + pp1[1] * n[1]) + (pp2[0] * n[0] + pp2[1] * n[1]))/2;
}


static inline OPS_FUN_PREFIX Real CalcSolidVOF(const Real* RotationCenter,
                                               const Real& Radii, const Real* P) {
    const Real x{P[0]};
    const Real y{P[1]};
    const Real CentralX{RotationCenter[0]};
    const Real CentralY{RotationCenter[1]};
    Real MeshSize{CS * TimeStep()};
    Real halfMesh{0.5 * MeshSize};
    Real dy_down = abs(y - CentralY) - halfMesh;
    Real dy_up = abs(y - CentralY) + halfMesh;
    Real dx_down = abs(x - CentralX) - halfMesh;
    Real dx_up = abs(x - CentralX) + halfMesh;
    // the edge closer to center is defined as the down
    Real XIntersectdy_up{1000 * Radii};
    Real XIntersectdy_down{1000 * Radii};
    Real YIntersectdx_up{1000 * Radii};
    Real YIntersectdx_down{1000 * Radii};
    Real ChordLength{0};
    Real Angle{0};
    Real AreaPolygon{0};
    Real AreaofBow{0};
    Real vof{0};
    // angle(0, 0) = 0;
    // ForcePoints(0,0,0)=0;
    // ForcePoints(1,0,0)=0;

    if ((x - CentralX) * (x - CentralX) + 
        (y - CentralY) * (y - CentralY) <=-1e-14+
        Radii * Radii) {
         vof = 1;
         // maybe the position can be emerged into geometry
    }
    if (abs(dy_down) <= Radii) {
        XIntersectdy_down = sqrt(Radii * Radii - dy_down * dy_down);
    }
    if (abs(dy_up) <= Radii) {
        XIntersectdy_up = sqrt(Radii * Radii - dy_up * dy_up);
    }
    if (abs(dx_down) <= Radii) {
        YIntersectdx_down = sqrt(Radii * Radii - dx_down * dx_down);
    }
    if (abs(dx_up) <= Radii) {
        YIntersectdx_up = sqrt(Radii * Radii - dx_up * dx_up);
    }
    if (XIntersectdy_down > dx_down && XIntersectdy_down < dx_up) {
        if (XIntersectdy_up > dx_down && XIntersectdy_up < dx_up) {
            ChordLength = sqrt((XIntersectdy_up - XIntersectdy_down) *
                               (XIntersectdy_up - XIntersectdy_down) +
                               MeshSize * MeshSize);
            AreaPolygon = (XIntersectdy_up + XIntersectdy_down 
                            - 2 * dx_down) * MeshSize / 2;
            // ForcePoints(0, 0, 0) = (XIntersectdy_up + XIntersectdy_down) / 2;
            // ForcePoints(1, 0, 0) = (dy_up + dy_down) / 2;
        } else if (YIntersectdx_down > dy_down && 
                   YIntersectdx_down < dy_up) {
            ChordLength = sqrt(
                        (YIntersectdx_down - dy_down) * 
                        (YIntersectdx_down - dy_down) +
                        (XIntersectdy_down - dx_down) * 
                        (XIntersectdy_down - dx_down));
            AreaPolygon = ((YIntersectdx_down - dy_down) *
                           (XIntersectdy_down - dx_down)) /2;
            // ForcePoints(0, 0, 0) = (XIntersectdy_down + dx_down) / 2;
            // ForcePoints(1, 0, 0) = (YIntersectdx_down + dy_down) / 2;
        }
        Angle = acos(1 - ChordLength * ChordLength / (2 * Radii * Radii));
        AreaofBow = (Angle - sin(Angle)) * Radii * Radii / 2;
        vof = (AreaPolygon + AreaofBow) / (MeshSize * MeshSize);
        // angle(0,0) =Angle;
    } else if (YIntersectdx_up > dy_down && YIntersectdx_up < dy_up) {
        if (YIntersectdx_down > dy_down && YIntersectdx_down < dy_up) {
            ChordLength = sqrt((YIntersectdx_up - YIntersectdx_down) *
                               (YIntersectdx_up - YIntersectdx_down) +
                               MeshSize * MeshSize);
            AreaPolygon = (YIntersectdx_up + YIntersectdx_down 
                            - 2 * dy_down) * MeshSize / 2;
            // ForcePoints(0, 0, 0) = (dx_down + dx_up) / 2;
            // ForcePoints(1, 0, 0) = (YIntersectdx_up + YIntersectdx_down) / 2;
        } else if (XIntersectdy_up > dx_down && XIntersectdy_up < dx_up) {
            ChordLength =sqrt(
                        (YIntersectdx_up - dy_up) * 
                        (YIntersectdx_up - dy_up) +
                        (XIntersectdy_up - dx_up) * 
                        (XIntersectdy_up - dx_up));
            AreaPolygon = MeshSize * MeshSize -
                            ((YIntersectdx_up - dy_up) * 
                             (XIntersectdy_up - dx_up)) / 2;
            // ForcePoints(0, 0, 0) = (XIntersectdy_up + dx_up) / 2;
            // ForcePoints(1, 0, 0) = (YIntersectdx_up + dy_up) / 2;
        }
        Angle = acos(1 - ChordLength * ChordLength / (2 * Radii * Radii));
        AreaofBow = (Angle - sin(Angle)) * Radii * Radii / 2;
        vof = (AreaPolygon + AreaofBow) / 
                            (MeshSize * MeshSize);
        // angle(0, 0) = Angle;
    }
    if (dx_down < 0) {
        if (XIntersectdy_down < abs(dx_down)) {
            ChordLength = XIntersectdy_down * 2;
            Angle = acos(1 - ChordLength * ChordLength 
                                / (2 * Radii * Radii));
            AreaofBow = (Angle - sin(Angle)) * Radii * Radii / 2;
            vof = AreaofBow / (MeshSize * MeshSize);
            
        } else if (XIntersectdy_up < abs(dx_down)) {
            ChordLength = XIntersectdy_up * 2;
            Angle = acos(1 - ChordLength * ChordLength 
                                / (2 * Radii * Radii));
            AreaofBow = (Angle - sin(Angle)) * Radii * Radii / 2;
            Real ChordLength2 = 
                        sqrt((YIntersectdx_up - YIntersectdx_down) *
                             (YIntersectdx_up - YIntersectdx_down) +
                             MeshSize * MeshSize);
            Real Angle2 =
                acos(1 - ChordLength2 * ChordLength2 
                              / (2 * Radii * Radii));
            Real AreaofBow2 = 
                    (Angle2 - sin(Angle2)) * Radii * Radii / 2;
            AreaPolygon = (YIntersectdx_up + YIntersectdx_down 
                            - 2 * dy_down) * MeshSize / 2;
            vof =
                (AreaPolygon + AreaofBow2 - AreaofBow) / 
                (MeshSize * MeshSize);
            // angle(0,0) =Angle2;
            // ForcePoints(0, 0, 0) = 0 ;
            // ForcePoints(1, 0, 0) = dy_up;
        }
    }
    if (dy_down < 0) {
        if (YIntersectdx_down < abs(dy_down)) {
            ChordLength = YIntersectdx_down * 2;
            Angle = acos(1 - ChordLength * ChordLength 
                                / (2 * Radii * Radii));
            AreaofBow = (Angle - sin(Angle)) * Radii * Radii / 2;
            vof = AreaofBow / (MeshSize * MeshSize);
        } else if (YIntersectdx_up < abs(dy_down)) {
            ChordLength = YIntersectdx_up * 2;
            Angle = acos(1 - ChordLength * ChordLength 
                                / (2 * Radii * Radii));
            AreaofBow = (Angle - sin(Angle)) * Radii * Radii / 2;
            Real ChordLength2 = 
                    sqrt((XIntersectdy_up - XIntersectdy_down) *
                         (XIntersectdy_up - XIntersectdy_down) +
                         MeshSize * MeshSize);
            Real Angle2 = acos(1 - ChordLength2 * ChordLength2 / 
                                            (2 * Radii * Radii));
            Real AreaofBow2 = 
                    (Angle2 - sin(Angle2)) * Radii * Radii / 2;
            AreaPolygon = (XIntersectdy_up + XIntersectdy_down 
                            - 2 * dx_down) * MeshSize / 2;
            vof = (AreaPolygon + AreaofBow2 - AreaofBow) 
                            / (MeshSize * MeshSize);
            // angle(0, 0) = Angle2;
            // ForcePoints(0, 0, 0) = dx_up;
            // ForcePoints(1, 0, 0) = 0;
        }
    }
    return vof;
}

//ray method, preparing for 3D case if on boundary plane
//TODO need to further improve
static inline OPS_FUN_PREFIX bool InPolygon(const Real* polygonVertex,
                                            const int number,
                                            const Real* P) {
    // const int number{sizeof(XpolygonVertex)/sizeof(XpolygonVertex[0])};
    // ops_printf("number of segments...%d\n",number);
    // const Real LargeP[]{P[0]+1000,P[1]};
    bool InPolygon{false};
    const int dim{SpaceDim()};
    for (int i = 0; i <= number; i++) {
        int j = i+1;
        if (number == i + 1) {
            j=0;
        }
        Real PointO[]{polygonVertex[i*dim],polygonVertex[i*dim+1]};
        Real PointE[]{polygonVertex[j*dim],polygonVertex[j*dim+1]};
        Real VP{VectorProduct(PointO,PointE,P)};
        Real SP{ScalarProduct(PointO,PointE,P)};
        if(abs(VP)<=1e-14&&SP<=0){
            return true;
        }//include the cases that P at solid vertex and on segment
        if (abs(PointE[1]-PointO[1])<=1e-14){
            continue;
        }
        if ((PointO[1]-P[1]>1e-14)!=(PointE[1]-P[1]>1e-14)){
            // the point X between Point0 and PointE
            if (VectorProduct(P, PointE, PointO) / (PointE[1] - PointO[1]) <
                0) {
                    InPolygon = !InPolygon;
            }
        } else if ((abs(PointO[1] - P[1]) <= 1e-14) &&
                   (PointO[1] - PointE[1] > 1e-14)) {
            InPolygon = !InPolygon;
        } else if ((abs(PointE[1] - P[1]) <= 1e-14) &&
                   (PointE[1] - PointO[1] > 1e-14)) {
            InPolygon = !InPolygon;
            // consider the ray line intersects with solid vertex
        }
    }
    return InPolygon;
}

static inline OPS_FUN_PREFIX void CalcSegmentNo(int& seriesnum, const Real* P,
                                            const Real* polygonVertex,
                                            const Real* PlaneNormal,
                                            const int polygon,
                                            const int NumofPol) {
    // ensure the normal is unit vector
    // const int number{sizeof(XpolygonVertex)};
    const int dim{SpaceDim()};
    for (int num = 0; num < NumofPol; num++) {
        int numofNearestSegement{0};
        Real distancefromsegmentold{100};
        for (int i = 0; i < polygon; i++) {
            // make the last vertex connect to the first node
            int j = i + 1;
            int I = i + num * polygon;
            int J = j + num * polygon;
            if (polygon == i + 1) {
                    J = num * polygon;
            }
            Real PointO[]{polygonVertex[I * dim], polygonVertex[I * dim + 1]};
            Real PointE[]{polygonVertex[J * dim], polygonVertex[J * dim + 1]};
            Real Normal[]{PlaneNormal[I * dim], PlaneNormal[I * dim + 1]};
            Real SP{ScalarProduct(PointO, PointE, P)};
            Real distancefromsegment{ScalarProduct(PointO, PointE, P, Normal)};
            // note:the distancefromsegment is signed.
            if (abs(distancefromsegment) <= sqrt(2) * CS * TimeStep()) {
                // if(abs(distancefromsegment) <= 0.5 *sqrt(2) * CS * TimeStep()){
                    numofNearestSegement++;
                // }
                if (SP <= 0) {
                    if (abs(distancefromsegment) < distancefromsegmentold) {
                        distancefromsegmentold = abs(distancefromsegment);
                        seriesnum = I;
                    }
                }
            }
        }
        if (numofNearestSegement > 1) {
            seriesnum = -1;
            // denote the case that node at corner.
        }
        if (num >= 1 && (num + 1) * polygon - 1 == seriesnum) {
            seriesnum =-10;
        }
    }
}


static inline OPS_FUN_PREFIX int InPolygon(const int& seriesnum,
                                            const Real* P,
                                            const Real* polygonVertex,
                                            const Real* PlaneNormal,
                                            const int polygon,
                                            const int NumofPol,
                                            Real& d,
                                            const int subgrids=1) {
    //***********************************************************//
    //Warning: the VP and distancefromsegment in function InPolygon()
    //must be less tolerance than those in function CalcIntersectIBB() 
    //distancefromsegment <= 1e-14 Vs. distancefromsegment <= 1e-12 in latter function                       
    // ensure the normal is unit vector
    // const int number{sizeof(XpolygonVertex)};
    const int dim{SpaceDim()};
    if (seriesnum>=0){
        int I = seriesnum;
        int J = seriesnum + 1;
        // make the last vertex connect to the first node
        if ((seriesnum + 1) % polygon == 0) {
            J = seriesnum + 1 - polygon;
        }
        Real PointO[]{polygonVertex[I * dim], polygonVertex[I * dim + 1]};
        Real PointE[]{polygonVertex[J * dim], polygonVertex[J * dim + 1]};
        Real Normal[]{PlaneNormal[I * dim], PlaneNormal[I * dim + 1]};
        Real VP{VectorProduct(PointO, PointE, P)};
        Real SP{ScalarProduct(PointO, PointE, P)};
        Real distancefromsegment{ScalarProduct(PointO, PointE, P, Normal)};
        // note:the distancefromsegment is signed.
        if (abs(VP) <= 1e-14 && SP <= 0) {
            if(seriesnum >= polygon){
                return -1;
            }else{
                return -1;
            }
        }
        if ((distancefromsegment > 0+1e-14 && seriesnum < polygon) ||
            (distancefromsegment <= 0-1e-14 && seriesnum >= polygon)) {
            return -1;//in solid
        } else if(distancefromsegment <= 0-1e-14 && seriesnum < polygon){
            return 2;//out chamber
        } else {
            return 1;//in chamber
        }
        
    }else {
        bool onSegmentof1stSolid{false};
        for (int num = 0; num < NumofPol; num++) {
            int numofInPolygon{0};
            // Real distancefromsegmentold{100};
            for (int i = 0; i < polygon; i++) {
                // make the last vertex connect to the first node
                    int j = i + 1;
                    int I = i + num * polygon;
                    int J = j + num * polygon;
                    if (polygon == i + 1) {
                        J = num * polygon;
                    }
                Real PointO[]{polygonVertex[I * dim],
                              polygonVertex[I * dim + 1]};
                Real PointE[]{polygonVertex[J * dim],
                              polygonVertex[J * dim + 1]};
                Real Normal[]{PlaneNormal[I * dim], PlaneNormal[I * dim + 1]};
                Real VP{VectorProduct(PointO, PointE, P)};
                Real SP{ScalarProduct(PointO, PointE, P)};
                Real distancefromsegment{
                    ScalarProduct(PointO, PointE, P, Normal)};
                // note:the distancefromsegment is signed.
                // if (abs(VP) <= 1e-15 && SP <= 0) {
                //     d = 0;
                //     if (polygon == i + 1) {
                //         if (num < 1) {
                //             onSegmentof1stSolid = true;
                //             continue;
                //         } else {
                //             if (SP < 0-1e-14) {
                //                 d = 1;
                //                 return 1;
                //             }
                //         }
                //     }
                //     return 2;    
                // } 
                if (abs(VP) <= 1e-14 && SP <= 0 &&
                    polygon == i + 1 && num == 0&&0!=NumofPol-1) {
                            onSegmentof1stSolid = true;
                            continue;
                }
                if (onSegmentof1stSolid&&polygon == i + 1
                    &&num >= 1&&abs(VP) <= 1e-14){
                    if (SP>0+1e-14){
                        d=0;
                        return -1;
                    }
                }
                if (abs(VP) <= 1e-14 && SP <= 0+1e-14) {
                        if (num >= 1) {
                            d=0;
                            if (polygon == i + 1) {
                                d=1;
                                return 1;

                            }
                            return -1;
                        }else{
                            d = 0;
                            return -1;
                        } 
                }

                if (distancefromsegment > 0+1e-14 ) {
                    numofInPolygon++;
                }
                if (polygon == i + 1) {
                    for (int snum = 1; snum < NumofPol; snum++) {
                        Real sPointO[]{polygonVertex[snum * polygon * dim],
                                       polygonVertex[snum * polygon * dim + 1]};
                        Real sPointE[]{
                            polygonVertex[((snum + 1) * polygon - 1) * dim],
                            polygonVertex[((snum + 1) * polygon - 1) * dim +
                                          1]};
                        if (ScalarProduct(sPointO, sPointE, P) < 0 ) {
                            distancefromsegment = 
                                0.5 * CS * TimeStep() / subgrids;
                        }
                    }
                }
                if (abs(distancefromsegment)/ 
                    (0.5 * CS * TimeStep() / subgrids) < d&&SP <= 0.0+1e-14) {
                    d = abs(distancefromsegment)/ 
                        (0.5 * CS * TimeStep() / subgrids);
                }
            }
            if (onSegmentof1stSolid){
                continue;
            }
            if (d > 1) {
                d = 1;
            }
            if ((num < 1 && numofInPolygon != polygon) ) {
                return 2;
            } else if(num >= 1 && numofInPolygon == polygon){
                return 1;
            }
        }
        return -1;
    }
}

static inline OPS_FUN_PREFIX void CalcIntersectIBB(const int& seriesnum,
                                            const Real& cx, const Real& cy,
                                            const Real* P,
                                            const Real* polygonVertex,
                                            const Real* PlaneNormal,
                                            const int polygon,
                                            const int NumofPol,
                                            Real* Root,
                                            Real& Distance) {
    // Real MeshSize{CS * TimeStep()};
    const int dim{SpaceDim()};
    Real Pneighbor[2]{0};
    Pneighbor[0]=P[0]+cx;
    Pneighbor[1]=P[1]+cy;
    if (seriesnum>=0){
        int I = seriesnum;
        int J = seriesnum + 1;
        // make the last vertex connect to the first node
        if ((seriesnum + 1) % polygon == 0) {
            J = seriesnum + 1 - polygon;
        }
        Real PointO[]{polygonVertex[I * dim], polygonVertex[I * dim + 1]};
        Real PointE[]{polygonVertex[J * dim], polygonVertex[J * dim + 1]};
        Real Normal[]{PlaneNormal[I * dim], PlaneNormal[I * dim + 1]};
        Real NormalDotLatticeVel{cx*Normal[0]+cy*Normal[1]};
        Real distancefromsegment{
                    ScalarProduct(PointO, PointE, P, Normal)};
        if (abs(distancefromsegment)<=1e-12) {
                Distance = 0;
                Root[0] = P[0];
                Root[1] = P[1];
        } else if(abs(ScalarProduct(PointO, PointE, Pneighbor, Normal))<=1e-12){
                Distance = 1;
                Root[0] = P[0] + cx;
                Root[1] = P[1] + cy;
        } else {
        if (abs(NormalDotLatticeVel)>1e-15){
            //do not perpendicular to each other n*(cx,cy)
            Real d = distancefromsegment /
                     NormalDotLatticeVel;
            // ops_printf("d1 = %e,seriesnum = %d, %e, %e\n",(d-1)*1e10,seriesnum,P[0],P[1]);
            if (d <= 1+1e-10  && d >= 0.0-1e-10) {
                Distance = d;  // check d>0
                Root[0] = P[0] + d * cx;
                Root[1] = P[1] + d * cy;
            }
        }//else{
        // }
        if (ScalarProduct(PointO, PointE, Root)>=0){
            // printf("d1 = %e,seriesnum = %d, %e, %e\n",(d-1)*1e10,seriesnum,P[0],P[1]);
            printf(
                    "Warning,there is almost no intersect,"
                    "segmentNo = %d, x =%e, y =%e, Rootx =%e,Rooty =%e, cx =%e, cy =%e\n",
                     seriesnum, P[0], P[1], Root[0],Root[1], cx, cy );
            assert(ScalarProduct(PointO, PointE, Root)<=1e-14);
        }
        }
    } else {
        for (int num = 0; num < NumofPol; num++) {
            
            for (int i = 0; i < polygon; i++) {
                if (num>=1&&polygon-1==i){
                    continue;
                }
                // make the last vertex connect to the first node
                int j = i+1;
                int I = i + num * polygon;
                int J = j + num * polygon;
                if (polygon == i + 1) {
                    J = num * polygon;
                }
                Real Intersect[SpaceDim()]{-100, -100};
                Real distance{-100};
                Real PointO[]{polygonVertex[I * dim],
                              polygonVertex[I * dim + 1]};
                Real PointE[]{polygonVertex[J * dim],
                              polygonVertex[J * dim + 1]};
                Real Normal[]{PlaneNormal[I * dim], PlaneNormal[I * dim + 1]};
                Real NormalDotLatticeVel{cx * Normal[0] + cy * Normal[1]};
                Real distancefromsegment{
                    ScalarProduct(PointO, PointE, P, Normal)};
                if (abs(ScalarProduct(PointO, PointE, Pneighbor, Normal))<=1e-12&&
                            ScalarProduct(PointO, PointE, Pneighbor)<=0){
                    distance = 1;
                    Intersect[0] = P[0] + cx;
                    Intersect[1] = P[1] + cy;
                } else if(abs(distancefromsegment)<=1e-12&&ScalarProduct(PointO, PointE, P)<=0) {
                    distance = 0;
                    Intersect[0] = P[0];
                    Intersect[1] = P[1];
                } else {    
                if (abs(NormalDotLatticeVel) > 1e-15) {
                    // ensure to not be parellel each other
                    Real d = distancefromsegment /
                             NormalDotLatticeVel;
                    // assert(d >= -0.0);
                    if (d <= 1+1e-10 && d >= 0.0-1e-10) {
                        distance = d;
                        Intersect[0] = P[0] + d * cx;
                        Intersect[1] = P[1] + d * cy;
                    } 
                } //else {
                //     Real VP{VectorProduct(PointO, PointE, P)};
                //     Real SP{ScalarProduct(PointO, PointE, P)};
                //     if (abs(VP) <= 1e-14) {
                //         if (SP <= 0) {
                //             distance = 0;
                //             Intersect[0] = P[0];
                //             Intersect[1] = P[1];
                //         } else {
                            
                //             distance =
                //                 sqrt((PointO[0] - P[0]) * (PointO[0] - P[0]) +
                //                      (PointO[1] - P[1]) * (PointO[1] - P[1])) /
                //                 (cx * cx + cy * cy);
                //             if (distance > 1 + 1e-14) {
                //                 distance = sqrt((PointE[0] - P[0]) *
                //                                     (PointE[0] - P[0]) +
                //                                 (PointE[1] - P[1]) *
                //                                     (PointE[1] - P[1])) /
                //                            (cx * cx + cy * cy);
                //             }
                //             Intersect[0] = P[0] + cx;
                //             Intersect[1] = P[1] + cy;
                //         }
                //     }
                // }
                }
                if (polygon-1==i){
                for (int num = 1; num < NumofPol; num++) {
                    int I = polygon - 1 + num * polygon;
                    int J = num * polygon;
                    // Real Intersect[SpaceDim()]{-100, -100};
                    // Real distance{-100};
                    Real PointO[]{polygonVertex[I * dim],
                              polygonVertex[I * dim + 1]};
                    Real PointE[]{polygonVertex[J * dim],
                              polygonVertex[J * dim + 1]};
                    Real Normal[]{PlaneNormal[I * dim], PlaneNormal[I * dim + 1]};
                    Real distancefromsegment{
                        ScalarProduct(PointO, PointE, Intersect, Normal)};
                    if (abs(distancefromsegment)<=1e-12&&
                            ScalarProduct(PointO, PointE, Intersect)<=0){
                        distance = -100;
                        Intersect[0] = -100;
                        Intersect[1] = -100;
                    } 
                }
                }
                if (ScalarProduct(PointO, PointE, Intersect)<=1e-12&&abs(distance)<=1+1e-10){
                    if (abs(distance) < abs(Distance)) {
                        Distance = distance;
                        Root[0] = Intersect[0];
                        Root[1] = Intersect[1];
                    }
                    // ops_printf("Distance = %e, %e, %e\n",Distance,Root[0],Root[1]);
                }
            //     if (abs(Distance) >= 1 + 1e-10&&num==NumofPol-1&&i==3) {
            //     ops_printf("d2 = %e,seriesnum = %d, Px=%e, Py=%e,cx =%e, cy =%e\n",
            //                (distance - 1), i, P[0], P[1], Normal[0], Normal[1]);
            //     assert(abs(Distance) <= 1 + 1e-10);
            // }
            }
            // if(abs(Distance) <= 1 + 1e-10){
            //     return;
            // }
        }
        if (abs(Distance) >= 1 + 1e-10) {
            printf("d2 = %e,seriesnum = %d, Px=%e, Py=%e,cx =%e, cy =%e\n",
                    (Distance - 1), seriesnum, P[0], P[1], cx,cy);
            }
        assert(abs(Distance) <= 1 + 1e-10);
    }
}

static inline OPS_FUN_PREFIX Real CalcSolidVOF(const Real* polygonVertex,
                                                const Real* PlaneNormal,
                                                const Real* P) {
    const Real MeshSize{CS * TimeStep()};
    const Real Normal[] {abs(PlaneNormal[0]), abs(PlaneNormal[1])};
    const Real normalization{Normal[0]+Normal[1]};
    Real m{0};
    Real vof{0};
    // determine the free surface plane constant
    // known surface function to solve the alpha
    Real PlaneConst =
        (PlaneNormal[0] * (polygonVertex[0] - (P[0] - 0.5 * MeshSize)) +
         PlaneNormal[1] * (polygonVertex[1] - (P[1] - 0.5 * MeshSize))) /
        MeshSize / normalization;
    if (PlaneNormal[0] < 0) {
            PlaneConst += -PlaneNormal[0] / normalization;
    }
    if (PlaneNormal[1] < 0) {
            PlaneConst += -PlaneNormal[1] / normalization;
    }
    if (PlaneConst > 1) {
        vof = 1;
    } else if (PlaneConst < 0) {
        vof = 0;
    } else {
        // bring to the standard case 0<m1<m2 in Reference
        // doi:10.1006/jcph.2000.6567
        if (Normal[0] <= 0 + 1e-14 && Normal[1] <= 0 + 1e-14) {
            ops_printf("Error,  SolidPlaneNormal must not be zero \n");
            assert(false);
        }
        if (Normal[0] < Normal[1]) {
            m = Normal[0] / normalization;
        } else {
            m = Normal[1] / normalization;
        }
        if (m >= 0+ 1e-14) {
            Real pc = PlaneConst;
            if (PlaneConst > 0.5) {
                pc = 1 - PlaneConst;
            } else {
                pc = PlaneConst;
            }
            const Real Criticalvof{m / (2 * (1 - m))};
            if (pc< m) {
                if (pc <= 0 + 1e-14) {
                    vof = 0;
                } else {
                    vof = pc*pc/(2*m*(1-m));
                }
            } else if (pc > m && pc <= 0.5) {
                vof = pc / (1 - m) -Criticalvof;
            }
            if (PlaneConst > 0.5) {
                vof = 1 - vof;
            }
        } else if (m < 0 + 1e-14) {
            vof = PlaneConst;
        }  // for nx=0 or ny=0
        if (isnan(vof) || isinf(vof)) {
            ops_printf("Error,  solidvof =%e \n", vof);
            assert(!(isnan(vof) || isinf(vof)));
        }
    }
    return vof;
}


static inline OPS_FUN_PREFIX PointPosition IfPointInPolygon
                            (const Real* point, const Real* polygon,
                                const long long polyVertexNum) {
    // const int number{sizeof(XpolygonVertex)/sizeof(XpolygonVertex[0])};
    // ops_printf("number of segments...%d\n",number);
    // const Real LargeP[]{P[0]+1000,P[1]};
   
    bool InPolygon{false};
    // bool Onsegment{false};
    const int dim{SpaceDim()};
    const Real epsilon{1e-13};
    for (int i = 0; i < polyVertexNum; i++) {
        int j = i+1;
        if (polyVertexNum == i + 1) {
            j=0;
        }
        Real PointO[]{polygon[i*dim],polygon[i*dim+1]};
        Real PointE[]{polygon[j*dim],polygon[j*dim+1]};
        Real VP{VectorProduct(PointO,PointE,point)};
        Real SP{ScalarProduct(PointO,PointE,point)};
        if ((abs(PointO[0]-point[0])<=epsilon&&
             abs(PointO[1]-point[1])<=epsilon)||
            (abs(PointE[0]-point[0])<=epsilon&&
             abs(PointE[1]-point[1])<=epsilon)){
            return IsVertex;//polygon vertex
        }else if(abs(VP)<=epsilon&&SP<=0){
            // Onsegment=true;
            return RelativelyInteriorToEdge;
            //include the cases that P is on segment
        }
        if (abs(PointE[1]-PointO[1])<=epsilon){
            continue;
        }
        if ((PointO[1]-point[1]>epsilon)!=(PointE[1]-point[1]>epsilon)){
            if (VectorProduct(point, PointE, PointO) / (PointE[1] - PointO[1]) <
                0) {
                    InPolygon = !InPolygon;
            }
             
        } else if ((abs(PointO[1] - point[1]) <= epsilon) &&
                   (PointO[1] - PointE[1] > epsilon)) {
            InPolygon = !InPolygon;
        } else if ((abs(PointE[1] - point[1]) <= epsilon) &&
                   (PointE[1] - PointO[1] > epsilon)) {
            InPolygon = !InPolygon;
            // consider the ray line intersects with solid vertex
        }
    }
    if (InPolygon){
        return StrictlyInterior;
    } else {
        return StrictlyExterior;
    }
}

#endif  //  FREESURFACE_HOST_DEVICE_H
