
#ifndef FREESURFACE_HOST_DEVICE_H
#define FREESURFACE_HOST_DEVICE_H
#ifndef OPS_FUN_PREFIX
#define OPS_FUN_PREFIX
#endif

static inline OPS_FUN_PREFIX Real CalcBGKFeqIncom(const int l, const Real rho, const Real u, const Real v,
                const Real T, const int polyOrder) {
    Real cu{(CS * XI[l * LATTDIM] * u + CS * XI[l * LATTDIM + 1] * v)};
    Real c2{(CS * XI[l * LATTDIM] * CS * XI[l * LATTDIM] +
             CS * XI[l * LATTDIM + 1] * CS * XI[l * LATTDIM + 1])};
    Real cu2{cu * cu};
    Real u2{u * u + v * v};
    Real res = rho + cu + 0.5 * (cu2 - u2 + (T - 1.0) * (c2 - LATTDIM));
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
    return WEIGHTS[l]  * res;
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
                    vof = 2*vofs[cx+1-cx][cy+1]-vofs[cx+1-cx-cx][cy+1];
            } break;
            case VG_IM: {
                    vof = 2*vofs[cx+1-cx][cy+1]-vofs[cx+1-cx-cx][cy+1];
            } break;
            case VG_JP: {
                    vof = 2*vofs[cx+1][cy+1-cy]-vofs[cx+1][cy+1-cy-cy];
            } break;
            case VG_JM: {
                    vof = 2*vofs[cx+1][cy+1-cy]-vofs[cx+1][cy+1-cy-cy];
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
    Real LengthofNormal {sqrt(NormalDirection[0]*NormalDirection[0]+
        NormalDirection[1]*NormalDirection[1])};
    if (LengthofNormal>0){
        NormalDirection[0] = NormalDirection[0]/LengthofNormal;
        NormalDirection[1] = NormalDirection[1]/LengthofNormal;
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
                const Real tau,const Real* f, const int* lattIdx,
                const Real Fx=0.0,const Real Fy=0.0) {
        Real tautotal;
        const Real MeshSize = CS*TimeStep();
        const Real C{0.2*0.2};
        const Real U[LATTDIM]{u,v};
        const Real uF[]{2*Fx*u,Fx*v+Fy*u,Fx*v+Fy*u,2*Fy*v};
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
            NonEquilibriumMoment[i*LATTDIM+j] += -(
                        rho*((i==j? 1.0 : 0.0))+U[i]*U[j]) + uF[l]*TimeStep()/2.0;
            NonEquilibriumMoment2 +=
                2 * NonEquilibriumMoment[l] * NonEquilibriumMoment[l];
        }
        tautotal = 0.5 * (tau - 0.5 * TimeStep() +
                     sqrt((tau + 0.5 * TimeStep()) * (tau + 0.5 * TimeStep()) +
                          2 * C * MeshSize* MeshSize *
                              sqrt(NonEquilibriumMoment2)/ rho));
        // printf("\nTotal Wall time %e\n", tautotal);
    return tautotal;
}

static inline OPS_FUN_PREFIX Real CalcIntersectFree(const Real* NormalDirection,
                                                const Real& VOF,
                                                const int& cx, const int& cy) {
    Real Root[4]{};
    Real PlaneConst{0};
    // determine the free surface plane constant
    Real nx{abs(NormalDirection[0])};
    Real ny{abs(NormalDirection[1])};
    const Real normalization{abs(nx)+abs(ny)};
    Real m{0};
    Real vof{0};
    // bring to the standard case 0<m1<m2 in Reference
    // doi:10.1006/jcph.2000.6567
    if (nx < ny) {
        m = nx / normalization;
    } else {
        m = ny / normalization;
    }
    if( VOF >= 1 || VOF <= 0) {
        vof = 0;
    } else if (VOF > 0.5){
        vof = 1 - VOF;
    } else {
        vof = VOF;
    }
    const Real Criticalvof{m / (2 * (1 - m))};
    // if (m >= 0 + 1e-14) {
        if (vof <= Criticalvof) {
            if (vof <= 0 + 1e-12) {
                    PlaneConst = 0;
            } else {
                    PlaneConst = sqrt(2 * m * (1 - m) * vof);
            }

        } else if (vof > Criticalvof && vof <= 0.5) {
            PlaneConst = vof * (1 - m) + 0.5 * m;
        }
        // if (VOF > 0.5) {
        //     PlaneConst = 1 - PlaneConst;
        // }
    // } else if (m < 0 + 1e-14) {
    //     PlaneConst = vof;
    // }  // for nx=0 or ny=0
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
    // printf("Error,  cx =%d cy = %d"  
    // "Error,  VOF = %e PlaneConst =%e nx = %e ny =%e \n", cx, cy, VOF, PlaneConst, Nx, Ny);
    //  *22*22*22*
    //  3        1
    //  3        1 
    //  3        1
    //  *00*00*00*
    int opp [4] {2,3,0,1};
    Real VofAtface[4]{};
    if (nx <= 0 + 1e-10 && ny <= 0+1e-10) {
       return VOF;
    }
// if (m >= 0 + 1e-14) {
    if (vof<Criticalvof){//vof=Criticalvof include nx or ny = 0 so must not
        VofAtface [0] = PlaneConst * normalization / nx;
        VofAtface [1] = 0;
        VofAtface [2] = 0;
        VofAtface [3] = PlaneConst * normalization / ny;
    }else if (vof >= Criticalvof && vof <= 0.5){
        if (nx >=ny){
            VofAtface [0] = PlaneConst * normalization / nx;
            VofAtface [1] = 0;
            VofAtface [2] = (PlaneConst * normalization - ny) / nx;
            VofAtface [3] = 1;
        }else{
            VofAtface [0] = 1;
            VofAtface [1] = (PlaneConst * normalization - nx) / ny;
            VofAtface [2] = 0;
            VofAtface [3] = (PlaneConst * normalization) / ny;
        }
    }
    if (VOF < 0.5) {
        for (int i = 0; i < 4; i++) {
            Root[i] = VofAtface[i];
        }
        if (Nx <= 0) {
            Root[1] = VofAtface[opp[1]];
            Root[3] = VofAtface[opp[3]];
        }
        if (Ny <= 0) {
            Root[0] = VofAtface[opp[0]];
            Root[2] = VofAtface[opp[2]];
        }
    } else {
        for (int i = 0; i < 4; i++) {
            Root[i] = 1 - VofAtface[opp[i]];
        }
        if (Nx <= 0) {
            Root[1] = 1 - VofAtface[1];
            Root[3] = 1 - VofAtface[3];
        }
        if (Ny <= 0) {
            Root[0] = 1 - VofAtface[0];
            Root[2] = 1 - VofAtface[2];
        }
    }
    for (int i = 0; i < 4; i++) {
        if ((isnan(Root[i]) || isinf(Root[i]))){
            ops_printf("Error,  nx =%e ny =%e\n", nx, ny);
            assert(!(isnan(Root[i]) || isinf(Root[i])));
        }
    }
    // printf("Error,  Nx =%e Ny =%e\n", Nx, Ny);
    if (cx == -1 && cy == 0){
        // printf("Error,  Root[3] =%e PlaneConst =%e\n", Root[3], PlaneConst);
        return Root[3];
    }else if (cy == -1 && cx == 0){
        // printf("Error,  Root[0] =%e PlaneConst =%e\n", Root[0], PlaneConst);
        return Root[0];
    }else if (cx == 1 && cy == 0){
        // printf("Error,  Root[1] =%e PlaneConst =%e\n", Root[1], PlaneConst);
        return Root[1];
    }else if (cy == 1 && cx == 0){
        // printf("Error,  Root[2] =%e PlaneConst =%e\n", Root[2], PlaneConst);
        return Root[2];
    }else if(cx == -1 && cy == -1) {
        if (cx*Nx+cy*Ny<=0){
            return (Root[3]+ Root[0])/2;
        }
    }else if(cx == -1 && cy == 1) {
        if (cx*Nx+cy*Ny<=0){
            return (Root[3]+ Root[2])/2;
        }
    }else if(cx == 1 && cy == -1) {
        if (cx*Nx+cy*Ny<=0){
            return (Root[1]+ Root[0])/2;
        }
    }else if(cx == 1 && cy == 1) {
        if (cx*Nx+cy*Ny<=0){
            return (Root[1]+ Root[2])/2;
        }
    }
    return VOF;
}
#endif  //  FREESURFACE_HOST_DEVICE_H
