
#ifndef FREESURFACE_HOST_DEVICE_H
#define FREESURFACE_HOST_DEVICE_H
#ifndef OPS_FUN_PREFIX
#define OPS_FUN_PREFIX
#endif
// static inline OPS_FUN_PREFIX Real CalcBGKFeq(const int l, const Real rho, const Real u, const Real v,
//                 const Real w, const Real T, const int polyOrder) {
//     Real cu{(CS * XI[l * LATTDIM] * u + CS * XI[l * LATTDIM + 1] * v +
//              CS * XI[l * LATTDIM + 2] * w)};
//     Real c2{(CS * XI[l * LATTDIM] * CS * XI[l * LATTDIM] +
//              CS * XI[l * LATTDIM + 1] * CS * XI[l * LATTDIM + 1] +
//              CS * XI[l * LATTDIM + 2] * CS * XI[l * LATTDIM + 2])};
//     Real cu2{cu * cu};
//     Real u2{u * u + v * v + w * w};
//     Real res = 1.0 + cu + 0.5 * (cu2 - u2 + (T - 1.0) * (c2 - LATTDIM));
//     if ((polyOrder) >= 3) {
//         res = res +
//               cu * (cu2 - 3.0 * u2 + 3.0 * (T - 1.0) * (c2 - LATTDIM - 2.0)) /
//                   6.0;
//     }
//     if ((polyOrder) >= 4) {
//         res =
//             res + (cu2 * cu2 - 6.0 * cu2 * u2 + 3.0 * u2 * u2) / 24.0 +
//             (T - 1.0) * ((c2 - (LATTDIM + 2)) * (cu2 - u2) - 2.0 * cu2) / 4.0 +
//             (T - 1.0) * (T - 1.0) *
//                 (c2 * c2 - 2.0 * (LATTDIM + 2) * c2 + LATTDIM * (LATTDIM + 2)) /
//                 8.0;
//     }
//     return WEIGHTS[l] * rho * res;
// }


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
        if (VG_IPJP_I == vg || VG_IPJM_I == vg||
            VG_IMJP_I == vg || VG_IMJM_I == vg){
            NormalDirection[0] = 0.5*(vofs[1][1]-vofs[1-cx][1])*(-cx);
            NormalDirection[1] = 0.5*(vofs[1][1]-vofs[1][1-cy])*(-cy);
            break;
        }else {
            NormalDirection[0] += -weight*cx*vof/8.0;
            NormalDirection[1] += -weight*cy*vof/8.0;
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
#endif  //  FREESURFACE_HOST_DEVICE_H