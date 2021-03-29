
#include <algorithm>

inline lane_f32 CosTheta( lane_v3 w) {return w.z;}
inline lane_f32 Cos2Theta( lane_v3 w) {return w.z * w.z;}
inline lane_f32 AbsCosTheta( lane_v3 w) {return Abs(w.z);}

inline lane_f32 Sin2Theta( lane_v3 w) {return std::max(0.0, 1.0 - Cos2Theta(w));}
inline lane_f32 SinTheta( lane_v3 w) {return std::sqrt(Sin2Theta(w));}

inline lane_f32 TanTheta( lane_v3 w) {return SinTheta(w) / CosTheta(w); }
inline lane_f32 Tan2Theta( lane_v3 w) {return Sin2Theta(w) / Cos2Theta(w); }

inline lane_f32 CosPhi( lane_v3 w) {
    lane_f32 sin = SinTheta(w);
    return (sin==0) ? 0 : Clamp(w.x/sin, LaneF32FromF32(-1.0f), LaneF32FromF32(1.0f));
}

inline lane_f32 SinPhi( lane_v3 w) {
    lane_f32 sin = SinTheta(w);
    return (sin==0) ? 0 : Clamp(w.y/sin, LaneF32FromF32(-1.0f), LaneF32FromF32(1.0f));
}

inline lane_f32 Cos2Phi( lane_v3 w) {
    return CosPhi(w) * CosPhi(w);
}

inline lane_f32 Sin2Phi( lane_v3 w) {
    return SinPhi(w) * SinPhi(w);
}

inline lane_f32 CosDPhi( lane_v3 wa,  lane_v3 wb) {
    return Clamp((wa.x * wb.x + wa.y * wb.y) /
                 std::sqrt((wa.x * wa.x + wa.y * wa.y) *
                           (wb.x * wb.x + wb.y * wb.y)), -1.0, 1.0);
}

lane_v3 SphericalDirection(f32 sinTheta, f32 cosTheta, f32 phi) {
    return V3(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}

inline lane_u32 SameHemisphere( lane_v3 &w,  lane_v3 &wp) {
    return w.z * wp.z > 0;
}    

//Microfacets

struct beckmann_distribution
{
    lane_f32 AlphaX;
    lane_f32 AlphaY;
    b32 SampleVisibleArea;
};

lane_f32 BeckmannDistribution_D(beckmann_distribution Distribution, lane_v3 wh) {
    // lane_f32 tan2Theta = Tan2Theta(wh);
    // if(std::isinf(tan2Theta)) return 0;
    // lane_f32 cos4theta = Cos2Theta(wh)*Cos2Theta(wh);
    // return std::exp(-tan2Theta * (Cos2Phi(wh)/(Distribution.AlphaX*Distribution.AlphaX)+Sin2Phi(wh)/(Distribution.AlphaY*Distribution.AlphaY)))    / 
    //         //-----------------------------------------------------------------------------
    //                             (Pi32 * Distribution.AlphaX * Distribution.AlphaY * cos4theta);

    f32 tan2Theta = Tan2Theta(wh);
    if (std::isinf(tan2Theta)) return 0.;
    f32 cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
    return std::exp(-tan2Theta * (Cos2Phi(wh) / (Distribution.AlphaX * Distribution.AlphaX) +
                                  Sin2Phi(wh) / (Distribution.AlphaY * Distribution.AlphaY))) /
           (Pi32 * Distribution.AlphaX * Distribution.AlphaY * cos4Theta);
}

lane_f32 BeckmannDistribution_Lambda(beckmann_distribution Distribution, lane_v3 w) {
    // lane_f32 absTanTheta = std::abs(TanTheta(w));
    // if(std::isinf(absTanTheta)) return 0;
    
    // lane_f32 Alpha = std::sqrt(Cos2Phi(w) * Distribution.AlphaX * Distribution.AlphaX + Sin2Phi(w) * Distribution.AlphaY * Distribution.AlphaY);

    // lane_f32 a = 1 / (Alpha * absTanTheta);
    // if(a >= 1.6f) {
    //     return 0;
    // }

    // return (1 - 1.259f * a + 0.396 * a * a) / 
    //       //-------------------------------
    //         (3.535f * a + 2.181f * a * a);

    f32 absTanTheta = std::abs(TanTheta(w));
    if (std::isinf(absTanTheta)) return 0.;
    // Compute _alpha_ for direction _w_
    f32 alpha =
        std::sqrt(Cos2Phi(w) * Distribution.AlphaX * Distribution.AlphaX + Sin2Phi(w) * Distribution.AlphaY * Distribution.AlphaY);
    f32 a = 1 / (alpha * absTanTheta);
    if (a >= 1.6f) return 0;
    return (1 - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
}

    
f32 BeckmannDistribution_G1(beckmann_distribution Distribution, lane_v3 w) {
    return 1 / (1 + BeckmannDistribution_Lambda(Distribution, w));
}

f32 BeckmannDistribution_G(beckmann_distribution Distribution, lane_v3 wo, lane_v3 wi) {
    return 1 / (1 + BeckmannDistribution_Lambda(Distribution, wo) + BeckmannDistribution_Lambda(Distribution, wi));
}

lane_v3 BeckmannDistribution_Sample_wh(beckmann_distribution Distribution, lane_v3 wo, lane_v2 u) {
    // if(!Distribution.SampleVisibleArea) {
    //     lane_f32 tan2Theta, phi;
    //     if(Distribution.AlphaX == Distribution.AlphaY) {
    //         lane_f32 logSample = Log(LaneF32FromF32(1.0f) - u.x);
    //         if(std::isinf(logSample)) logSample=0;
    //         tan2Theta = -Distribution.AlphaX * Distribution.AlphaX * logSample;
    //         phi = u.y * LaneF32FromF32(2.0f * Pi32);
    //     } else {
    //         lane_f32 logSample = Log(u.x);
    //         phi = std::atan(Distribution.AlphaY / Distribution.AlphaX *
    //                         std::tan(2 * Pi32 * u.y + 0.5f * Pi32));
    //         if (u.y > 0.5f)
    //             phi += Pi32;
    //         lane_f32 sinPhi = std::sin(phi), cosPhi = std::cos(phi);
    //         lane_f32 AlphaX2 = Distribution.AlphaX * Distribution.AlphaX;
    //         lane_f32 AlphaY2 = Distribution.AlphaY * Distribution.AlphaY;
    //         tan2Theta = -logSample /
    //             (cosPhi * cosPhi / AlphaX2 + sinPhi * sinPhi / AlphaY2);
    //     }

    //     lane_f32 cosTheta = 1 / std::sqrt(1 + tan2Theta);
    //     lane_f32 sinTheta = std::sqrt(std::max((lane_f32)0, 1-cosTheta*cosTheta));
    //     lane_v3 wh = SphericalDirection(sinTheta, cosTheta, phi);
    //     if(!SameHemisphere(wo, wh)) wh = -wh;

    //     return wh;
    // } else {

    // }

    if (!Distribution.SampleVisibleArea) {
        // Sample full distribution of normals for Beckmann distribution

        // Compute $\tan^2 \theta$ and $\phi$ for Beckmann distribution sample
        f32 tan2Theta, phi;
        if (Distribution.AlphaX == Distribution.AlphaY) {
            f32 logSample = std::log(1 - u.x);
            tan2Theta = -Distribution.AlphaX * Distribution.AlphaX * logSample;
            phi = u.y * 2 * Pi32;
        } else {
            // Compute _tan2Theta_ and _phi_ for anisotropic Beckmann
            // distribution
            f32 logSample = std::log(1 - u.x);
            phi = std::atan(Distribution.AlphaY / Distribution.AlphaX *
                            std::tan(2 * Pi32 * u.y + 0.5f * Pi32));
            if (u.y > 0.5f) phi += Pi32;
            f32 sinPhi = std::sin(phi), cosPhi = std::cos(phi);
            f32 alphax2 = Distribution.AlphaX * Distribution.AlphaX;
            f32 alphay2 = Distribution.AlphaY * Distribution.AlphaY;
            tan2Theta = -logSample /
                        (cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
        }

        // Map sampled Beckmann angles to normal direction _wh_
        f32 cosTheta = 1 / std::sqrt(1 + tan2Theta);
        f32 sinTheta = std::sqrt(std::max((f32)0, 1 - cosTheta * cosTheta));
        lane_v3 wh = SphericalDirection(sinTheta, cosTheta, phi);
        if (!SameHemisphere(wo, wh)) wh = -wh;
        return wh;
    } else {
        
    }
}

lane_f32 BeckmannDistribution_Pdf(beckmann_distribution Distribution, lane_v3 wo, lane_v3 wh) {
    if(!Distribution.SampleVisibleArea) {
        return BeckmannDistribution_D(Distribution, wh) * BeckmannDistribution_G1(Distribution, wo) * Abs(Inner(wo, wh)) / AbsCosTheta(wo);
    }else {

    }
}

//////
struct fresnel_dielectric {
    lane_f32 etaI, etaT, k;
};

lane_f32 EvaluateFresnelDielectric(fresnel_dielectric *Fresnel, lane_f32 cosThetaI) 
{
    cosThetaI = Clamp(cosThetaI, -1.0, 1.0);

    bool entering = cosThetaI > 0;
    if(!entering) {
        lane_f32 tmp = Fresnel->etaI;
        Fresnel->etaI = Fresnel->etaT;
        Fresnel->etaT = tmp;
        cosThetaI = Abs(cosThetaI);
    }

    f32 sinThetaI = SquareRoot(Max(0.0f, 1.0f - cosThetaI*cosThetaI));
    f32 sinThetaT = Fresnel->etaI / Fresnel->etaT * sinThetaI;

    if(sinThetaT >= 1) return 1;

    f32 cosThetaT = SquareRoot(std::max(0.0, 1.0 - sinThetaT * sinThetaT));


    f32 rParl = ((Fresnel->etaT * cosThetaI) - (Fresnel->etaI * cosThetaT)) / 
                  ((Fresnel->etaT * cosThetaI) + (Fresnel->etaI * cosThetaT));

    f32 rPerp = ((Fresnel->etaI * cosThetaI) - (Fresnel->etaT * cosThetaT)) / 
                  ((Fresnel->etaI * cosThetaI) + (Fresnel->etaT * cosThetaT));

    return 0.5 * (rParl * rParl + rPerp * rPerp);    
}


//////
struct microfacet_reflection
{
    beckmann_distribution Distribution;
    
    lane_v3 R;
    fresnel_dielectric FresnelDielectric;
};

microfacet_reflection MicrofacetReflection(lane_v3 DiffuseColor, lane_f32 Roughness=0.1f, lane_f32 etaI=1.0f, lane_f32 etaT = 1.5f)
{
    microfacet_reflection Result = {};
    Result.R = DiffuseColor;

    Result.Distribution = {};
    Result.Distribution.AlphaX = Roughness;
    Result.Distribution.AlphaY = Roughness;
    Result.Distribution.SampleVisibleArea = false;

    Result.FresnelDielectric = {};
    Result.FresnelDielectric.etaI = etaI;
    Result.FresnelDielectric.etaT = etaT;
    Result.FresnelDielectric.k = 1.0f;

    return Result;
}

lane_v3 MicroFacetReflection_f(microfacet_reflection *MicroFacetReflection,  lane_v3 &wo,  lane_v3 &wi)  {
    f32 cosThetaO = AbsCosTheta(wo);
    f32 cosThetaI = AbsCosTheta(wi);

    lane_v3 wh = wi+wo;

    if(cosThetaI ==0 || cosThetaO==0) return V3(0.0f,0.0f,0.0f);
    if(wh.x ==0 && wh.y==0 && wh.z==0) return V3(0.0f,0.0f,0.0f);

    wh = NOZ(wh);
    lane_f32 F = EvaluateFresnelDielectric(&MicroFacetReflection->FresnelDielectric, Inner(wi, wh));
	
	return  MicroFacetReflection->R
			* BeckmannDistribution_D(MicroFacetReflection->Distribution, wh) 
			* BeckmannDistribution_G(MicroFacetReflection->Distribution, wo, wi) 
			* F 
			* (1.0f /  (4 * cosThetaI * cosThetaO));
}


lane_v3 MicroFacetReflection_Sample_f(microfacet_reflection *MicroFacetReflection,  lane_v3 &wo, lane_v3 *wi,  lane_v2 sample, f32 *pdf)  {
    lane_v3 wh = BeckmannDistribution_Sample_wh(MicroFacetReflection->Distribution, wo, sample);
    *wi = Reflect(-wo, wh);
    // if(!SameHemisphere(wo, *wi)) return V3(0.0f, 0.0f, 0.0f);

    *pdf = BeckmannDistribution_Pdf( MicroFacetReflection->Distribution, wo, wh)/(4*Inner(wo, wh));    
    
    return MicroFacetReflection_f(MicroFacetReflection, wo, *wi);
}

f32 MicroFacetReflection_Pdf(microfacet_reflection *MicroFacetReflection,  lane_v3& wo,  lane_v3& wi)  {
    if(!SameHemisphere(wo, wi)) return 0;
    lane_v3 wh = NOZ(wo+wi);

    return  BeckmannDistribution_Pdf(MicroFacetReflection->Distribution, wo, wh) * (1.0f / (4.0f*Inner(wo, wh)));    
}
