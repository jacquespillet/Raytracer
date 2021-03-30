
#include <algorithm>

inline lane_f32 CosTheta( lane_v3 w) {return w.z;}
inline lane_f32 Cos2Theta( lane_v3 w) {return w.z * w.z;}
inline lane_f32 AbsCosTheta( lane_v3 w) {return Abs(w.z);}

inline lane_f32 Sin2Theta( lane_v3 w) {return Max(LaneF32FromF32(0.0), LaneF32FromF32(1.0) - Cos2Theta(w));}
inline lane_f32 SinTheta( lane_v3 w) {return SquareRoot(Sin2Theta(w));}

inline lane_f32 TanTheta( lane_v3 w) {return SinTheta(w) / CosTheta(w); }
inline lane_f32 Tan2Theta( lane_v3 w) {return Sin2Theta(w) / Cos2Theta(w); }

inline lane_f32 CosPhi( lane_v3 w) {
    lane_f32 sin = SinTheta(w);
    
    lane_u32 SinIs0 = (sin==LaneF32FromF32(0));

    lane_f32 Result = LaneF32FromF32(0);
    ConditionalAssign(&Result, SinIs0, Clamp(w.x/sin, LaneF32FromF32(-1.0f), LaneF32FromF32(1.0f)));

    return Result;
}

inline lane_f32 SinPhi( lane_v3 w) {
    lane_f32 sin = SinTheta(w);
    lane_u32 SinIs0 = (sin==LaneF32FromF32(0));

    lane_f32 Result = LaneF32FromF32(0);
    ConditionalAssign(&Result, SinIs0, Clamp(w.y/sin, LaneF32FromF32(-1.0f), LaneF32FromF32(1.0f)));
    
    return Result;
}

inline lane_f32 Cos2Phi( lane_v3 w) {
    return CosPhi(w) * CosPhi(w);
}

inline lane_f32 Sin2Phi( lane_v3 w) {
    return SinPhi(w) * SinPhi(w);
}

inline lane_f32 CosDPhi( lane_v3 wa,  lane_v3 wb) {
    return Clamp((wa.x * wb.x + wa.y * wb.y) /
                 SquareRoot((wa.x * wa.x + wa.y * wa.y) *
                           (wb.x * wb.x + wb.y * wb.y)), LaneF32FromF32(-1.0f), LaneF32FromF32(1.0f));
}

lane_v3 SphericalDirection(lane_f32 sinTheta, lane_f32 cosTheta, lane_f32 phi) {
    return LaneV3(sinTheta * Cos(phi), sinTheta * Sine(phi), cosTheta);
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
    lane_f32 tan2Theta = Tan2Theta(wh);
    
    lane_u32 laneMask = (tan2Theta < LaneF32FromF32(FLT_MAX) & tan2Theta > LaneF32FromF32(-FLT_MAX));
    // lane_u32 laneMask = (isfinite(tan2Theta));
    
    lane_f32 cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
    lane_f32 FinalValue = Exp(-tan2Theta * (Cos2Phi(wh) / (Distribution.AlphaX * Distribution.AlphaX) +
                                  Sin2Phi(wh) / (Distribution.AlphaY * Distribution.AlphaY))) /
           (Pi32 * Distribution.AlphaX * Distribution.AlphaY * cos4Theta);

    lane_f32 Result = LaneF32FromF32(0.0f);
    ConditionalAssign(&Result, laneMask, FinalValue);
    return Result;
}

lane_f32 BeckmannDistribution_Lambda(beckmann_distribution Distribution, lane_v3 w) {

    lane_f32 absTanTheta = Abs(TanTheta(w));
    
    lane_u32 laneMask = (absTanTheta < LaneF32FromF32(FLT_MAX) & absTanTheta > LaneF32FromF32(-FLT_MAX));
    // Compute _alpha_ for direction _w_
    lane_f32 alpha = SquareRoot(Cos2Phi(w) * Distribution.AlphaX * Distribution.AlphaX + Sin2Phi(w) * Distribution.AlphaY * Distribution.AlphaY);
    
    lane_f32 a = LaneF32FromF32(1.0f) / (alpha * absTanTheta);
    
    laneMask &= a < LaneF32FromF32(1.6f);
    
    lane_f32 FinalValue = (1 - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);

    lane_f32 Result = LaneF32FromF32(0.0f);
    ConditionalAssign(&Result, laneMask, FinalValue);
    return Result;
}

    
lane_f32 BeckmannDistribution_G1(beckmann_distribution Distribution, lane_v3 w) {
    return LaneF32FromF32(1.0f) / (LaneF32FromF32(1.0f) + BeckmannDistribution_Lambda(Distribution, w));
}

lane_f32 BeckmannDistribution_G(beckmann_distribution Distribution, lane_v3 wo, lane_v3 wi) {
    return LaneF32FromF32(1.0f) / (LaneF32FromF32(1.0f) + BeckmannDistribution_Lambda(Distribution, wo) + BeckmannDistribution_Lambda(Distribution, wi));
}
lane_v3 BeckmannDistribution_Sample_wh(beckmann_distribution Distribution, lane_v3 wo, lane_v2 u) {
    if (!Distribution.SampleVisibleArea) {
        // Sample full distribution of normals for Beckmann distribution

        // Compute $\tan^2 \theta$ and $\phi$ for Beckmann distribution sample
        lane_f32 tan2Theta, phi;
        
        //TODO(Jacques): Treat this case
        //if (Distribution.AlphaX == Distribution.AlphaY) {
        
            lane_f32 logSample = Log(LaneF32FromF32(1.0f) - u.x);
            tan2Theta = -Distribution.AlphaX * Distribution.AlphaX * logSample;
            phi = u.y * 2 * Pi32;
        
        // } else {
        //     // Compute _tan2Theta_ and _phi_ for anisotropic Beckmann
        //     // distribution
        //     lane_f32 logSample = std::log(1 - u.x);
        //     phi = std::atan(Distribution.AlphaY / Distribution.AlphaX *
        //                     std::tan(2 * Pi32 * u.y + 0.5f * Pi32));
        //     if (u.y > 0.5f) phi += Pi32;
        //     lane_f32 sinPhi = std::sin(phi), cosPhi = std::cos(phi);
        //     lane_f32 alphax2 = Distribution.AlphaX * Distribution.AlphaX;
        //     lane_f32 alphay2 = Distribution.AlphaY * Distribution.AlphaY;
        //     tan2Theta = -logSample /
        //                 (cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
        // }

        // Map sampled Beckmann angles to normal direction _wh_
        lane_f32 cosTheta = 1 / SquareRoot(1 + tan2Theta);
        lane_f32 sinTheta = SquareRoot(Max(LaneF32FromF32(0.0f), LaneF32FromF32(1.0f) - cosTheta * cosTheta));
        lane_v3 wh = SphericalDirection(sinTheta, cosTheta, phi);
        
        lane_u32 SameHemisphereMask = AndNot(SameHemisphere(wo, wh), LaneU32FromU32(0xffffffff));

        lane_v3 Result = wh;

        ConditionalAssign(&Result, SameHemisphereMask, -wh);
        return Result;
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
    cosThetaI = Clamp(cosThetaI, LaneF32FromF32(-1.0f), LaneF32FromF32(1.0f));

    lane_u32 entering = cosThetaI > LaneF32FromF32(0.0f);
    

    //TODO(Jacques): HANDLE THAT
    // if(!entering) {
    //     lane_f32 tmp = Fresnel->etaI;
    //     Fresnel->etaI = Fresnel->etaT;
    //     Fresnel->etaT = tmp;
    //     cosThetaI = Abs(cosThetaI);
    // }


    lane_f32 sinThetaI = SquareRoot(Max(LaneF32FromF32(0.0f), LaneF32FromF32(1.0f) - cosThetaI*cosThetaI));
    lane_f32 sinThetaT = Fresnel->etaI / Fresnel->etaT * sinThetaI;

    lane_u32 SinThetaGTOneMask = (sinThetaT >= 1);
    
    //if(sinThetaT >= 1) return 1;

    lane_f32 cosThetaT = SquareRoot(Max(LaneF32FromF32(0.0f), LaneF32FromF32(1.0f) - sinThetaT * sinThetaT));


    lane_f32 rParl = ((Fresnel->etaT * cosThetaI) - (Fresnel->etaI * cosThetaT)) / 
                  ((Fresnel->etaT * cosThetaI) + (Fresnel->etaI * cosThetaT));

    lane_f32 rPerp = ((Fresnel->etaI * cosThetaI) - (Fresnel->etaT * cosThetaT)) / 
                  ((Fresnel->etaI * cosThetaI) + (Fresnel->etaT * cosThetaT));

    lane_f32 Result = LaneF32FromF32(0.5f) * (rParl * rParl + rPerp * rPerp); 
    ConditionalAssign(&Result, SinThetaGTOneMask, LaneF32FromF32(1.0f));

    return Result;
}


//////
struct microfacet_reflection
{
    beckmann_distribution Distribution;
    
    lane_v3 R;
    fresnel_dielectric FresnelDielectric;
};

microfacet_reflection MicrofacetReflection(lane_v3 DiffuseColor, lane_f32 Roughness=LaneF32FromF32(0.1f), lane_f32 etaI=LaneF32FromF32(1.0f), lane_f32 etaT =LaneF32FromF32( 1.5f))
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
    lane_f32 cosThetaO = AbsCosTheta(wo);
    lane_f32 cosThetaI = AbsCosTheta(wi);

    lane_v3 wh = wi+wo;


    lane_u32 laneMask = (cosThetaI ==LaneF32FromF32(0) | cosThetaO==LaneF32FromF32(0));
    laneMask |= (wh.x ==LaneF32FromF32(0) & wh.y==LaneF32FromF32(0) & wh.z==LaneF32FromF32(0));
    

    wh = NOZ(wh);
    lane_f32 F = EvaluateFresnelDielectric(&MicroFacetReflection->FresnelDielectric, Inner(wi, wh));
	
    lane_v3 Result = MicroFacetReflection->R
			* BeckmannDistribution_D(MicroFacetReflection->Distribution, wh) 
			* BeckmannDistribution_G(MicroFacetReflection->Distribution, wo, wi) 
			* F 
			* (1.0f /  (4 * cosThetaI * cosThetaO));

    ConditionalAssign(&Result, laneMask, LaneV3(0.0f,0.0f,0.0f));
	return  Result;
}


lane_v3 MicroFacetReflection_Sample_f(microfacet_reflection *MicroFacetReflection,  lane_v3 &wo, lane_v3 *wi,  lane_v2 sample, lane_f32 *pdf)  {
    lane_v3 wh = BeckmannDistribution_Sample_wh(MicroFacetReflection->Distribution, wo, sample);
    *wi = Reflect(-wo, wh);
    
    
    *pdf = BeckmannDistribution_Pdf( MicroFacetReflection->Distribution, wo, wh)/(4*Inner(wo, wh));    
    
    lane_v3 Result =MicroFacetReflection_f(MicroFacetReflection, wo, *wi);
    
    //TODO(Jacques)
    lane_u32 SameHemisphereMask = AndNot(SameHemisphere(wo, *wi), LaneU32FromU32(0xFFFFFFFF));
    ConditionalAssign(&Result, SameHemisphereMask, LaneV3(0.0f, 0.0f, 0.0f));

    return Result;
}

lane_f32 MicroFacetReflection_Pdf(microfacet_reflection *MicroFacetReflection,  lane_v3& wo,  lane_v3& wi)  {

    lane_v3 wh = NOZ(wo+wi);
    lane_f32 Result =  BeckmannDistribution_Pdf(MicroFacetReflection->Distribution, wo, wh) * (LaneF32FromF32(1.0f) / (LaneF32FromF32(4.0f)*Inner(wo, wh))); 

    lane_u32 SameHemisphereMask = AndNot(SameHemisphere(wo, wi), LaneU32FromU32(0xFFFFFFFF));
    ConditionalAssign(&Result, SameHemisphereMask, LaneF32FromF32(0.0f));

    return Result;  
}
