
lane_v2 ConcentricSampleDisk(lane_v2 &u) {
    lane_v2 Result = {};
    
    lane_v2 uOffset = LaneF32FromF32(2.0f) * u - LaneV2(LaneF32FromF32(1.0f), LaneF32FromF32(1.0f));
    lane_u32 OffsetIs0 = (uOffset.x == LaneF32FromF32(0.0f) & uOffset.y==LaneF32FromF32(0.0f));
    
    lane_u32 OffsetXGtOffsetY = Abs(uOffset.x) > Abs(uOffset.y);

    lane_f32 theta = LaneF32FromF32(PiOver2) - LaneF32FromF32(PiOver4) * (uOffset.x / uOffset.y);
    lane_f32 r =  r = uOffset.y;

    ConditionalAssign(&r, OffsetXGtOffsetY, uOffset.x);
    ConditionalAssign(&theta, OffsetXGtOffsetY, PiOver4 * (uOffset.y / uOffset.x));
    
    Result = r * LaneV2(Cos(theta), Sine(theta));
    ConditionalAssign(&Result, OffsetIs0, LaneV2(LaneF32FromF32(0.0f), LaneF32FromF32(0.0f)));
    
    return Result;
}


lane_v3 CosineSampleHemisphere(lane_v2 &u) {
    lane_v2 d = ConcentricSampleDisk(u);
    lane_f32 z = SquareRoot(Max(LaneF32FromF32(0.0f), LaneF32FromF32(1.0f) -d.x*d.x - d.y*d.y));

    return LaneV3(d.x, d.y, z);
}
lane_f32 CosineHemispherePdf() {
    return LaneF32FromF32(1.0f)/(LaneF32FromF32(2.0 * Pi32));
}


struct lambertian_reflection
{
    lane_v3 R;
};

lambertian_reflection LambertianReflection(lane_v3 Color = LaneV3(1,1,1))
{
    lambertian_reflection Result = {};
    Result.R = Color;
    return Result;
}

lane_v3 LambertianReflection_f(lambertian_reflection *LambertianReflection, lane_v3 wo, lane_v3 wi)
{
    return LambertianReflection->R * LaneF32FromF32(InvPi32);
}

lane_f32 LambertianReflection_Pdf(lambertian_reflection *LambertianReflection,lane_v3 wo, lane_v3 wi)
{
    lane_f32 Result = AbsCosTheta(wi) * LaneF32FromF32(InvPi32);
    lane_u32 NotSameHemisphere = AndNot(SameHemisphere(wo, wi), LaneU32FromU32(0xFFFFFFFF));
    ConditionalAssign(&Result, NotSameHemisphere, LaneF32FromF32(0.0f));
    return Result; 
    
}

lane_v3 LambertianReflection_Sample_f(lambertian_reflection *LambertianReflection, lane_v3 wo,lane_v3 *wi,  lane_v2 sample, lane_f32 *pdf)
{
    *wi = CosineSampleHemisphere(sample);
    lane_u32 zNegative = wo.z <0;
    lane_v3 negZ = *wi; negZ.z *= -1.0f;
    ConditionalAssign(wi, zNegative, negZ);
    

    *pdf = LambertianReflection_Pdf(LambertianReflection, wo, *wi);
    return LambertianReflection_f(LambertianReflection, wo, *wi);
}