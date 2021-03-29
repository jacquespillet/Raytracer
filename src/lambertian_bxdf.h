
lane_v2 ConcentricSampleDisk(lane_v2 &u) {
    lane_v2 uOffset = 2.0f * u - V2(1.0f, 1.0f);
    if(uOffset.x == 0 && uOffset.y==0) {
        return V2(0.0f,0.0f);
    }
    lane_f32 theta, r;
    if (Abs(uOffset.x) > Abs(uOffset.y)) {
    r = uOffset.x;
    theta = PiOver4 * (uOffset.y / uOffset.x);
    } else {
        r = uOffset.y;
        theta = PiOver2 - PiOver4 * (uOffset.x / uOffset.y);
    }
    return r * V2(Cosine(theta), Sine(theta));

}


lane_v3 CosineSampleHemisphere(lane_v2 &u) {
    lane_v2 d = ConcentricSampleDisk(u);
    lane_f32 z = SquareRoot(Max((lane_f32)0.0f, 1-d.x*d.x - d.y*d.y));

    return V3(d.x, d.y, z);
}
lane_f32 CosineHemispherePdf() {
    return 1.0/(2.0 * Pi32);
}


struct lambertian_reflection
{
    lane_v3 R;
};

lambertian_reflection LambertianReflection(lane_v3 Color = V3(1,1,1))
{
    lambertian_reflection Result = {};
    Result.R = Color;
    return Result;
}

lane_v3 LambertianReflection_f(lambertian_reflection *LambertianReflection, lane_v3 wo, lane_v3 wi)
{
    return LambertianReflection->R * InvPi32;
}

f32 LambertianReflection_Pdf(lambertian_reflection *LambertianReflection,lane_v3 wo, lane_v3 wi)
{
    return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi32 : 0;
}

lane_v3 LambertianReflection_Sample_f(lambertian_reflection *LambertianReflection, lane_v3 wo,lane_v3 *wi,  lane_v2 sample, lane_f32 *pdf)
{
    *wi = CosineSampleHemisphere(sample);
    if(wo.z <0) wi->z *= -1;

    *pdf = LambertianReflection_Pdf(LambertianReflection, wo, *wi);
    return LambertianReflection_f(LambertianReflection, wo, *wi);
}