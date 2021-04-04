
struct plastic_material
{
    lambertian_reflection Lambertian;
    microfacet_reflection Microfacets;
};

plastic_material PlasticMaterial()
{
    plastic_material Result = {};
    
    // Result.Microfacets = (microfacet_reflection*) malloc(sizeof(microfacet_reflection));
    // Result.Lambertian  = (lambertian_reflection*) malloc(sizeof(lambertian_reflection));
    
    Result.Microfacets = MicrofacetReflection(LaneV3(0.7,0.7,0.7), LaneF32FromF32(0.1f));
    Result.Lambertian = LambertianReflection(LaneV3(0.7, 0.7, 0.7));

    return Result;
}

void DeletePlasticMaterial(plastic_material *PlasticMaterial)
{
    // free(PlasticMaterial.Microfacets);
    // free(PlasticMaterial.Lambertian);
}

lane_f32 PlasticMaterial_Pdf(plastic_material *PlasticMaterial, lane_v3 *wo, lane_v3 * wi)  {
    lane_u32 ZIsZero = wo->z == LaneF32FromF32(0.0f);
    
    lane_f32 pdf = LaneF32FromF32(0.0f);
    lane_u32 matchingComps = LaneU32FromU32(2);
    
     
    pdf = pdf + LambertianReflection_Pdf(&PlasticMaterial->Lambertian, *wo, *wi);
    pdf = pdf + MicroFacetReflection_Pdf(&PlasticMaterial->Microfacets, *wo, *wi);
    
    lane_f32 v = pdf / LaneF32FromU32(matchingComps);

    ConditionalAssign(&v, ZIsZero, LaneF32FromF32(0.0f));
    
    return v;
}

lane_v3 PlasticMaterial_f(plastic_material *PlasticMaterial, lane_v3 *wo, lane_v3 * wi)  {
    lane_v3 Result = LaneV3(0,0,0);
     
    Result = Result + LambertianReflection_f(&PlasticMaterial->Lambertian, *wo, *wi);
    Result = Result + MicroFacetReflection_f(&PlasticMaterial->Microfacets, *wo, *wi);

    return Result;
}


lane_v3 PlasticMaterial_Sample_f(plastic_material *PlasticMaterial, lane_v3 &wo, lane_v3 * wi,  lane_v2 &u, lane_f32 BxDFSample, lane_f32 *pdf)  {
    //Check if there is at least one matching bxdf in the stored bxdfs
    lane_u32 matchingComps = LaneU32FromU32(2);
    lane_f32 matchingCompsFloat = LaneF32FromF32(2.0f);
    
    //If there is one matching, takes a random one in the list of all the matching bxdfs
    // lane_u32 comp = LaneU32FromF32(Min(Floor(u.x * matchingCompsFloat), matchingCompsFloat-LaneF32FromF32(1.0f)));
    lane_u32 comp = LaneU32FromF32((Lane_Min(Lane_Floor((BxDFSample) * matchingCompsFloat), matchingCompsFloat - LaneF32FromF32(1.0f))));
    
    *pdf=LaneF32FromF32(0.0f);

    lane_v3 LambertianWi, MicrofacetWi;
    lane_f32 LambertianPd = LaneF32FromF32(0.0f);
    lane_f32 MicrofacetPdf = LaneF32FromF32(0.0f);

    lane_v3 Result = {};

    lane_f32 pdfl, pdfm;
    lane_v3 wil, wim;
    lane_v3 brdfl, brdfm;
    brdfl = LambertianReflection_Sample_f(&PlasticMaterial->Lambertian, wo, &wil, u, &pdfl);
    brdfm = MicroFacetReflection_Sample_f(&PlasticMaterial->Microfacets, wo, &wim, u, &pdfm);

    lane_u32 IsLambertianMask = (comp == LaneU32FromU32(0));
    lane_u32 IsMicrofacetMask = (comp == LaneU32FromU32(1));
    
    *wi = wim;
    ConditionalAssign(wi, IsLambertianMask, wil);
    
    lane_u32 NullPdfmMask = (pdfm==LaneF32FromF32(0.0f)) & IsMicrofacetMask;
    lane_u32 NullPdflMask = (pdfl==LaneF32FromF32(0.0f)) & IsLambertianMask;
    lane_u32 NullPdfMask = NullPdfmMask |  NullPdflMask;
    
    //Calculating pdf for this direction 
    #if 1
    *pdf = pdfl + pdfm;
    *pdf = *pdf / matchingCompsFloat;
    
    //Calculating brdf for this direction 
    Result = brdfl + brdfm;  
    #else
    *wi = wil;
    *pdf =  pdfl;
    //Calculating brdf for this direction 
    Result =  brdfl;  
    #endif
    
    ConditionalAssign(&Result, NullPdfMask, LaneV3(0,0,0));

    return Result;        
}