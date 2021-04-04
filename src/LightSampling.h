lane_v2 UniformSampleTriangle(lane_v2 u) {
    lane_f32 su0 = Lane_SquareRoot(u.x);
    return LaneV2(1 - su0, u.y * su0);
}


hit SampleTriangle(triangle *Triangle, lane_v2 u) {
    lane_v2 b = UniformSampleTriangle(u);
    
    lane_v3 p0 = LaneV3FromV3(Triangle->V1);
    lane_v3 p1 = LaneV3FromV3(Triangle->V2);
    lane_v3 p2 = LaneV3FromV3(Triangle->V3);

    hit OutHit;
    OutHit.Position = b.x * p0 + b.y * p1 + (1-b.x-b.y)*p2;

    OutHit.Normal = Lane_NOZ(
        b.x * LaneV3FromV3(Triangle->N1) + 
        b.y * LaneV3FromV3(Triangle->N2) + 
        (1-b.x-b.y) * LaneV3FromV3(Triangle->N3)
    );

    return OutHit;
}

lane_f32 TriangleArea(triangle *Triangle)
{
    lane_f32 Result = LaneF32FromF32(0.0f);
    lane_v3 p0 = LaneV3FromV3(Triangle->V1);
    lane_v3 p1 = LaneV3FromV3(Triangle->V2);
    lane_v3 p2 = LaneV3FromV3(Triangle->V3);
    
    Result = 0.5f * Lane_Length(Lane_Cross(p1-p0, p2-p0));    

    return Result;
}


lane_f32 ShapePDF(shape *Triangle, hit *Hit, lane_v3 Wi){
    //TODO: Handle the case where the point is not visible from the shape
    //Ray ray = ref.SpawnRay(wi);
    lane_v3 RayOrigin = Hit->Position;
    lane_v3 RayDirection = Wi;
    f32 tHit;
    hit IntersectLight = {};
    lane_f32 Tolerance = LaneF32FromF32(0.0001f);
    lane_f32 MinHitDistance = LaneF32FromF32(0.001f);
        
    
    IntersectLight.Distance = LaneF32FromF32(F32Max); //Closest hit
    IntersectLight.MaterialIndex = LaneU32FromU32(0); //Index of the hit material
    IntersectLight.Position = {};
    IntersectLight.Normal = {};

    HitTriangle(Triangle, &IntersectLight, RayOrigin, RayDirection, Tolerance, MinHitDistance);
    
    lane_u32 ZeroPdfMask = (IntersectLight.MaterialIndex ==  LaneU32FromU32(0));
    
    
    //Probability that we hit the triangle at this point from that point
    
    lane_f32 Numerator = Lane_DistanceSq(Hit->Position, IntersectLight.Position);
    lane_f32 Area = TriangleArea(&Triangle->Fields.TriangleFields.Triangle);
    lane_f32 Denominator = (Lane_Abs(Lane_Inner(IntersectLight.Normal, -Wi)));
    lane_f32 pdf = Numerator / (Denominator * Area);
    
    ConditionalAssign(&pdf, ZeroPdfMask, LaneF32FromF32(0.0f));
    
    return pdf;
}

lane_v3 SampleEmitter(shape *EmitterShape, hit *Hit, lane_v2 &u, lane_v3 *wi, lane_f32 *pdf)
{
    hit SampledPosition = SampleTriangle(&EmitterShape->Fields.TriangleFields.Triangle, u);
    
    //Vector from the hit position to the light
    *wi=  Lane_NOZ(SampledPosition.Position - Hit->Position);
    
    //Probability to hit the shape from the point    
    *pdf = ShapePDF(EmitterShape, Hit, *wi);
    
    //TODO: Use the material emittance here instead of arbitrary 10
    lane_v3 EmittedLight = LaneV3(10, 10, 10);
    
    lane_u32 LightIsNullMask = Lane_Inner(SampledPosition.Normal, -*wi) < LaneF32FromF32(0.0f);
    ConditionalAssign(&EmittedLight, LightIsNullMask, LaneV3(0,0,0));

    return EmittedLight;
}

lane_u32 IsBlack(lane_v3 Value)
{
    lane_u32 Result = (Value.x == LaneF32FromF32(0.0f) & Value.y == LaneF32FromF32(0.0f) & Value.z == LaneF32FromF32(0.0f));
    return Result;
}

lane_f32 PowerHeuristic(u32 nf, lane_f32 fPdf, u32 ng, lane_f32 gPdf) {
    lane_f32 f = nf * fPdf;
    lane_f32 g = ng * gPdf;
    return (f*f) / (f*f + g*g);
}

lane_v3 SampleLights(hit *Hit, world *World, random_series* Series)
{
    if(World->LightsCount == 0) return LaneV3(0,0,0);

    lane_v3 WoLocal =  Lane_TransformDirection(Hit->InverseTransform, Hit->Wo);

    lane_v3 Ld = LaneV3(0.0f, 0.0f, 0.0f);

    lane_v3 wi;
    lane_f32 lightPdf = LaneF32FromF32(0.0f);
    lane_f32 scatteringPdf = LaneF32FromF32(0.0f);

    // VisibilityTester visibility;

    //Find one light
    lane_f32 RandomLaneF32 = Lane_Min((RandomUnilateral(Series)  * LaneF32FromU32(World->LightsCount)), LaneF32FromU32(World->LightsCount-1));
    f32 Random = Extract0(RandomLaneF32);
    shape *SampledTriangle = &World->Lights[(u32)Random];

    //Sample the light. Wi contains the direction from hitPoint to light, lightPdf contains the probability of sampling that light
    lane_v3 Li = SampleEmitter(SampledTriangle, Hit, LaneV2(RandomUnilateral(Series),RandomUnilateral(Series)), &wi, &lightPdf);
    
    lane_u32 LaneMask = lightPdf > 0;
#if 1
    LaneMask &= AndNot(IsBlack(Li), LaneU32FromU32(0xFFFFFFFF)); 

    if(!MaskIsZeroed(LaneMask)) {

        lane_v3 f = LaneV3(0,0,0);

        lane_v3 WiLocal =  Lane_TransformDirection(Hit->InverseTransform, wi);
        lane_v3 WoLocal =  Lane_TransformDirection(Hit->InverseTransform, Hit->Wo);
        
        
        lane_v3 MatRefColor  = GatherLaneV3(World->Materials, Hit->MaterialIndex, ReflectionColor);
        lane_f32 MatSpecular  = GatherF32(World->Materials, Hit->MaterialIndex, Specular);
        //compute the brdf value and its pdf for the generated direction
        
        lane_u32 MaterialType = GatherU32(World->Materials, Hit->MaterialIndex, Type); 

        lane_u32 LambertianMaterialMask = (MaterialType == LaneU32FromU32(material_types::Lambertian));
        if(!MaskIsZeroed(LambertianMaterialMask)) 
        {
            lambertian_reflection Lambertian = LambertianReflection(MatRefColor);
            f = LambertianReflection_f(&Lambertian, WoLocal, WiLocal);
            scatteringPdf = LambertianReflection_Pdf(&Lambertian, WoLocal, WiLocal);
        } 
        lane_u32 MicrofacetMaterialMask = (MaterialType == LaneU32FromU32(material_types::Microfacets));
        if(!MaskIsZeroed(MicrofacetMaterialMask)) 
        {
            microfacet_reflection Microfacets = MicrofacetReflection(MatRefColor);
            f = MicroFacetReflection_f(&Microfacets, WoLocal, WiLocal);
            scatteringPdf = MicroFacetReflection_Pdf(&Microfacets, WoLocal, WiLocal);
        } 
        lane_u32 PlasticMaterialMask = (MaterialType == LaneU32FromU32(material_types::Plastic));
        if(!MaskIsZeroed(PlasticMaterialMask)) 
        {
            plastic_material Plastic = PlasticMaterial();
            plastic_material Reflection = PlasticMaterial();
            Reflection.Lambertian.R = MatRefColor;
            Reflection.Microfacets.R = MatRefColor;
            
            f = PlasticMaterial_f(&Plastic, &WoLocal, &WiLocal);
            scatteringPdf = PlasticMaterial_Pdf(&Plastic, &WoLocal, &WiLocal);
        }
    
        //If the bxdf is not null
    
        LaneMask &= AndNot(IsBlack(f), LaneU32FromU32(0xFFFFFFFF));
        //TODO
        // if(!visibility.Unoccluded(scene)){ //If the path from point to light is occluded -> Li is 0
        //     Li = lane_v3(0);
        // }

        //At this point, if Li is black, the surface is occluded from light so no shading needed.
        LaneMask &= AndNot(IsBlack(Li), LaneU32FromU32(0xFFFFFFFF));
        
        lane_f32 weight = PowerHeuristic(1, lightPdf, 1, scatteringPdf);
        Ld +=  Lane_Hadamard(f,  Li) * weight * (1.0f / lightPdf);
    }

#endif

#if 1
    //SAMPLE BRDF

    lane_v3 f = LaneV3(0,0,0);
    

    lane_u32 MaterialType = GatherU32(World->Materials, Hit->MaterialIndex, Type); 
    //Generate a direction from the hit surface pdf

    lane_u32 LambertianMaterialMask = (MaterialType == LaneU32FromU32(material_types::Lambertian));
    if(!MaskIsZeroed(LambertianMaterialMask)) 
    {
        lambertian_reflection Lambertian = LambertianReflection(LaneV3(1,1,1));
        //Samples a scattering direction from the hit point    
        f = LambertianReflection_Sample_f(&Lambertian, WoLocal, &wi, LaneV2(RandomUnilateral(Series),RandomUnilateral(Series)), &scatteringPdf);
        wi = Lane_TransformDirection(Hit->Transform, wi);
        f = f * Lane_Abs(Lane_Inner(wi, Hit->Normal));
    } 
    lane_u32 MicrofacetMaterialMask = (MaterialType == LaneU32FromU32(material_types::Microfacets));
    if(!MaskIsZeroed(MicrofacetMaterialMask)) 
    {
        microfacet_reflection Microfacet = MicrofacetReflection(LaneV3(1,1,1));
        //Samples a scattering direction from the hit point    
        f = MicroFacetReflection_Sample_f(&Microfacet, WoLocal, &wi, LaneV2(RandomUnilateral(Series),RandomUnilateral(Series)), &scatteringPdf);
        wi = Lane_TransformDirection(Hit->Transform, wi);
        f = f * Lane_Abs(Lane_Inner(wi, Hit->Normal));
    } 
    lane_u32 PlasticMaterialMask = (MaterialType == LaneU32FromU32(material_types::Plastic));
    if(!MaskIsZeroed(PlasticMaterialMask)) 
    {
        plastic_material Plastic = PlasticMaterial();
        //Samples a scattering direction from the hit point    
        f = PlasticMaterial_Sample_f(&Plastic, WoLocal, &wi, LaneV2(RandomUnilateral(Series),RandomUnilateral(Series)), RandomUnilateral(Series), &scatteringPdf);
        wi = Lane_TransformDirection(Hit->Transform, wi);
        f = f * Lane_Abs(Lane_Inner(wi, Hit->Normal));
    }

    LaneMask &= AndNot(IsBlack(f), LaneU32FromU32(0xFFFFFFFF)); 
    LaneMask &= scatteringPdf > LaneF32FromF32(0.0f);

    if(!MaskIsZeroed(LaneMask)) {
        lane_f32 weight = LaneF32FromF32(1.0f);
        
        //Check if the generated direction hits the light
        lightPdf = ShapePDF(SampledTriangle, Hit, wi);
        
        lane_u32 HitLightMask = (lightPdf!=LaneF32FromF32(0.0f));
        
        weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
    
        
        
        
        //TODO
        // Ray ray = it.SpawnRay(wi);
        // lane_v3 Tr(1.0);
        // // bool foundSurfaceInteraction = handleMedia ? scene.IntersectTr(ray, sampler, &lightIsect, &Tr) : scene.Intersect(ray,  &lightIsect);
        bool foundSurfaceInteraction = true;
        

        lane_v3 Li = LaneV3(0.0f, 0.0f, 0.0f);
        // if(foundSurfaceInteraction) {
            Li = LaneV3(10, 10, 10);
        // }

        lane_v3 newLd = Ld + Lane_Hadamard(f, Li) * weight * (1.0f / scatteringPdf);

        lane_u32 LiIsNotBlackMask = AndNot(IsBlack(Li), LaneU32FromU32(0xFFFFFFFF));
        LiIsNotBlackMask &= HitLightMask;
        ConditionalAssign(&Ld, LiIsNotBlackMask , newLd);
    }      


#endif


    //free(Emitters);

    lane_v3 Result = LaneF32FromU32(World->LightsCount) * Ld;
    
    lane_u32 MissMask = AndNot(LaneMask, LaneU32FromU32(0xFFFFFFFF));
    ConditionalAssign(&Result, MissMask, LaneV3(0,0,0));
    return Result;
}
