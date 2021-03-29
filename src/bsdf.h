
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
    
    Result.Microfacets = MicrofacetReflection(V3(1,0,0), 0.01f);
    Result.Lambertian = LambertianReflection(V3(0.4, 0.2, 0.8));

    return Result;
}

void DeletePlasticMaterial(plastic_material *PlasticMaterial)
{
    // free(PlasticMaterial.Microfacets);
    // free(PlasticMaterial.Lambertian);
}

lane_v3 PlasticMaterial_Sample_f(plastic_material *PlasticMaterial, lane_v3 &wo, lane_v3 * wi,  lane_v2 &u, f32 *pdf)  {
    //Check if there is at least one matching bxdf in the stored bxdfs
    u32 matchingComps = 2;
    
    //If there is one matching, takes a random one in the list of all the matching bxdfs
    u32 comp = Min((int)std::floor(u.x * matchingComps), matchingComps-1);
    
    *pdf=0;
    lane_v3 f;
    if(comp==0) //Lambertian
    {
        lane_v3 f = LambertianReflection_Sample_f(&PlasticMaterial->Lambertian, wo, wi, u, pdf);
    }
    else //Microfacets
    {
        lane_v3 f = MicroFacetReflection_Sample_f(&PlasticMaterial->Microfacets, wo, wi, u, pdf);
    }   
    
    if(*pdf==0) return V3(0,0,0);
    
    
    *pdf += LambertianReflection_Pdf(&PlasticMaterial->Lambertian, wo, *wi);
    *pdf += MicroFacetReflection_Pdf(&PlasticMaterial->Microfacets, wo, *wi);
    *pdf /= matchingComps;
    
    f=V3(0,0,0);
    f += LambertianReflection_f(&PlasticMaterial->Lambertian, wo, *wi);            
    f += MicroFacetReflection_f(&PlasticMaterial->Microfacets, wo, *wi);            
    
    return f;        
}