
#include <stdio.h>
#include <stdlib.h>

#include "Ray.h"
#include <time.h>
#include <assert.h>

#include <windows.h>

//Resources for next session :

//Books
//file:///C:/Users/jacqu/Google%20Drive/Projets%20persos/Books/Graphics/Offline/Ray%20Tracing_%20the%20Rest%20of%20Your%20Life.pdf
//file:///C:/Users/jacqu/Google%20Drive/Projets%20persos/Books/Graphics/Offline/Kevin%20Suffern%20-%20Ray%20Tracing%20from%20the%20Ground%20Up-A%20K%20Peters%20(2007).pdf

//Monte carlo and importance sampling
//https://computergraphics.stackexchange.com/questions/4979/what-is-importance-sampling
//https://www.gamedev.net/blogs/entry/2261086-importance-sampling/
//http://www.pbr-book.org/3ed-2018/Light_Transport_I_Surface_Reflection/Sampling_Reflection_Functions.html

//Random number generation
//http://www.pbr-book.org/3ed-2018/Sampling_and_Reconstruction.html
//https://blog.demofox.org/2020/11/25/multiple-importance-sampling-in-1d/

//PBR :
//https://www.reddit.com/r/gamedev/comments/4wjfbv/cooktorrance_brdf/
//https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
//https://graphicscompendium.com/gamedev/15-pbr
//http://www.codinglabs.net/article_physically_based_rendering_cook_torrance.aspx
//https://schuttejoe.github.io/post/disneybsdf/


internal b32 RenderTile(work_queue *Queue);

//THREADING
internal u64 LockedAddAndReturnPreviousValue(u64 volatile *Value, u64 Addend) {
    u64 Result = InterlockedExchangeAdd64((volatile LONG64*) Value, Addend);
    return Result;
}


internal DWORD WINAPI WorkerThread(void *lpParameter) {
    work_queue *Queue = (work_queue*)lpParameter;
    while(RenderTile(Queue)){};
    return(0);
} 

internal void CreateWorkThread(void *Parameter) {
    //Creates a thread and make it run the function WorkerThread, and pass the work_queue parameter
    DWORD ThreadID;
    HANDLE ThreadHandle = CreateThread(
        0,
        0,
        WorkerThread,
        Parameter,
        0,
        &ThreadID
    );
    CloseHandle(ThreadHandle);
}

internal u32 GetCPUCoreCount(void) 
{
    SYSTEM_INFO Info;
    GetSystemInfo(&Info);
    u32 Result = Info.dwNumberOfProcessors;

#if defined DEBUG
    Result = 1;
#endif
    return Result;
}


//RAYTRACING
lane_u32 SampleMaterial(material *Materials, hit Hit,random_series *Series, lane_v3 *Attenuation, lane_v3 *RayDirection, lane_v2 BrdfSample) {
    //Get material properties
    lane_mat4 Transform = Hit.Transform;
    lane_mat4 InverseTransform = Hit.InverseTransform;
    lane_v3 LocalRayDirection = Lane_TransformDirection(InverseTransform, *RayDirection);

    lane_u32 MaterialType = GatherU32(Materials, Hit.MaterialIndex, Type); 

    lane_u32 DiffuseMaterialMask = (MaterialType == LaneU32FromU32(material_types::Diffuse));
    if(!MaskIsZeroed(DiffuseMaterialMask)) 
    {
        lane_v3 MatRefColor  = GatherLaneV3(Materials, Hit.MaterialIndex, ReflectionColor);
        
        lane_v3 RandomBounce = Lane_NOZ(LaneV3(0,0,1) + RandomInSphere(Series));
        lane_f32 CosineTerm = Lane_Abs(Lane_Inner(RandomBounce, LaneV3(0,0,1)));
        *RayDirection = Lane_NOZ(Lane_TransformDirection(Transform, RandomBounce));

        ConditionalAssign(Attenuation, DiffuseMaterialMask, CosineTerm *  MatRefColor);
    }

    
    lane_u32 MetalicMaterialMask = (MaterialType == LaneU32FromU32(material_types::Metalic));
    if(!MaskIsZeroed(MetalicMaterialMask)) 
    {
        lane_v3 MatRefColor  = GatherLaneV3(Materials, Hit.MaterialIndex, ReflectionColor);
        lane_f32 MatSpecular  = GatherF32(Materials, Hit.MaterialIndex, Specular);
        lane_v3 ReflectedDirection = Lane_Reflect(LocalRayDirection, LaneV3(0,0,1));
        ReflectedDirection =Lane_NOZ(ReflectedDirection + MatSpecular * RandomInSphere(Series));
        lane_f32 CosineTerm = Lane_Abs(Lane_Inner(ReflectedDirection, LaneV3(0,0,1)));
        *RayDirection = Lane_TransformDirection(Transform, ReflectedDirection);
        ConditionalAssign(Attenuation, MetalicMaterialMask, CosineTerm *  MatRefColor);
    }
    
    lane_u32 PlasticMaterialMask = (MaterialType == LaneU32FromU32(material_types::Plastic));
    if(!MaskIsZeroed(PlasticMaterialMask)) 
    {
        lane_v3 MatRefColor  = GatherLaneV3(Materials, Hit.MaterialIndex, ReflectionColor);
        lane_f32 MatSpecular  = GatherF32(Materials, Hit.MaterialIndex, Specular);

        
        //Sample BSDF to get the next path direction
        plastic_material Reflection = PlasticMaterial();
        Reflection.Lambertian.R = MatRefColor;
        Reflection.Microfacets.R = MatRefColor;
        
        lane_v3 OutputDirection  = LaneV3(0.0f, 0.0f, 0.0f);
        lane_f32 pdf=LaneF32FromF32(0.0f);
        // lane_v2 Random = BrdfSample;
        lane_v2 Random = {RandomUnilateral(Series), RandomUnilateral(Series)};
        lane_f32 RandomBxDF = RandomUnilateral(Series);
        lane_v3 brdf = PlasticMaterial_Sample_f(&Reflection, -LocalRayDirection, &OutputDirection,Random,RandomBxDF,  &pdf);


        lane_u32 pdfGT0 = pdf > 0;
        lane_u32 localMask = PlasticMaterialMask & pdfGT0;

        lane_v3 PlasticAttenuation = brdf * Lane_Abs(Lane_Inner(OutputDirection, LaneV3(0,0,1))) * (1.0f / pdf);            
        ConditionalAssign(RayDirection, localMask, Lane_TransformDirection(Transform, OutputDirection));
        ConditionalAssign(Attenuation, localMask, PlasticAttenuation);
    }

    lane_u32 MicrofacetsMaterialMask = (MaterialType == LaneU32FromU32(material_types::Microfacets));
    if(!MaskIsZeroed(MicrofacetsMaterialMask)) 
    {
        lane_v3 MatRefColor  = GatherLaneV3(Materials, Hit.MaterialIndex, ReflectionColor);
        lane_f32 MatSpecular  = GatherF32(Materials, Hit.MaterialIndex, Specular);
        
        //Sample BSDF to get the next path direction
        microfacet_reflection Reflection = MicrofacetReflection(MatRefColor);
        lane_v3 OutputDirection  = LaneV3(0.0f, 0.0f, 0.0f);
        lane_f32 pdf=LaneF32FromF32(0.0f);
        // lane_v2 Random = BrdfSample;
        lane_v2 Random = {RandomUnilateral(Series), RandomUnilateral(Series)};
        
        lane_v3 brdf = MicroFacetReflection_Sample_f(&Reflection, -LocalRayDirection, &OutputDirection,Random, &pdf);

        
        lane_u32 pdfGT0 = pdf > 0;
        lane_u32 localMask = MicrofacetsMaterialMask & pdfGT0;

        lane_v3 MicrofacetsAttenuation = brdf * Lane_Abs(Lane_Inner(OutputDirection, LaneV3(0,0,1))) * (1.0f / pdf);            
        ConditionalAssign(RayDirection, localMask, Lane_TransformDirection(Transform, OutputDirection));
        ConditionalAssign(Attenuation, localMask, MicrofacetsAttenuation);
    }


    lane_u32 LambertianMaterialMask = (MaterialType == LaneU32FromU32(material_types::Lambertian));
    if(!MaskIsZeroed(LambertianMaterialMask)) 
    {
        lane_v3 MatRefColor  = GatherLaneV3(Materials, Hit.MaterialIndex, ReflectionColor);
        lane_f32 MatSpecular  = GatherF32(Materials, Hit.MaterialIndex, Specular);
        
        //Sample BSDF to get the next path direction
        lambertian_reflection Reflection = LambertianReflection(MatRefColor);
        lane_v3 OutputDirection  = LaneV3(0.0f, 0.0f, 0.0f);
        lane_f32 pdf=LaneF32FromF32(0.0f);
        // lane_v2 Random = BrdfSample;
        lane_v2 Random = {RandomUnilateral(Series), RandomUnilateral(Series)};

        lane_v3 brdf = LambertianReflection_Sample_f(&Reflection, -LocalRayDirection, &OutputDirection,Random, &pdf);
        
        lane_u32 pdfGT0 = pdf > 0;
        lane_u32 localMask = LambertianMaterialMask & pdfGT0;

        lane_v3 LambertianAttenuation = brdf * Lane_Abs(Lane_Inner(OutputDirection, LaneV3(0,0,1))) * (1.0f / pdf);            
        ConditionalAssign(RayDirection, localMask, Lane_TransformDirection(Transform, OutputDirection));
        ConditionalAssign(Attenuation, localMask, LambertianAttenuation);
    }
    
    
    lane_u32 DielectricMaterialMask = (MaterialType == LaneU32FromU32(material_types::Dielectric));
    if(!MaskIsZeroed(DielectricMaterialMask)) {
        lane_v3 MatRefColor  = GatherLaneV3(Materials, Hit.MaterialIndex, ReflectionColor);
        
        lane_v3 ReflectedVector = Lane_Reflect(*RayDirection, Hit.Normal);

        lane_f32 IndexOfRefraction = GatherF32(Materials, Hit.MaterialIndex, IndexOfRefraction);

        lane_v3 OutNormal = Hit.Normal;
        lane_f32 NiOverNt = 1.0f / IndexOfRefraction;
        lane_f32 CosineTerm = -Lane_Inner(*RayDirection, Hit.Normal) / Lane_Length(*RayDirection); 
        
        lane_u32 IsInside = (Lane_Inner(*RayDirection, Hit.Normal) > 0);
        ConditionalAssign(&OutNormal, IsInside, -Hit.Normal);
        ConditionalAssign(&NiOverNt, IsInside, IndexOfRefraction);
        ConditionalAssign(&CosineTerm, IsInside,  IndexOfRefraction * Lane_Inner(*RayDirection, Hit.Normal) / Lane_Length(*RayDirection));

        lane_v3 Reflected = Lane_Reflect(*RayDirection, Hit.Normal);
        lane_v3 Refracted;
        lane_u32 DidRefract = Lane_Refract(*RayDirection, OutNormal, NiOverNt, &Refracted);
        
        lane_f32 ReflectionProbability = LaneF32FromF32(1.0f);
        ConditionalAssign(&ReflectionProbability, DidRefract, Lane_ShlickFresnelApproximation(CosineTerm, IndexOfRefraction));

        *RayDirection = Lane_NOZ(Refracted);
        lane_u32 ScatteredMask = RandomUnilateral(Series) < ReflectionProbability;
        ConditionalAssign(RayDirection, ScatteredMask, Reflected);
        
        ConditionalAssign(Attenuation, DielectricMaterialMask, MatRefColor);
    }

    return LaneU32FromU32(0xffffffff);
}

internal void GetCameraRay(camera Camera, lane_f32 FilmX, lane_f32 FilmY, lane_v2 SubPixelSample, random_series *Series, lane_v3 *OutRayOrigin, lane_v3 *OutRayDirection) {
    //Offset within the pixel
    lane_f32 OffX = (2.0f * SubPixelSample.x - 1.0f) * Camera.HalfPixW + FilmX + Camera.HalfPixW;
    lane_f32 OffY = (2.0f * SubPixelSample.y - 1.0f) * Camera.HalfPixH + FilmY + Camera.HalfPixH;
    
    lane_v3 FilmP = Camera.FilmCenter + OffX * LaneV3(1,0,0) * Camera.HalfFilmW + OffY * LaneV3(0,1,0) * Camera.HalfFilmH; //World position on the film to compute the ray
    
    lane_v3 PositionOnSensor = RandomInDisk(Series) * Camera.Aperture* LaneF32FromF32(0.5f);
    
    
    *OutRayOrigin = LaneV3(1,0,0) * PositionOnSensor.x + LaneV3(0,1,0) * PositionOnSensor.y ;
    *OutRayOrigin = Lane_TransformPosition(Camera.transform, *OutRayOrigin);
    
    *OutRayDirection = Lane_NOZ(FilmP);
    *OutRayDirection = Lane_TransformDirection(Camera.transform, *OutRayDirection);
}

camera CreateCamera(u32 Width, u32 Height)
{
    f32 FilmDistance =  (1.0f);
    lane_v3 CameraPosition = LaneV3(5, 5, 5);
    lane_v3 FilmCenter = LaneV3(0, 0, FilmDistance);

    camera Result = {};
    Result.transform = Lane_LookAt(CameraPosition, LaneV3(0,0,0), LaneV3(0, 1, 0));
    Result.Aperture = 0.00f;

    //Film
    Result.FilmW = 1.0f;
    Result.FilmH = 1.0f;
    
    //Aspect ratio
    if(Width > Height) {
        Result.FilmH = Result.FilmW * ((f32)Height / (f32)Width);
    }
    if(Height > Width) {
        Result.FilmH = Result.FilmH * ((f32)Width / (f32)Height);
    }   


    Result.HalfFilmH = Result.FilmH * 0.5f;
    Result.HalfFilmW = Result.FilmW * 0.5f;
    Result.FilmCenter = FilmCenter;
    Result.HalfPixW = 0.5f / Width;
    Result.HalfPixH = 0.5f / Height;

    return Result;
}


internal void CastRays(cast_state *State)
{
    lane_v3 ZeroVector = LaneV3(0,0,0);
    
    //Unpack the state
    world *World = State->World;
    bvh *BVH = World->BVH;
    u32 RaysPerPixel = State->RaysPerPixel;
    u32 BounceCount = State->BounceCount;
    lane_f32 FilmX = LaneF32FromF32(State->FilmX);
    lane_f32 FilmY = LaneF32FromF32(State->FilmY);
    random_series Series = State->Series;
    
    lane_u32 BouncesComputed = LaneU32FromU32(0);
    lane_v3 FinalColor = {};

    u32 LaneWidth = LANE_WIDTH;
    u32 LaneRayCount = (RaysPerPixel / LaneWidth);
    
    lane_f32 RayContrib = LaneF32FromF32( 1.0f / (f32)(RaysPerPixel)); //We accumulate within the pixel and divide at the end


    
    MultiJitter((lane_v2*)State->SubPixelSamples, LaneRayCount, &Series);
    
    MultiJitter(State->BrdfSamples[0], LaneRayCount, &Series);
    // ShuffleArray(BrdfSamples[0], LaneRayCount, &Series);
    MultiJitter(State->BrdfSamples[1], LaneRayCount, &Series);
    ShuffleArray(State->BrdfSamples[1], LaneRayCount, &Series);

    for(u32 RayIndex=0; RayIndex < LaneRayCount; ++RayIndex) {

        lane_f32 Tolerance = LaneF32FromF32(0.0001f);
        lane_f32 MinHitDistance = LaneF32FromF32(0.001f);
        
        lane_v3 Attenuation = LaneV3(1,1,1);
        lane_v3 Sample = {};
        
        lane_u32 LaneMask = LaneU32FromU32(0xFFFFFFFF);


        lane_v3 RayOrigin, RayDirection;
        GetCameraRay(State->Camera, FilmX, FilmY, State->SubPixelSamples[RayIndex], &State->Series, &RayOrigin, &RayDirection);
        
        //Casting the ray through the scene
        for(u32 RayCount=0; RayCount < BounceCount; RayCount++) {

            //Initialize values
            hit Hit = {};
            Hit.Distance = LaneF32FromF32(F32Max); //Closest hit
            Hit.MaterialIndex=LaneU32FromU32(0); //Index of the hit material
            Hit.Position = {};
            Hit.Normal = {};
            lane_u32 LaneIncrement = LaneU32FromU32(1);

            BouncesComputed += (LaneIncrement & LaneMask);

            HitBVH(RayOrigin, RayDirection, Tolerance, MinHitDistance, BVH, &Hit, 0);

            //When we hit a light, we multiply its intensity with the attenuation
            lane_v3 MatEmitColor = LaneMask & GatherLaneV3(World->Materials, Hit.MaterialIndex, EmitionColor); //Get the emition colors of the rays that have hit an emitter
            Sample += Lane_Hadamard(Attenuation, MatEmitColor); //Add to the current pixel sample the value of the emitted color that we have hit times the current attenuation
            
            //Sample 1 light from the point of intersection
            //Sample += Sample1Light(Hit.Point);    

            //When we hit the material 0, we have missed. We kill the ray at that point, it's not gonna continue.
            LaneMask = LaneMask & (Hit.MaterialIndex != LaneU32FromU32(0));
            
            //Same if we directly hit an emitted, we can stop tracing.
            LaneMask = LaneMask & (MatEmitColor == ZeroVector);
            
            //If all the rays in the lane are done, we break out of the loop
            if(MaskIsZeroed(LaneMask)) {
                break;
            } else {
                //Compute lighting
                // lane_v3 HitPointAttenuation = ComputeHitPointAttenuation(World, Hit.MaterialIndex, Hit, &RayDirection, &Series);
                
                lane_v2 brdfSample;
                if(RayCount < 2) {
                    //brdfSample = State->BrdfSamples[RayCount][RayIndex];
                    brdfSample = {RandomUnilateral(&Series), RandomUnilateral(&Series)};
                }
                else
                {
                    brdfSample = {RandomUnilateral(&Series), RandomUnilateral(&Series)};
                }

                lane_v3 HitPointAttenuation; 
                LaneMask &= SampleMaterial(World->Materials, Hit, &Series, &HitPointAttenuation, &RayDirection,  brdfSample);
                
                RayOrigin = Hit.Position; 
                Attenuation = Lane_Hadamard(Attenuation, HitPointAttenuation);
            }
        }
        FinalColor += RayContrib * Sample;
    }
    
    State->FinalColor = HorizontalAdd(FinalColor);
    State->BouncesComputed = HorizontalAdd(BouncesComputed);
    State->Series= Series;
}

internal b32 RenderTile(work_queue *Queue) {
    //This function processes one work order. It's called by the threads.

    //Get the index of the work order to process
    u64 WorkOrderIndex =  LockedAddAndReturnPreviousValue(&Queue->NextWorkOrderIndex, 1);
    
    //If we reached the total work order count we stop
    if(WorkOrderIndex >= Queue->WorkOrderCount) {
        return false;
    }

    work_order *Order = Queue->WorkOrders + WorkOrderIndex;

    //Get information about the order : Image, tile coordinates...    
    image_32 Image = Order->Image;
    u32 XMin = Order->MinX;
    u32 YMin = Order->MinY;
    u32 OnePassXMax = Order->OnePastXMax;
    u32 OnePassYMax = Order->OnePastYMax; 
    
    //Structure holding information about the current ray casting state.
    cast_state State={};

    State.RaysPerPixel=Queue->RaysPerPixel;
    State.BounceCount = Queue->MaxBounceCount;
    State.World = Order->World;

    //Camera
    State.Camera = CreateCamera(Image.Width, Image.Height);

    State.Series = Order->Entropy;


    u32 LaneWidth = LANE_WIDTH;
    u32 LaneRayCount = (State.RaysPerPixel / LaneWidth);
    State.SubPixelSamples = (lane_v2 *) malloc(LaneRayCount * sizeof(lane_v2));
    State.BrdfSamples[0] = (lane_v2*) malloc(LaneRayCount * sizeof(lane_v2));
    State.BrdfSamples[1] = (lane_v2*) malloc(LaneRayCount * sizeof(lane_v2));

    //Loop through all the pixels of the tile
    u64 BouncesComputed=0;
    for(u32 y =YMin; y<OnePassYMax; y++) {
        u32 *Out = GetPixelPointer(Image, XMin, y);
        State.FilmY = 2.0f * ((f32)y / (f32)Image.Height) - 1.0f;
        for(u32 x =XMin; x<OnePassXMax; x++) {
            State.FilmX = 2.0f * ((f32)x / (f32)Image.Width) - 1.0f;
            State.FinalColor = {0,0,0};
            CastRays(&State);


            f32 R = 255.0f * LinearToSRGB(State.FinalColor.x);
            f32 G = 255.0f * LinearToSRGB(State.FinalColor.y);
            f32 B = 255.0f * LinearToSRGB(State.FinalColor.z);
            f32 A = 255.0f;
            u32 BMPValue = ((Roundf32ToUint32(B) << 0)  |
                            (Roundf32ToUint32(G) << 8)  |
                            (Roundf32ToUint32(R) << 16) |
                            (Roundf32ToUint32(A) << 24));
                            

            *Out++ = BMPValue;
        }

    }
    LockedAddAndReturnPreviousValue(&Queue->BouncesComputed, State.BouncesComputed);
    LockedAddAndReturnPreviousValue(&Queue->TileRetiredCount, 1);

    free(State.SubPixelSamples);
    free(State.BrdfSamples[0]);
    free(State.BrdfSamples[1]);

    return true;
}



internal aabb SurroundAABBs(aabb A, aabb B) {
    v3 SmallBound = { Min(A.min.x, B.min.x), Min(A.min.y, B.min.y), Min(A.min.z, B.min.z) };
    v3 BigBound = { Max(A.max.x, B.max.x), Max(A.max.y, B.max.y), Max(A.max.z, B.max.z) };
    return { SmallBound, BigBound };
}

internal int BoxXCompare(const void *a, const void *b) {
    shape *lhs = (shape*)a;
    shape *rhs = (shape*)b;
    
    if(lhs->AABB.min.x - rhs->AABB.min.x < 0.0f)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

internal int BoxYCompare(const void *a, const void *b) {
    shape *lhs = (shape*)a;
    shape *rhs = (shape*)b;
    
    if(lhs->AABB.min.y - rhs->AABB.min.y < 0.0f)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

internal int BoxZCompare(const void *a, const void *b) {
    shape* lhs = (shape*)a;
    shape* rhs = (shape*)b;
    
    if(lhs->AABB.min.z - rhs->AABB.min.z < 0.0f)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

void AddShapesToWorld(world* World, shape *Shapes, u32 ShapesCount)
{
    World->Shapes = (shape*) realloc(World->Shapes, (World->ShapesCount + ShapesCount) * sizeof(shape)); 
    memcpy(World->Shapes + World->ShapesCount, Shapes, ShapesCount * sizeof(shape));
    World->ShapesCount += ShapesCount;
    free(Shapes);
}


clock_t TimeElapsedBVH=0;
internal bvh *BuildBVH(shape *Objects, int NumObjects, u32 level) {
    bvh *Result = new bvh();

    // u32 axis = u32(3 * RandomUnilateralSlow());
    u32 axis = level % 3;


    clock_t StartClock = clock();
    if(axis==0){
        qsort(Objects, NumObjects, sizeof(shape), BoxXCompare);
    }
    else if(axis==1){
        qsort(Objects, NumObjects, sizeof(shape), BoxYCompare);
    }
    else {
        qsort(Objects, NumObjects, sizeof(shape), BoxZCompare);
    }
    
    clock_t EndClock = clock();
    TimeElapsedBVH += EndClock - StartClock;

    if(NumObjects==1) {
        Result->Left = new bvh();
        Result->Right = new bvh();
        Result->Left->Shapes = Result->Right->Shapes = Objects;
        Result->Left->AABB = Result->Left->Shapes->AABB;
        Result->Right->AABB = Result->Right->Shapes->AABB;
    }
    else if(NumObjects==2) {
        Result->Left = new bvh();
        Result->Right = new bvh();
        Result->Left->Shapes = Objects;
        Result->Right->Shapes = Objects + 1;
        Result->Left->AABB = Result->Left->Shapes->AABB;
        Result->Right->AABB = Result->Right->Shapes->AABB;
    } else {
        level ++;
            // std::nth_element(v.begin(), v.begin() + v.size()/2, v.end());
            // std::nth_element(v.begin(), v.begin()+1, v.end(), std::greater<int>());
            // std::nth_element(Objects, Objects + NumObjects/2,
            //         Objects + NumObjects,
            //         [axis](shape a,
            //             shape b) {
            //             if(axis==0)
            //             {
            //                 return AABBCentroid(a.AABB).x < AABBCentroid(b.AABB).x;
            //             }
            //             else if(axis==1)
            //             {
            //                 return AABBCentroid(a.AABB).y < AABBCentroid(b.AABB).y;
            //             }
            //             else
            //             {
            //                 return AABBCentroid(a.AABB).z < AABBCentroid(b.AABB).z;
            //             }
            //         });
        Result->Left = BuildBVH(Objects, NumObjects/2, level);
        Result->Right = BuildBVH(Objects + NumObjects/2, NumObjects - NumObjects/2, level);
    }
    Result->AABB = SurroundAABBs(Result->Left->AABB, Result->Right->AABB);


    return Result;
}

#include "mesh.h"

#define TEST_STUFF 0

int main(int argCount, char **args) {

#if TEST_STUFF
    int size = 512;
    int numSamples = NextPowerOfTwo(1024);

    random_series Entropy = {LaneU32FromU32(rand())}; 
    lane_v2 *Samples = (lane_v2*) malloc(numSamples * sizeof(lane_v2));
    // MultiJitter(Samples, numSamples, &Entropy);
    Jitter(Samples, numSamples, &Entropy);
    //ShuffleArray(Samples, numSamples, &Entropy);


    image_32 Image = AllocateImage(size, size);

    for(int y=0; y<size; y++)
    {
        for(int x=0; x<size; x++)
        {
               u32 *Pixel = GetPixelPointer(Image, x, y); 
                *Pixel = 0;
        }
    }

    

    for(int i=0; i<numSamples; i++)
    {
        lane_v2 Sample = Samples[i];
        u32 x = Max(0, Min(size, Sample.x * size));
        u32 y = Max(0, Min(size, Sample.y * size));
        u32 *Pixel = GetPixelPointer(Image, x, y); 
        *Pixel = 0xffffffff;
    }

    WriteImage(Image, "random.bmp");
    free(Samples);
#else

    //Specular, emition, reflection
    material Materials[] = {
        EmissiveMaterial({0.7,0.7,0.7}), //Reserved for the sky
        DiffuseMaterial({0.5f, 0.5f, 0.5f}),
        MetallicMaterial({0.2f, 0.8f, 0.2f}, 0.25F),
        DielectricMaterial(1.2f),
        LambertianMaterial({0.0f, 1.0f, 0.0f}),
        MicrofacetsMaterial({0.1f, 0.8f, 0.2f}),
        PlasticMaterial({1.0f, 0.0f, 0.0f}),
        EmissiveMaterial({10,10,10})
    };
   
 
    // u32 Cube1Count;
    // shape* Cube1 = MeshFromFile("models/cube/cube.obj", &Cube1Count, 2, Translate(Identity(), V3(0, 0, 0)));
 
    u32 Cube2Count;
    shape* Cube2 = MeshFromFile("models/cube/cube.obj", &Cube2Count, 5, Translate(Identity(), V3(1, 0, 0)));
    shape *Sph2 = Sphere(V3(4, 4, 0), 1, 7);

    shape *Sph = Sphere(V3(0, 1, 0), 1, 4);
    
    u32 GroundCount;
    shape* Ground = MeshFromFile("models/quad/quad.obj", &GroundCount, 6, Scale(Identity(), V3(5, 5, 5)));

    
    world World = {};
    World.MaterialCount = ArrayCount(Materials);
    World.Materials = Materials;

    // AddShapesToWorld(&World, Cube1, Cube1Count);
    // AddShapesToWorld(&World, Cube2, Cube2Count);
    AddShapesToWorld(&World, Ground, GroundCount);
    AddShapesToWorld(&World, Sph, 1);
    // AddShapesToWorld(&World, Sph2, 1);

    printf("Building BVH....\n");
    World.BVH =  BuildBVH(World.Shapes, World.ShapesCount, 0);
    printf("Finished BVH %d \n ", TimeElapsedBVH);

    image_32 Image = AllocateImage(1080, 720);
    
    u32 CoreCount = GetCPUCoreCount();
    u32 TileWidth;
    u32 TileHeight = TileWidth = 64;

    u32 TileCountX = (Image.Width + TileWidth-1) / TileWidth;
    u32 TileCountY = (Image.Height + TileHeight-1) / TileHeight;
    u32 TotalTileCount = TileCountX * TileCountY;
    
    work_queue Queue = {};
    Queue.MaxBounceCount = 8;
    Queue.RaysPerPixel = Max(LANE_WIDTH, NextPowerOfTwo(128));
    if(argCount==2) {
        Queue.RaysPerPixel = atoi(args[1]);
    }
    

    //Allocate some work orders for each tile.
    //The threads will process these orders
    Queue.WorkOrders = (work_order*) malloc(TotalTileCount * sizeof(work_order));
    work_order* WorkOrder = Queue.WorkOrders; //A pointer to the start of the queue
    
    printf("Confioguration %d cores with %dx%d tiles (%dk/tile), %d-wide lanes \n", CoreCount, TileWidth, TileHeight, TileWidth * TileHeight * 4 / 1024, LANE_WIDTH);
    printf("Quality %d Rays per pixel, %d max bounces \n", Queue.RaysPerPixel, Queue.MaxBounceCount);

    //Loop that builds the work orders
    for(u32 TileY = 0; TileY < TileCountY; ++TileY) 
    {
        u32 MinY = TileY * TileHeight;
        u32 OnePastMaxY = MinY + TileHeight;
        if(OnePastMaxY > Image.Height) OnePastMaxY = Image.Height;
        for(u32 TileX = 0; TileX < TileCountX; ++TileX) 
        {
            u32 MinX = TileX * TileWidth;
            u32 OnePastMaxX = MinX + TileWidth;
            if(OnePastMaxX > Image.Width) OnePastMaxX = Image.Width;
            
            //Create the actual order
            work_order *Order = Queue.WorkOrders + Queue.WorkOrderCount++; //Get the pointer to the current order
            assert(Queue.WorkOrderCount <= TotalTileCount);

            //Give information about the tile to render
            Order->World = &World;
            Order->Image = Image;
            Order->MinX = MinX;
            Order->MinY = MinY;
            Order->OnePastXMax = OnePastMaxX;
            Order->OnePastYMax = OnePastMaxY;

#if (LANE_WIDTH==8)
            random_series Entropy = {LaneU32FromU32(
                                                     rand()
                                                    ,rand()
                                                    ,rand()
                                                    ,rand()
                                                    ,rand()
                                                    ,rand()
                                                    ,rand()
                                                    ,rand()
                                                    )
                                    }; 
#elif (LANE_WIDTH==4)
            random_series Entropy = {LaneU32FromU32(
                                                    9823742 + TileX * 12584 + TileY * 23423
                                                    ,68564 + TileX * 69826 + TileY * 146841
                                                    ,36454 + TileX * 671789 + TileY * 1566
                                                    ,6178+ TileX * 167498 + TileY * 324984
                                                    )
                                    }; 
#else
      random_series Entropy = {LaneU32FromU32(9823742 + TileX * 12584 + TileY * 23423)};
#endif
            Order->Entropy= Entropy;
        }        
    }

    assert(Queue.WorkOrderCount == TotalTileCount);

    //Fence to sync all the threads
    LockedAddAndReturnPreviousValue(&Queue.NextWorkOrderIndex, 0);



    clock_t StartClock = clock();
    
    //For all the cores, start running
    for(u32 CoreIndex =1; CoreIndex < CoreCount; CoreIndex++) {
        CreateWorkThread(&Queue);
    }

    //Main thread render
    while(Queue.TileRetiredCount < TotalTileCount) {
        if(RenderTile(&Queue)) {
            printf("\rRaycasting tile %d....\n", 100 * (u32)Queue.TileRetiredCount / TotalTileCount);
            fflush(stdout);
        }   
    }

    
    clock_t EndClock = clock();
    clock_t TimeElapsed = EndClock - StartClock;
    printf("\nTime spent %dms\n", TimeElapsed);
    printf("\nTotal bounces %llu\n", Queue.BouncesComputed);
    printf("\nMS Per bounce%f\n", (f64)TimeElapsed / (f64)Queue.BouncesComputed);
        
    WriteImage(Image, "test.bmp");    

    free(World.Shapes);

#endif
    return 0;
}