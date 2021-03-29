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


#define COSINE_WEIGHTED 1

//RAYTRACING
lane_u32 SampleMaterial(material *Materials, hit Hit,random_series *Series, lane_v3 *Attenuation, lane_v3 *RayDirection) {
    //Get material properties
    lane_mat4 Transform = Hit.Transform;
    lane_mat4 InverseTransform = Hit.InverseTransform;
    lane_v3 LocalRayDirection = TransformDirection(InverseTransform, *RayDirection);

    lane_u32 MaterialType = GatherU32(Materials, Hit.MaterialIndex, Type); 

    lane_u32 DiffuseMaterialMask = (MaterialType == LaneU32FromU32(material_types::Diffuse));
    if(!MaskIsZeroed(DiffuseMaterialMask)) 
    {
        lane_v3 MatRefColor  = GatherV3(Materials, Hit.MaterialIndex, ReflectionColor);
        
        lane_v3 RandomBounce = NOZ(V3(0,0,1) + RandomInSphere(Series));
        lane_f32 CosineTerm = Abs(Inner(RandomBounce, V3(0,0,1)));
        *RayDirection = NOZ(TransformDirection(Transform, RandomBounce));

        ConditionalAssign(Attenuation, DiffuseMaterialMask, CosineTerm *  MatRefColor);
    }

    
    lane_u32 MetalicMaterialMask = (MaterialType == LaneU32FromU32(material_types::Metalic));
    if(!MaskIsZeroed(MetalicMaterialMask)) 
    {
        lane_v3 MatRefColor  = GatherV3(Materials, Hit.MaterialIndex, ReflectionColor);
        lane_f32 MatSpecular  = GatherF32(Materials, Hit.MaterialIndex, Specular);
        
#if 0
        lane_v3 ReflectedDirection = Reflect(LocalRayDirection, V3(0,0,1));
        ReflectedDirection =NOZ(ReflectedDirection + MatSpecular * RandomInSphere(Series));
        lane_f32 CosineTerm = Abs(Inner(ReflectedDirection, V3(0,0,1)));
        *RayDirection = TransformDirection(Transform, ReflectedDirection);
        ConditionalAssign(Attenuation, MetalicMaterialMask, CosineTerm *  MatRefColor);
#else

#define MICROFACET 0
#define LAMBERTIAN 0
#define PLASTIC 1   
#if MICROFACET
        //Sample BSDF to get the next path direction
        microfacet_reflection Reflection = MicrofacetReflection(V3(1,1,1));
        lane_v3 OutputDirection  = V3(0.0f, 0.0f, 0.0f);
        lane_f32 pdf=LaneF32FromF32(0.0f);
        lane_v2 Random = {RandomUnilateral(Series), RandomUnilateral(Series)};
        
        lane_v3 brdf = MicroFacetReflection_Sample_f(&Reflection, -LocalRayDirection, &OutputDirection,Random, &pdf);
#endif

#if LAMBERTIAN
        //Sample BSDF to get the next path direction
        lambertian_reflection Reflection = LambertianReflection(V3(1,1,1));
        lane_v3 OutputDirection  = V3(0.0f, 0.0f, 0.0f);
        lane_f32 pdf=LaneF32FromF32(0.0f);
        lane_v2 Random = {RandomUnilateral(Series), RandomUnilateral(Series)};

        lane_v3 brdf = LambertianReflection_Sample_f(&Reflection, -LocalRayDirection, &OutputDirection,Random, &pdf);
#endif

#if PLASTIC
        //Sample BSDF to get the next path direction
        plastic_material Reflection = PlasticMaterial();
        
        lane_v3 OutputDirection  = V3(0.0f, 0.0f, 0.0f);
        lane_f32 pdf=LaneF32FromF32(0.0f);
        lane_v2 Random = {RandomUnilateral(Series), RandomUnilateral(Series)};

        lane_v3 brdf = PlasticMaterial_Sample_f(&Reflection, -LocalRayDirection, &OutputDirection,Random, &pdf);

        // DeletePlasticMaterial(&Reflection);
#endif

        if(pdf>0)
        {
            lane_v3 MetalicAttenuation = brdf * Abs(Inner(OutputDirection, V3(0,0,1))) * (1.0f / pdf);
            *RayDirection = TransformDirection(Transform, OutputDirection);
            ConditionalAssign(Attenuation, MetalicMaterialMask, MetalicAttenuation);
        }

#endif
        // lane_f32 CosineTerm = LaneF32FromF32(1.0f);

    }
    
    
    lane_u32 DielectricMaterialMask = (MaterialType == LaneU32FromU32(material_types::Dielectric));
    if(!MaskIsZeroed(DielectricMaterialMask)) {
        lane_v3 MatRefColor  = GatherV3(Materials, Hit.MaterialIndex, ReflectionColor);
        
        lane_v3 ReflectedVector = Reflect(*RayDirection, Hit.Normal);

        lane_f32 IndexOfRefraction = GatherF32(Materials, Hit.MaterialIndex, IndexOfRefraction);

        lane_v3 OutNormal = Hit.Normal;
        lane_f32 NiOverNt = 1.0f / IndexOfRefraction;
        lane_f32 CosineTerm = -Inner(*RayDirection, Hit.Normal) / Length(*RayDirection); 
        
        lane_u32 IsInside = (Inner(*RayDirection, Hit.Normal) > 0);
        ConditionalAssign(&OutNormal, IsInside, -Hit.Normal);
        ConditionalAssign(&NiOverNt, IsInside, IndexOfRefraction);
        ConditionalAssign(&CosineTerm, IsInside,  IndexOfRefraction * Inner(*RayDirection, Hit.Normal) / Length(*RayDirection));

        lane_v3 Reflected = Reflect(*RayDirection, Hit.Normal);
        lane_v3 Refracted;
        lane_u32 DidRefract = Refract(*RayDirection, OutNormal, NiOverNt, &Refracted);
        
        lane_f32 ReflectionProbability = LaneF32FromF32(1.0f);
        ConditionalAssign(&ReflectionProbability, DidRefract, ShlickFresnelApproximation(CosineTerm, IndexOfRefraction));

        *RayDirection = NOZ(Refracted);
        lane_u32 ScatteredMask = RandomUnilateral(Series) < ReflectionProbability;
        ConditionalAssign(RayDirection, ScatteredMask, Reflected);
        
        ConditionalAssign(Attenuation, DielectricMaterialMask, MatRefColor);
    }

    return LaneU32FromU32(0xffffffff);
}

internal void GetCameraRay(camera Camera, lane_f32 FilmX, lane_f32 FilmY, lane_v3 SubPixelSample, random_series *Series, lane_v3 *OutRayOrigin, lane_v3 *OutRayDirection) {
    //Offset within the pixel
    lane_f32 OffX = (2.0f * SubPixelSample.x - 1.0f) * Camera.HalfPixW + FilmX + Camera.HalfPixW;
    lane_f32 OffY = (2.0f * SubPixelSample.y - 1.0f) * Camera.HalfPixH + FilmY + Camera.HalfPixH;
    
    lane_v3 FilmP = Camera.FilmCenter + OffX * V3(1,0,0) * Camera.HalfFilmW + OffY * V3(0,1,0) * Camera.HalfFilmH; //World position on the film to compute the ray
    
    lane_v3 PositionOnSensor = RandomInDisk(Series) * Camera.Aperture* LaneF32FromF32(0.5f);
    
    
    *OutRayOrigin = V3(1,0,0) * PositionOnSensor.x + V3(0,1,0) * PositionOnSensor.y ;
    *OutRayOrigin = TransformPosition(Camera.transform, *OutRayOrigin);
    
    *OutRayDirection = NOZ(FilmP);
    *OutRayDirection = TransformDirection(Camera.transform, *OutRayDirection);
}

camera CreateCamera(u32 Width, u32 Height)
{
    f32 FilmDistance =  (1.0f);
    lane_v3 CameraPosition = V3(0, 3, 10);
    lane_v3 FilmCenter = V3(0, 0, FilmDistance);

    camera Result = {};
    Result.transform = LookAt(CameraPosition, V3(0,0,0), V3(0, 1, 0));
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
    lane_v3 ZeroVector = V3(0,0,0);
    
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


    lane_v3 *SubPixelSamples = (lane_v3 *) malloc(LaneRayCount * sizeof(lane_v3));
    MultiJitter((lane_v3*)SubPixelSamples, RaysPerPixel, &State->Series);

    
    for(u32 RayIndex=0; RayIndex < LaneRayCount; ++RayIndex) {

        lane_f32 Tolerance = LaneF32FromF32(0.0001f);
        lane_f32 MinHitDistance = LaneF32FromF32(0.001f);
        
        lane_v3 Attenuation = V3(1,1,1);
        lane_v3 Sample = {};
        
        lane_u32 LaneMask = LaneU32FromU32(0xFFFFFFFF);


        lane_v3 RayOrigin, RayDirection;
        GetCameraRay(State->Camera, FilmX, FilmY, *(SubPixelSamples + RayIndex), &State->Series, &RayOrigin, &RayDirection);
        
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

            //PLANES
            for(u32 PlaneIndex=0; PlaneIndex < World->PlaneCount; ++PlaneIndex) {
                plane Plane = World->Planes[PlaneIndex];
                HitPlane(&Plane, &Hit, RayOrigin, RayDirection, Tolerance, MinHitDistance);
            }

            //SPHERES
            HitBVH(RayOrigin, RayDirection, Tolerance, MinHitDistance, BVH, &Hit, 0);

            //VOLUMES
            // for(u32 VolumeIndex=0; VolumeIndex < World->VolumeCount; ++VolumeIndex) {
            //     volume Volume = World->Volumes[VolumeIndex];
            //     HitVolume(&Volume, &Hit, RayOrigin, RayDirection, Tolerance, MinHitDistance, &Series);
            // }



            //When we hit a light, we multiply its intensity with the attenuation
            lane_v3 MatEmitColor = LaneMask & GatherV3(World->Materials, Hit.MaterialIndex, EmitionColor); //Get the emition colors of the rays that have hit an emitter
            Sample += Hadamard(Attenuation, MatEmitColor); //Add to the current pixel sample the value of the emitted color that we have hit times the current attenuation
            
            //Sample 1 light from the point of intersection
            //Sample += Sample1Light(Hit.Point);    

            //When we hit the material 0, we have missed. We kill the ray at that point, it's not gonna continue.
            LaneMask = LaneMask & (Hit.MaterialIndex != LaneU32FromU32( 0));
            
            //Same if we directly hit an emitted, we can stop tracing.
            LaneMask = LaneMask & (MatEmitColor == ZeroVector);
            
            //If all the rays in the lane are done, we break out of the loop
            if(MaskIsZeroed(LaneMask)) {
                break;
            } else {
                //Compute lighting
                // lane_v3 HitPointAttenuation = ComputeHitPointAttenuation(World, Hit.MaterialIndex, Hit, &RayDirection, &Series);
                
                lane_v3 HitPointAttenuation; 
                LaneMask &= SampleMaterial(World->Materials, Hit, &Series, &HitPointAttenuation, &RayDirection);
                
                RayOrigin = Hit.Position; 
                Attenuation = Hadamard(Attenuation, HitPointAttenuation);
            }
        }
        FinalColor += RayContrib * Sample;
    }
    free(SubPixelSamples);
    
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

    return true;
}



internal aabb SurroundAABBs(aabb A, aabb B) {
    lane_v3 SmallBound = { Min(A.min.x, B.min.x), Min(A.min.y, B.min.y), Min(A.min.z, B.min.z) };
    lane_v3 BigBound = { Max(A.max.x, B.max.x), Max(A.max.y, B.max.y), Max(A.max.z, B.max.z) };
    return { SmallBound, BigBound };
}

internal int BoxXCompare(const void *a, const void *b) {
    sphere *lhs = (sphere*)a;
    sphere *rhs = (sphere*)b;
    
    lane_u32 Result = (lhs->AABB.min.x - rhs->AABB.min.x < 0.0f);
    if(!MaskIsZeroed(Result)){
        return -1;
    } else {
        return 1;
    }
}

internal int BoxYCompare(const void *a, const void *b) {
    sphere *lhs = (sphere*)a;
    sphere *rhs = (sphere*)b;
    
    lane_u32 Result = (lhs->AABB.min.y - rhs->AABB.min.y < 0.0f);
    if(!MaskIsZeroed(Result)){
        return -1;
    } else {
        return 1;
    }
}

internal int BoxZCompare(const void *a, const void *b) {
    sphere* lhs = (sphere*)a;
    sphere* rhs = (sphere*)b;
    

    lane_u32 Result = (lhs->AABB.min.z - rhs->AABB.min.z < 0.0f);
    if(!MaskIsZeroed(Result)){
        return -1;
    } else {
        return 1;
    }
}

internal bvh *BuildBVH(sphere *Objects, int NumObjects, u32 level) {
    bvh *Result = new bvh();

    // u32 axis = u32(3 * RandomUnilateralSlow());
    u32 axis = level % 3;

    if(axis==0){
        qsort(Objects, NumObjects, sizeof(sphere), BoxXCompare);
    }
    else if(axis==1){
        qsort(Objects, NumObjects, sizeof(sphere), BoxYCompare);
    }
    else {
        qsort(Objects, NumObjects, sizeof(sphere), BoxZCompare);
    }

    if(NumObjects==1) {
        Result->Left = new bvh();
        Result->Right = new bvh();
        Result->Left->Spheres = Result->Right->Spheres = Objects;
        Result->Left->AABB = Result->Left->Spheres->AABB;
        Result->Right->AABB = Result->Right->Spheres->AABB;
    }
    else if(NumObjects==2) {
        Result->Left = new bvh();
        Result->Right = new bvh();
        Result->Left->Spheres = Objects;
        Result->Right->Spheres = Objects + 1;
        Result->Left->AABB = Result->Left->Spheres->AABB;
        Result->Right->AABB = Result->Right->Spheres->AABB;
    } else {
        level ++;
        Result->Left = BuildBVH(Objects, NumObjects/2, level);
        Result->Right = BuildBVH(Objects + NumObjects/2, NumObjects - NumObjects/2, level);
    }
    Result->AABB = SurroundAABBs(Result->Left->AABB, Result->Right->AABB);


    return Result;
}

#define TEST_STUFF 0

int main(int argCount, char **args) {

#if TEST_STUFF

    // mat4 lookat = LookAt(V3(0,0,0), NOZ(V3(0.8, 0.1, 0.2)), V3(0, 1, 0));
    // mat4_ lookat_ = LookAt_(V3(0,0,0), NOZ(V3(0.8, 0.1, 0.2)), V3(0, 1, 0));
    
    
    mat4 basis = OrthoBasisFromNormal(NOZ(V3(0.4, 0.7, 0.2)));
    mat4_ basis_ = OrthoBasisFromNormal_(NOZ(V3(0.4, 0.7, 0.2)));







#else

    //Specular, emition, reflection
    material Materials[] = {
        EmissiveMaterial({0.5f, 0.6f, 0.8f}),
        DiffuseMaterial({0.5f, 0.5f, 0.5f}),
        DiffuseMaterial({0.7f, 0.1f, 0.1f}),
        EmissiveMaterial({3.0f, 2.0f, 2.0f}),
        MetallicMaterial({0.2f, 0.8f, 0.2f}, 0.25F),
        DielectricMaterial(1.2f),
        VolumetricMaterial(1.2f)
    };
   
    //Normal, distance, matIndex
    plane Planes[] = {
        {{0, 1, 0},{1.2f}, 4}
        // ,{{0, 0, -1},{2}, 2} 
    };
 
 
    sphere Spheres[] {
        Sphere({-1,  -0.7, 0}, 0.5, 4),
		Sphere({0,  -0.7,0}, 0.5, 1)
        ,Sphere({0.5,  -0.7,1.5}, 0.5, 3)
        //,Sphere({0,  1,3}, 0.5, 5)
    };
    
    volume Volumes[] {
        {6 ,{{-0.5,-0.5,-0.5},{0.5,0.5,0.5}}}
    };
    
    
    world World = {};
    World.MaterialCount = ArrayCount(Materials);
    World.Materials = Materials;

    World.PlaneCount = ArrayCount(Planes);
    World.Planes = Planes;

    World.VolumeCount = ArrayCount(Volumes);
    World.Volumes = Volumes;
    
    World.SphereCount = ArrayCount(Spheres);
    World.Spheres = Spheres;

    World.BVH =  BuildBVH((sphere*)&Spheres, ArrayCount(Spheres), 0);

    image_32 Image = AllocateImage(1080, 720);
    
    u32 CoreCount = GetCPUCoreCount();
    u32 TileWidth;
    u32 TileHeight = TileWidth = 64;

    u32 TileCountX = (Image.Width + TileWidth-1) / TileWidth;
    u32 TileCountY = (Image.Height + TileHeight-1) / TileHeight;
    u32 TotalTileCount = TileCountX * TileCountY;
    
    work_queue Queue = {};
    Queue.MaxBounceCount = 8;
    Queue.RaysPerPixel = NextPowerOfTwo(128);
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


#endif
    return 0;
}