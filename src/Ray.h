#include <stdio.h>
#include <stdint.h>
#define internal static
#define global static


///////////////////
//TYPES
///////////////////
typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int16_t  s16;
typedef int32_t  s32;
typedef int64_t  s64;

typedef float f32;
typedef double f64;

typedef bool b32;

#define ArrayCount(Array) (sizeof(Array) / sizeof((Array)[0]))
#define Assert(Expression) if(!(Expression)) {*(volatile int *)0 = 0;}

#define U32Max ((u32)-1)

#include "Maths.h"

#include "ray_lane.h"

struct random_series 
{
    lane_u32 State;
};

#include "Image.h"
#include "Random.h"

///////////////////
//3D OBJECTS
///////////////////

enum material_types {
    Diffuse,
    Metalic,
    Dielectric,
    Emissive,
    Volumetric,
    NumberOfMaterialTypes
};

struct hit;

struct material 
{
    f32 Specular;
    f32 Density;
    f32 IndexOfRefraction;
    v3 EmitionColor;
    v3 ReflectionColor;
    material_types Type;
};

material MetallicMaterial(v3 ReflectionColor, f32 Specular) {
    material Result = {};
    Result.Specular = Specular;
    Result.EmitionColor={};
    Result.ReflectionColor = ReflectionColor;
    Result.Type = material_types::Metalic;
    return Result;
}

material DiffuseMaterial(v3 ReflectionColor) {
    material Result = {};
    Result.Specular = 0;
    Result.EmitionColor={};
    Result.ReflectionColor = ReflectionColor;
    Result.Type = material_types::Diffuse;
    return Result;    
}

material EmissiveMaterial(v3 EmitionColor) {
    material Result = {};
    Result.Specular = 0;
    Result.EmitionColor=EmitionColor;
    Result.ReflectionColor = {};
    Result.Type = material_types::Diffuse;
    return Result;    
}

material DielectricMaterial(f32 IndexOfRefraction) {
    material Result = {};
    Result.Specular = 0;
    Result.IndexOfRefraction = IndexOfRefraction;
    Result.EmitionColor={};
    Result.ReflectionColor = {1,1,1};
    Result.Type = material_types::Dielectric;
    return Result;    
}

material VolumetricMaterial(f32 Density) {
    material Result = {};
    Result.Specular = 0;
    Result.IndexOfRefraction = 0;
    Result.Density = Density;
    Result.EmitionColor={};
    Result.ReflectionColor = {1,1,1};
    Result.Type = material_types::Volumetric;
    return Result;    
}



struct aabb{
    lane_v3 min;
    lane_v3 max;
};


struct plane
{
public:
    v3 N;
    f32 d;
    u32 MatIndex;
    
    aabb AABB;
};

struct sphere
{
public:
    lane_mat4 Transform;
    
    f32 r;
    u32 MatIndex;
    aabb AABB;
};

struct triangle
{
public:
    lane_mat4 Transform;

    v3 V1;
    v3 LaneV2;
    v3 LaneV3;

    v3 Normal;

    u32 MatIndex;

    aabb AABB;
};

struct volume
{
public:
    u32 MatIndex;
    aabb AABB;
};

internal sphere Sphere(v3 P, f32 r, u32 MatIndex) {
    sphere Result = {};
    
    Result.Transform = Translate(Identity(), LaneV3(P.x, P.y, P.z));
    
    
    Result.r = r;
    Result.MatIndex = MatIndex;

    Result.AABB = {
        LaneV3(-r + P.x, -r + P.y, -r + P.z),
        LaneV3(r + P.x, r + P.y, r + P.z)
    };

    return Result;
}

internal triangle Triangle(v3 V1, v3 LaneV2, v3 LaneV3, lane_mat4 Transform, u32 MatIndex) {
    triangle Result = {};
    
    Result.Transform = Transform;
    
    Result.V1 = V1;
    Result.LaneV2 = LaneV2;
    Result.LaneV3 = LaneV3;

    Result.MatIndex = MatIndex;

    // v3 WorldV1 = TransformPosition(Transform, V1);
    // v3 WorldLaneV2 = TransformPosition(Transform, LaneV2);
    // v3 WorldLaneV3 = TransformPosition(Transform, LaneV3);

    // v3 MinPosition = Min(WorldV1, Min(WorldLaneV2, WorldLaneV3));
    // v3 MaxPosition = Max(WorldV1, Max(WorldLaneV2, WorldLaneV3));

    // Result.AABB = {
    //     LaneLaneV3FromLaneV3(MinPosition), LaneLaneV3FromLaneV3(MaxPosition)
    // };

    return Result;    
}


struct bvh {
public:
    bvh *Left;
    bvh *Right;
    aabb AABB;
    sphere *Spheres;
};

struct world 
{
    u32 MaterialCount;
    material *Materials;

    u32 PlaneCount;
    plane *Planes;

    u32 SphereCount;
    sphere *Spheres;

    u32 VolumeCount;
    volume *Volumes;

    bvh *BVH;
};

#include "Hit.h"

//THREADING

struct work_order 
{
    world *World;
    image_32 Image;
    u32 MinX;
    u32 MinY;
    u32 OnePastXMax;
    u32 OnePastYMax;

    random_series Entropy;
};


struct work_queue
{
    u32 WorkOrderCount;
    work_order *WorkOrders;
    
    volatile u64 NextWorkOrderIndex;
    volatile u64 BouncesComputed;
    volatile u64 TileRetiredCount;

    u32 RaysPerPixel;
    u32 MaxBounceCount;
};

//Ray tracing

struct camera {
    lane_mat4 transform;
    
    //Film
    lane_v3  FilmCenter; 
    lane_f32 HalfPixW; 
    lane_f32 HalfPixH; 
    lane_f32 HalfFilmW; 
    lane_f32 HalfFilmH; 
    lane_f32 FilmW;
    lane_f32 FilmH;

    lane_f32 Aperture;

    lane_f32 LensRadius;
};

struct cast_state 
{
    //Static information
    world *World; 
    u32 RaysPerPixel; 
    u32 BounceCount; 
    
    camera Camera;


    //variable
    f32 FilmX; 
    f32 FilmY; 

    //Random
    random_series Series;
    
    u64 BouncesComputed;
 
    v3 FinalColor;
};

#include "microfacets_bxdf.h"
#include "lambertian_bxdf.h"
#include "bsdf.h"
