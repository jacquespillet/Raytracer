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
    Lambertian,
    Microfacets,
    Plastic,
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
    Result.Type = material_types::Emissive;
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

material LambertianMaterial(v3 ReflectionColor) {
    material Result = {};
    Result.Specular = 0;
    Result.IndexOfRefraction = 1.0f;
    Result.EmitionColor={};
    Result.ReflectionColor = ReflectionColor;
    Result.Type = material_types::Lambertian;
    return Result;    
}

material MicrofacetsMaterial(v3 ReflectionColor) {
    material Result = {};
    Result.Specular = 0;
    Result.IndexOfRefraction = 1.0f;
    Result.EmitionColor={};
    Result.ReflectionColor = ReflectionColor;
    Result.Type = material_types::Microfacets;
    return Result;    
}

material PlasticMaterial(v3 ReflectionColor) {
    material Result = {};
    Result.Specular = 0;
    Result.IndexOfRefraction = 1.0f;
    Result.EmitionColor={};
    Result.ReflectionColor = ReflectionColor;
    Result.Type = material_types::Plastic;
    return Result;    
}


#define MAX_SHAPE_STRUCT_SIZE 2048

enum shape_type
{
    sphereType,
    triangleType,
    numShapeTypes
};

struct aabb{
    v3 min;
    v3 max;
};

v3 AABBCentroid(aabb AABB)
{
    return (AABB.min + AABB.max) * 0.5f;
}


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
    f32 r;
};

struct triangle
{
public:
    v3 V1;
    v3 V2;
    v3 V3;

    v3 N1;
    v3 N2;
    v3 N3;
};


struct shape
{
    union fields
    {
        struct triangle_fields
        {
            triangle Triangle;
            uint8_t padding[MAX_SHAPE_STRUCT_SIZE - sizeof(triangle)];
        } TriangleFields;
        
        struct sphere_fields
        {
            sphere Sphere;
            uint8_t padding[MAX_SHAPE_STRUCT_SIZE - sizeof(sphere)];
        } SphereFields;
    } Fields;
    
    shape_type Type;
    aabb AABB;
    u32 MatIndex;
    mat4 Transform;
};

struct volume
{
public:
    u32 MatIndex;
    aabb AABB;
};

internal shape* Sphere(v3 P, f32 r, u32 MatIndex) {
    shape *Result = (shape*)malloc(sizeof(shape));
    
    Result->Transform = Translate(Identity(), V3(P.x, P.y, P.z));
    
    
    Result->Fields.SphereFields.Sphere.r = r;
    Result->MatIndex = MatIndex;

    Result->AABB = {
        V3(-r + P.x, -r + P.y, -r + P.z),
        V3(r + P.x, r + P.y, r + P.z)
    };

    Result->Type = shape_type::sphereType;

    return Result;
}

internal shape Triangle(v3 Vertex1, v3 Vertex2, v3 Vertex3, v3 Normal1, v3 Normal2, v3 Normal3, mat4 Transform, u32 MatIndex) {
    shape Result = {};
    
    Result.Transform = Transform;

    v3 WorldV1 = TransformPosition(Transform, Vertex1);
    v3 WorldV2 = TransformPosition(Transform, Vertex2);
    v3 WorldV3 = TransformPosition(Transform, Vertex3);
    
    Result.Fields.TriangleFields.Triangle.V1 = WorldV1;
    Result.Fields.TriangleFields.Triangle.V2 = WorldV2;
    Result.Fields.TriangleFields.Triangle.V3 = WorldV3;
    
    Result.Fields.TriangleFields.Triangle.N1 = Normal1;
    Result.Fields.TriangleFields.Triangle.N2 = Normal2;
    Result.Fields.TriangleFields.Triangle.N3 = Normal3;

    Result.MatIndex = MatIndex;


    v3 MinPosition = Min(WorldV1, Min(WorldV2, WorldV3));
    v3 MaxPosition = Max(WorldV1, Max(WorldV2, WorldV3));

    Result.AABB = {
        MinPosition, MaxPosition
    };
    
    Result.Type = shape_type::triangleType;

    return Result;    
}


struct bvh {
public:
    bvh *Left;
    bvh *Right;
    aabb AABB;
    shape *Shapes;
};

struct world 
{
    u32 MaterialCount;
    material *Materials;

    u32 ShapesCount;
    shape *Shapes;

    // u32 VolumeCount;
    // volume *Volumes;

    bvh *BVH;
};

#include "Hit.h"

//THREADING

struct work_order 
{
    world *World;
    image_32 
    Image;
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

    //Memory arenas
    lane_v2 *BrdfSamples[2];
    lane_v2 *SubPixelSamples;
    
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
