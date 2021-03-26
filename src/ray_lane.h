#if !defined(LANE_WIDTH)
#define LANE_WIDTH 1
#endif

struct v3 
{
    f32 x, y, z;
};


#if (LANE_WIDTH == 1)
#include "ray_lane1.h"

#elif (LANE_WIDTH == 8)

#if COMPILER_MSVC
#include <intrin.h>
#elif COMPILER_LLVM
#include <x86intrin.h>
#else
#error "SEE/NEON optimizations are not available for this compiler"
#endif


#include "ray_lane8.h"

#elif (LANE_WIDTH == 4)

#if COMPILER_MSVC
#include <intrin.h>
#elif COMPILER_LLVM
#include <x86intrin.h>
#else
#error "SEE/NEON optimizations are not available for this compiler"
#endif

#include "ray_lane4.h"

#else
#error Lane width must be 1 or 4
#endif


#if (LANE_WIDTH != 1)

internal lane_f32 operator/(lane_f32 A, f32 B) {
    lane_f32 Result = A / LaneF32FromF32(B);
    return Result;
}

internal lane_f32 operator/(f32 A, lane_f32 B) {
    lane_f32 Result = LaneF32FromF32(A)/B;
    return Result;
}

internal lane_f32 operator+(lane_f32 A, f32 B) {
    lane_f32 Result = A + LaneF32FromF32(B);
    return Result;
}

internal lane_f32 operator+(f32 A, lane_f32 B) {
    lane_f32 Result = LaneF32FromF32(A)+B;
    return Result;
}

internal lane_f32 operator-(lane_f32 A, f32 B) {
    lane_f32 Result = A - LaneF32FromF32(B);
    return Result;
}

internal lane_f32 operator-(f32 A, lane_f32 B) {
    lane_f32 Result = LaneF32FromF32(A)-B;
    return Result;
}

internal lane_f32 operator*(lane_f32 A, f32 B) {
    lane_f32 Result = A * LaneF32FromF32(B);
    return Result;
}

internal lane_f32 operator*(f32 A, lane_f32 B) {
    lane_f32 Result = LaneF32FromF32(A)*B;
    return Result;
}



internal lane_f32 operator+=(lane_f32& A, f32 B) {
    A = A + B;
    return A;
}

internal lane_u32 operator+=(lane_u32& A, lane_u32 B) {
    A = A + B;
    return A;
}

internal lane_f32 operator-=(lane_f32& A, f32 B) {
    A = A - B;
    return A;
}

internal lane_f32 operator*=(lane_f32& A, f32 B) {
    A = A * B;
    return A;
}

internal lane_f32 operator/=(lane_f32& A, f32 B) {
    A = A / B;
    return A;
}



internal lane_u32 operator&=(lane_u32& A, lane_u32 B) {
    A = A & B;
    return A;
}

internal lane_u32 operator^=(lane_u32& A, lane_u32 B) {
    A = A ^ B;
    return A;
}

internal lane_u32 operator|=(lane_u32& A, lane_u32 B) {
    A = A | B;
    return A;
}

lane_u32 &lane_u32::operator=(u32 B) {
    *this = LaneU32FromU32(B);
    return *this;
}

lane_f32 &lane_f32::operator=(f32 B) {
    *this = LaneF32FromF32(B);
    return *this;
}

internal lane_f32 operator-(lane_f32 A) {
    lane_f32 Result = LaneF32FromF32(0) - A;
    return Result;
}


internal lane_u32 operator>(lane_f32 A, f32 B) {
    lane_u32 Result = (A > LaneF32FromF32(B));
    return Result;
}

internal lane_u32 operator>=(lane_f32 A, f32 B) {
    lane_u32 Result = (A >= LaneF32FromF32(B));
    return Result;
}

internal lane_u32 operator>(f32 A, lane_f32 B) {
    lane_u32 Result = (LaneF32FromF32(A) > B);
    return Result;
}

internal lane_u32 operator<(lane_f32 A, f32 B) {
    lane_u32 Result = (A < LaneF32FromF32(B));
    return Result;
}

internal lane_u32 operator<=(lane_f32 A, f32 B) {
    lane_u32 Result = (A <= LaneF32FromF32(B));
    return Result;
}

internal lane_u32 operator!=(lane_v3 A, lane_v3 B) {
    lane_u32 Result;
    Result = ((A.x != B.x)) & (A.y != B.y) & (A.z != B.z);
    return Result;
}


internal lane_u32 operator<(f32 A, lane_f32 B) {
    lane_u32 Result = (LaneF32FromF32(A) < B);
    return Result;
}

internal void ConditionalAssign(lane_u32 *Dest, lane_u32 Mask, lane_u32 Source) {
    *Dest = (AndNot(Mask, *Dest)) | (Mask & Source);
}


inline lane_f32 Clamp01(lane_f32 Value) {
    lane_f32 Result = Min( Max(Value, LaneF32FromF32(0.0f)), LaneF32FromF32(1.0f));
    return Result;
}


internal lane_v3
operator&(lane_u32 A, lane_v3 B)
{
    lane_v3 Result;
    
    Result.x = A & B.x;
    Result.y = A & B.y;
    Result.z = A & B.z;
    
    return(Result);
}


internal lane_f32 V3Component(lane_v3 Vector, u32 ComponentIndex) {
    return *( ((lane_f32*)(&Vector)) + ComponentIndex);
}

#endif

internal f32 V3Component(v3 Vector, u32 ComponentIndex) {
    return *( ((f32*)(&Vector)) + ComponentIndex);
}


internal lane_v3 GatherV3_(void *BasePtr, u32 Stride, lane_u32 Indices) {
    lane_v3 Result;
    
    Result.x = GatherF32_((f32*)BasePtr + 0, Stride, Indices);
    Result.y = GatherF32_((f32*)BasePtr + 1, Stride, Indices);
    Result.z = GatherF32_((f32*)BasePtr + 2, Stride, Indices);

    return Result;
}


#define GatherF32(BasePtr, Index, Member) GatherF32_(&(BasePtr->Member), sizeof(*(BasePtr)), Index)
#define GatherU32(BasePtr, Index, Member) GatherU32_(&(BasePtr->Member), sizeof(*(BasePtr)), Index)
#define GatherV3(BasePtr, Index, Member) GatherV3_(&(BasePtr->Member), sizeof(*(BasePtr)), Index)


internal lane_u32 operator==(lane_v3 A, lane_v3 B) {
    lane_u32 Result;
    Result = (A.x == B.x) & (A.y == B.y) & (A.z == B.z);
    return Result;
}

v3 Extract0(lane_v3 A) {
    v3 Result;
    Result.x = *(f32 *)&A.x;
    Result.y = *(f32 *)&A.y;
    Result.z = *(f32 *)&A.z;

    return Result;
}

u32 Extract0(lane_u32 A) {
    u32 Result;
    Result = *(u32 *)&A;
    return Result;
}

f32 Extract0(lane_f32 A) {
    f32 Result;
    Result = *(f32 *)&A;
    return Result;
}


v3 ExtractAtIndex(lane_v3 A, u32 Index) {
    v3 Result;
    Result.x = *(((f32 *)&A.x)+Index);
    Result.y = *(((f32 *)&A.y)+Index);
    Result.z = *(((f32 *)&A.z)+Index);
    return Result;
}

internal void ConditionalAssign(lane_v3 *Dest, lane_u32 Mask, lane_v3 Source) {
    ConditionalAssign(&Dest->x, Mask, Source.x);
    ConditionalAssign(&Dest->y, Mask, Source.y);
    ConditionalAssign(&Dest->z, Mask, Source.z);
}


inline lane_v3 LaneV3FromV3(v3 V) {
    lane_v3 Result;
    Result.x = LaneF32FromF32(V.x);
    Result.y = LaneF32FromF32(V.y);
    Result.z = LaneF32FromF32(V.z);
    return Result;
}

inline lane_v3
V3(f32 X, f32 Y, f32 Z)
{
    lane_v3 Result;
    
    Result.x = X;
    Result.y = Y;
    Result.z = Z;
    
    return(Result);
}

inline lane_v3
operator+(lane_v3 A, lane_v3 B)
{
    lane_v3 Result;
    
    Result.x = A.x + B.x;
    Result.y = A.y + B.y;
    Result.z = A.z + B.z;
    
    return(Result);
}

inline lane_v3
operator-(lane_v3 A, lane_v3 B)
{
    lane_v3 Result;
    
    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    Result.z = A.z - B.z;
    
    return(Result);
}


inline lane_v3
operator-(lane_v3 A)
{
    lane_v3 Result;
    
    Result.x = -A.x;
    Result.y = -A.y;
    Result.z = -A.z;
    
    return(Result);
}



inline lane_v3
Cross(lane_v3 A, lane_v3 B)
{
    lane_v3 Result;
    
    Result.x = A.y*B.z - A.z*B.y;
    Result.y = A.z*B.x - A.x*B.z;
    Result.z = A.x*B.y - A.y*B.x;
    
    return(Result);
}


inline lane_v3
Hadamard(lane_v3 A, lane_v3 B)
{
    lane_v3 Result = {A.x*B.x, A.y*B.y, A.z*B.z};
    
    return(Result);
}


inline lane_v3 &
operator+=(lane_v3 &A, lane_v3 B)
{
    A = A + B;
    
    return(A);
}




internal lane_v3 operator*(lane_v3 A, lane_f32 B) {
    lane_v3 Result;

    Result.x = A.x * B;
    Result.y = A.y * B;
    Result.z = A.z * B;
    return Result;
}

internal lane_v3 operator*(lane_f32 A, lane_v3 B) {
    lane_v3 Result = B * A;
    return Result;
}

inline lane_f32
Inner(lane_v3 A, lane_v3 B)
{
    lane_f32 Result =  A.x*B.x + A.y*B.y + A.z*B.z;
    
    return(Result);
}

inline lane_f32
LengthSq(lane_v3 A)
{
    lane_f32 Result = Inner(A, A);
    
    return(Result);
}

inline lane_f32
Length(lane_v3 A)
{
    lane_f32 Result = SquareRoot(LengthSq(A));
    return(Result);
}

//V3
inline lane_v3
NOZ(lane_v3 A)
{
    lane_v3 Result = {};
    
    lane_f32 LenSq = LengthSq(A);
    lane_u32 Mask = (LenSq > LaneF32FromF32(Square(0.0001f)));
    ConditionalAssign(&Result, Mask, A * (1.0f / SquareRoot(LenSq)));
    
    return(Result);
}


internal v3 HorizontalAdd(lane_v3 A) {
    v3 Result;
    Result.x = HorizontalAdd(A.x);
    Result.y = HorizontalAdd(A.y);
    Result.z = HorizontalAdd(A.z);
    return Result;
}



inline lane_v3 LaneV3(lane_f32 X, lane_f32 Y, lane_f32 Z) {
    lane_v3 Result;
    
    Result.x = X;
    Result.y = Y;
    Result.z = Z;

    return Result;
}

inline lane_v3
Lerp(lane_v3 A, lane_f32 t, lane_v3 B)
{
    lane_v3 Result = (1.0f - t)*A + t*B;
    
    return(Result);
}


internal lane_f32 Pow5(lane_f32 A) {
    lane_f32 Result = A * A * A * A * A;
    return Result;
}

inline lane_v3 Reflect(lane_v3 A, lane_v3 Normal) {
    lane_v3 Result;
    Result = A - 2.0f * Inner(A, Normal) * Normal;
    return Result;
}

inline lane_u32 Refract(lane_v3 IncidentVector, lane_v3 Normal, lane_f32 IndexOfRefraction, lane_v3 *OutRefractedVector) {
    lane_v3 IncidentVectorNormalized = NOZ(IncidentVector);
    lane_f32 IncidentDotNormal = Inner(IncidentVectorNormalized, Normal);

    lane_f32 Discriminant = LaneF32FromF32(1.0f) - IndexOfRefraction * IndexOfRefraction * (LaneF32FromF32(1.0f) - IncidentDotNormal * IncidentDotNormal);

    
    lane_u32 DiscriminantPositiveMask = (Discriminant>0);
    lane_v3 RefractedVector = IndexOfRefraction * (IncidentVectorNormalized - Normal * IncidentDotNormal) - Normal * SquareRoot(Discriminant);
    ConditionalAssign(OutRefractedVector, DiscriminantPositiveMask, RefractedVector);
    
    return DiscriminantPositiveMask;
}


inline lane_f32 ShlickFresnelApproximation(lane_f32 Cosine, lane_f32 IndexOfRefraction) {
    lane_f32 One = LaneF32FromF32(1.0);
    lane_f32 R0 = (One- IndexOfRefraction) / (One + IndexOfRefraction);
    R0 = R0 * R0;
    return R0 + (One - R0) * Pow5(One - Cosine);
}

