#if !defined(LANE_WIDTH)
#define LANE_WIDTH 1
#endif

struct v2 
{
    f32 x, y;
};

struct v3 
{
    f32 x, y, z;
};

struct v4 
{
    f32 x, y, z, w;
};

#if (LANE_WIDTH == 1)

#include "ray_lane1.h"

#elif (LANE_WIDTH == 8)

#if COMPILER_MSVC
#include <immintrin.h>
#include <intrin.h>
#elif COMPILER_LLVM
#include <x86intrin.h>
#else
#error "SEE/NEON optimizations are not available for this compiler"
#endif


#include "ray_lane8.h"

#elif (LANE_WIDTH == 4)

#if COMPILER_MSVC

#include <immintrin.h>
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



internal lane_v3
operator&(lane_u32 A, lane_v3 B)
{
    lane_v3 Result;
    
    Result.x = A & B.x;
    Result.y = A & B.y;
    Result.z = A & B.z;
    
    return(Result);
}


internal lane_f32 LaneV3Component(lane_v3 Vector, u32 ComponentIndex) {
    return *( ((lane_f32*)(&Vector)) + ComponentIndex);
}


inline lane_v3 LaneV3(f32 X, f32 Y, f32 Z)
{
    lane_v3 Result;
    
    Result.x = LaneF32FromF32(X);
    Result.y = LaneF32FromF32(Y);
    Result.z = LaneF32FromF32(Z);
    
    return(Result);
}


inline lane_v4
LaneV4(f32 X, f32 Y, f32 Z, f32 W)
{
    lane_v4 Result;
    
    Result.x = LaneF32FromF32(X);
    Result.y = LaneF32FromF32(Y);
    Result.z = LaneF32FromF32(Z);
    Result.w = LaneF32FromF32(W);
    
    return(Result);
}




// //Normal functions
// inline f32 Clamp01(f32 Value) {
//     f32 Result = Min( Max(Value, 0.0f), 1.0f);
//     return Result;
// }

// inline f32 Clamp(f32 Value,f32 MinValue, f32 MaxValue) {
//     f32 Result = Min( Max(Value, MinValue), MaxValue);
//     return Result;
// }




#endif

// //Normal functions

inline v3
operator*(f32 A, v3 B)
{
    v3 Result;
    
    Result.x = A*B.x;
    Result.y = A*B.y;
    Result.z = A*B.z;
    
    return(Result);
}

inline v3
operator*(v3 B, f32 A)
{
    v3 Result = A*B;
    
    return(Result);
}

inline v3
operator-(v3 A, v3 B)
{
    v3 Result;
    
    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    Result.z = A.z - B.z;
    
    return(Result);
}


inline lane_v3
Lane_Cross(lane_v3 A, lane_v3 B)
{
    lane_v3 Result;
    
    Result.x = A.y*B.z - A.z*B.y;
    Result.y = A.z*B.x - A.x*B.z;
    Result.z = A.x*B.y - A.y*B.x;
    
    return(Result);
}


inline v3 V3(f32 X, f32 Y, f32 Z)
{
    v3 Result;
    
    Result.x = X;
    Result.y = Y;
    Result.z = Z;
    
    return(Result);
}

internal f32 Min(f32 A, f32 B) {
    f32 Result = ((A<B) ? A : B);
    return Result;
}

internal f32 Max(f32 A, f32 B) {
    f32 Result = ((A>B) ? A : B);
    return Result;
}

internal v3 Min(v3 A, v3 B) {
    v3 Result = V3(Min(A.x, B.x), Min(A.y, B.y), Min(A.z, B.z));
    return Result;
}


internal v3 Max(v3 A, v3 B) {
    v3 Result = {Max(A.x, B.x), Max(A.y, B.y), Max(A.z, B.z)};
    return Result;
}

inline v3
Cross(v3 A, v3 B)
{
    v3 Result;
    
    Result.x = A.y*B.z - A.z*B.y;
    Result.y = A.z*B.x - A.x*B.z;
    Result.z = A.x*B.y - A.y*B.x;
    
    return(Result);
}

internal f32 V3Component(v3 Vector, u32 ComponentIndex) {
    return *( ((f32*)(&Vector)) + ComponentIndex);
}
inline v3 operator+(v3 A, v3 B)
{
    v3 Result;
    
    Result.x = A.x + B.x;
    Result.y = A.y + B.y;
    Result.z = A.z + B.z;
    
    return(Result);
}



#include <algorithm>
f32 Abs(f32 In)
{
    return abs(In);
}


v3 Abs(v3 Value)
{
    v3 Result;
    
    Result.x = Abs(Value.x);
    Result.y = Abs(Value.y);
    Result.z = Abs(Value.z);

    return Result;
}



//Vector functions

inline lane_f32 Lane_Clamp01(lane_f32 Value) {
    lane_f32 Result = Lane_Min( Lane_Max(Value, LaneF32FromF32(0.0f)), LaneF32FromF32(1.0f));
    return Result;
}

inline lane_f32 Lane_Clamp(lane_f32 Value,lane_f32 MinValue, lane_f32 MaxValue) {
    lane_f32 Result = Lane_Min( Lane_Max(Value, MinValue), MaxValue);
    return Result;
}

internal lane_v3 GatherLaneV3_(void *BasePtr, u32 Stride, lane_u32 Indices) {
    lane_v3 Result;
    
    Result.x = GatherF32_((f32*)BasePtr + 0, Stride, Indices);
    Result.y = GatherF32_((f32*)BasePtr + 1, Stride, Indices);
    Result.z = GatherF32_((f32*)BasePtr + 2, Stride, Indices);

    return Result;
}

#define GatherF32(BasePtr, Index, Member) GatherF32_(&(BasePtr->Member), sizeof(*(BasePtr)), Index)
#define GatherU32(BasePtr, Index, Member) GatherU32_(&(BasePtr->Member), sizeof(*(BasePtr)), Index)
#define GatherLaneV3(BasePtr, Index, Member) GatherLaneV3_(&(BasePtr->Member), sizeof(*(BasePtr)), Index)

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

f32 ExtractAt(lane_f32 A, u32 Index) {
    f32 Result;
    Result = *((f32 *)&A + Index);
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

internal void ConditionalAssign(lane_v2 *Dest, lane_u32 Mask, lane_v2 Source) {
    ConditionalAssign(&Dest->x, Mask, Source.x);
    ConditionalAssign(&Dest->y, Mask, Source.y);
}


inline lane_v3 LaneV3FromV3(v3 V) {
    lane_v3 Result;
    Result.x = LaneF32FromF32(V.x);
    Result.y = LaneF32FromF32(V.y);
    Result.z = LaneF32FromF32(V.z);
    return Result;
}

inline lane_v2 LaneV2(lane_f32 X, lane_f32 Y)
{
    lane_v2 Result;
    
    Result.x = X;
    Result.y = Y;
    
    return(Result);
}

#if LANE_WIDTH == 1
#pragma optimize( "", off )
#endif
inline lane_v3 LaneV3(lane_f32 X, lane_f32 Y, lane_f32 Z)
{
    lane_v3 Result;
    
    Result.x = X;
    Result.y = Y;
    Result.z = Z;
    
    return(Result);
}
#if LANE_WIDTH == 1
#pragma optimize( "", on )
#endif

inline lane_v4
LaneV4(lane_f32 X, lane_f32 Y, lane_f32 Z, lane_f32 W)
{
    lane_v4 Result;
    
    Result.x = X;
    Result.y = Y;
    Result.z = Z;
    Result.w = W;
    
    return(Result);
}



inline lane_v3 operator+(lane_v3 A, lane_v3 B)
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
Lane_Hadamard(lane_v3 A, lane_v3 B)
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
Lane_Inner(lane_v3 A, lane_v3 B)
{
    lane_f32 Result =  A.x*B.x + A.y*B.y + A.z*B.z;
    
    return(Result);
}

inline f32
Inner(v3 A, v3 B)
{
    f32 Result =  A.x*B.x + A.y*B.y + A.z*B.z;
    
    return(Result);
}

inline lane_f32
Lane_LengthSq(lane_v3 A)
{
    lane_f32 Result = Lane_Inner(A, A);
    
    return(Result);
}

inline lane_f32
Lane_DistanceSq(lane_v3 A, lane_v3 B)
{
    lane_f32 Result = Lane_LengthSq(B-A);
    
    return(Result);
}

inline f32
LengthSq(v3 A)
{
    f32 Result = Inner(A, A);
    
    return(Result);
}

inline lane_f32
Lane_Length(lane_v3 A)
{
    lane_f32 Result = Lane_SquareRoot(Lane_LengthSq(A));
    return(Result);
}

//LaneV3
inline lane_v3
Lane_NOZ(lane_v3 A)
{
    lane_v3 Result = {};
    
    lane_f32 LenSq = Lane_LengthSq(A);
    lane_u32 Mask = (LenSq > LaneF32FromF32(Square(0.0001f)));
    ConditionalAssign(&Result, Mask, A * (1.0f / Lane_SquareRoot(LenSq)));
    
    return(Result);
}

//LaneV3
inline v3
NOZ(v3 A)
{
    v3 Result = {};
    
    f32 LenSq = LengthSq(A);
    
    if(LenSq > Square(0.0001f))
        Result = A * (1.0f / SquareRoot(LenSq));
    
    return(Result);
}


internal v3 HorizontalAdd(lane_v3 A) {
    v3 Result;
    Result.x = HorizontalAdd(A.x);
    Result.y = HorizontalAdd(A.y);
    Result.z = HorizontalAdd(A.z);
    return Result;
}



inline lane_v3 LaneLaneV3(lane_f32 X, lane_f32 Y, lane_f32 Z) {
    lane_v3 Result;
    
    Result.x = X;
    Result.y = Y;
    Result.z = Z;

    return Result;
}

inline lane_v3
Lane_Lerp(lane_v3 A, lane_f32 t, lane_v3 B)
{
    lane_v3 Result = (1.0f - t)*A + t*B;
    
    return(Result);
}


internal lane_f32 Lane_Pow5(lane_f32 A) {
    lane_f32 Result = A * A * A * A * A;
    return Result;
}

inline lane_v3 Lane_Reflect(lane_v3 A, lane_v3 Normal) {
    lane_v3 Result;
    Result = A - 2.0f * Lane_Inner(A, Normal) * Normal;
    return Result;
}

inline lane_u32 Lane_Refract(lane_v3 IncidentVector, lane_v3 Normal, lane_f32 IndexOfRefraction, lane_v3 *OutRefractedVector) {
    lane_v3 IncidentVectorNormalized = Lane_NOZ(IncidentVector);
    lane_f32 IncidentDotNormal = Lane_Inner(IncidentVectorNormalized, Normal);

    lane_f32 Discriminant = LaneF32FromF32(1.0f) - IndexOfRefraction * IndexOfRefraction * (LaneF32FromF32(1.0f) - IncidentDotNormal * IncidentDotNormal);

    
    lane_u32 DiscriminantPositiveMask = (Discriminant>0);
    lane_v3 RefractedVector = IndexOfRefraction * (IncidentVectorNormalized - Normal * IncidentDotNormal) - Normal * Lane_SquareRoot(Discriminant);
    ConditionalAssign(OutRefractedVector, DiscriminantPositiveMask, RefractedVector);
    
    return DiscriminantPositiveMask;
}


inline lane_f32 Lane_ShlickFresnelApproximation(lane_f32 Cosine, lane_f32 IndexOfRefraction) {
    lane_f32 One = LaneF32FromF32(1.0);
    lane_f32 R0 = (One- IndexOfRefraction) / (One + IndexOfRefraction);
    R0 = R0 * R0;
    return R0 + (One - R0) * Lane_Pow5(One - Cosine);
}


internal lane_v3 Lane_Min(lane_v3 A, lane_v3 B) {
    lane_v3 Result = LaneV3(Lane_Min(A.x, B.x), Lane_Min(A.y, B.y), Lane_Min(A.z, B.z));
    return Result;
}

internal lane_v3 Lane_Max(lane_v3 A, lane_v3 B) {
    lane_v3 Result = LaneV3(Lane_Max(A.x, B.x), Lane_Max(A.y, B.y), Lane_Max(A.z, B.z));
    return Result;
}


inline v4
V4(f32 X, f32 Y, f32 Z, f32 W)
{
    v4 Result;
    
    Result.x = X;
    Result.y = Y;
    Result.z = Z;
    Result.w = W;
    
    return(Result);
}



inline u32
RGBAPack4x8(v4 Unpacked)
{
    u32 Result = ((Roundf32ToUint32(Unpacked.x) << 24) |
                  (Roundf32ToUint32(Unpacked.y) << 16) |
                  (Roundf32ToUint32(Unpacked.z) << 8) |
                  (Roundf32ToUint32(Unpacked.w) << 0));
    
    return(Result);
}

inline u32
ARGBPack4x8(v4 Unpacked)
{
    u32 Result = (
                  (Roundf32ToUint32(Unpacked.z) << 0)  |
                  (Roundf32ToUint32(Unpacked.y) << 8)  |
                  (Roundf32ToUint32(Unpacked.x) << 16) |
                  (Roundf32ToUint32(Unpacked.w) << 24)
                );
    
    return(Result);
}


inline v4
Linear1ToSRGB255(v4 C)
{
    v4 Result;
    
    f32 One255 = 255.0f;
    
    Result.x = One255*SquareRoot(C.x);
    Result.y = One255*SquareRoot(C.y);
    Result.z = One255*SquareRoot(C.z);
    Result.w = One255*C.w;
    
    return(Result);
}

inline v4 V4(v3 xyz, f32 W)
{
    v4 Result;
    
    Result.x = xyz.x;
    Result.y = xyz.y;
    Result.z = xyz.z;
    Result.w = W;
    
    return(Result);
}

inline v4
operator+(v4 A, v4 B)
{
    v4 Result;
    
    Result.x = A.x + B.x;
    Result.y = A.y + B.y;
    Result.z = A.z + B.z;
    Result.w = A.w + B.w;
    
    return(Result);
}

inline v4
operator-(v4 A, v4 B)
{
    v4 Result;
    
    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    Result.z = A.z - B.z;
    Result.w = A.w - B.w;
    
    return(Result);
}

inline v4
operator*(f32 A, v4 B)
{
    v4 Result;
    
    Result.x = A*B.x;
    Result.y = A*B.y;
    Result.z = A*B.z;
    Result.w = A*B.w;
    
    return(Result);
}



inline v4
operator*(v4 B, f32 A)
{
    v4 Result = A*B;
    
    return(Result);
}

inline lane_v2
operator*(lane_f32 A, lane_v2 B)
{
    lane_v2 Result;
    
    Result.x = A*B.x;
    Result.y = A*B.y;
    
    return(Result);
}


inline lane_v2
operator-(lane_v2 A, lane_v2 B)
{
    lane_v2 Result;
    
    Result.x = A.x - B.x;
    Result.y = A.y - B.y;
    
    return(Result);
}

inline lane_v2
operator*(lane_v2 B, lane_f32 A)
{
    lane_v2 Result = A*B;
    
    return(Result);
}


inline f32
Lane_Inner(v4 A, v4 B)
{
    f32 Result = A.x*B.x + A.y*B.y + A.z*B.z+ A.w*B.w;
    
    return(Result);
}


inline f32
Lane_LengthSq(v4 A)
{
    f32 Result = Lane_Inner(A, A);
    
    return(Result);
}

inline v4
Lane_NOZ(v4 A)
{
    v4 Result = {};
    
    f32 LenSq = Lane_LengthSq(A);
    if(LenSq > Square(0.0001f))
    {
        Result = A * (1.0f / SquareRoot(LenSq));
    }
    
    return(Result);
}

lane_v3 Lane_Abs(lane_v3 Value)
{
    lane_v3 Result;
    
    Result.x = Lane_Abs(Value.x);
    Result.y = Lane_Abs(Value.y);
    Result.z = Lane_Abs(Value.z);

    return Result;
}


#if LANE_WIDTH != 1
#endif
#include "matrix.h"
#include "lane_matrix.h"

