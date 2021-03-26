///////////////////
//TYPES
///////////////////
#include <float.h>

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int16_t  s16;
typedef int32_t  s32;
typedef int64_t  s64;

typedef float f32;
typedef double f64;

#define Pi32 3.14159265359f
#define Tau32 6.28318530717958647692f
#define F32Max FLT_MAX
#define F32Min -FLT_MAX


///////////////////
//SCALAR
///////////////////
inline f32
Square(f32 A)
{
    f32 Result = A*A;
    
    return(Result);
}

#include <math.h>
inline f32 SquareRoot(f32 A) {
    f32 Result = (f32)sqrt(A);
    return Result;
}

inline f32 Pow(f32 A,f32 B) {
    f32 Result = (f32)pow(A, B);
    return Result;
}

inline u32 Roundf32ToUint32(f32 F) 
{
    u32 Result = (u32)(F + 0.5f);
    return Result;
}

///////////////////
//VECTOR4
///////////////////
struct v4 {
  f32 x, y, z, w;  
};

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

// inline v4
// V4(lane_v3 xyz, f32 W)
// {
//     v4 Result;
    
//     Result.x = xyz.x;
//     Result.y = xyz.y;
//     Result.z = xyz.z;
//     Result.w = W;
    
//     return(Result);
// }

// inline v4
// operator+(v4 A, v4 B)
// {
//     v4 Result;
    
//     Result.x = A.x + B.x;
//     Result.y = A.y + B.y;
//     Result.z = A.z + B.z;
//     Result.w = A.w + B.w;
    
//     return(Result);
// }

// inline v4
// operator-(v4 A, v4 B)
// {
//     v4 Result;
    
//     Result.x = A.x - B.x;
//     Result.y = A.y - B.y;
//     Result.z = A.z - B.z;
//     Result.w = A.w - B.w;
    
//     return(Result);
// }

// inline v4
// operator*(f32 A, v4 B)
// {
//     v4 Result;
    
//     Result.x = A*B.x;
//     Result.y = A*B.y;
//     Result.z = A*B.z;
//     Result.w = A*B.w;
    
//     return(Result);
// }

// inline v4
// operator*(v4 B, f32 A)
// {
//     v4 Result = A*B;
    
//     return(Result);
// }


// inline f32
// Inner(v4 A, v4 B)
// {
//     f32 Result = A.x*B.x + A.y*B.y + A.z*B.z+ A.w*B.w;
    
//     return(Result);
// }


// inline f32
// LengthSq(v4 A)
// {
//     f32 Result = Inner(A, A);
    
//     return(Result);
// }

// inline v4
// NOZ(v4 A)
// {
//     v4 Result = {};
    
//     f32 LenSq = LengthSq(A);
//     if(LenSq > Square(0.0001f))
//     {
//         Result = A * (1.0f / SquareRoot(LenSq));
//     }
    
//     return(Result);
// }


///////////////////
//Packing Operations
///////////////////

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

inline u32 NextPowerOfTwo(u32 A) {
    A--;
    A |= A >> 1;
    A |= A >> 2;
    A |= A >> 4;
    A |= A >> 8;
    A |= A >> 16;
    A++;
    return A;
}
