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

#define InvPi32 0.3183098861f
#define Pi32 3.14159265359f
#define Tau32 6.28318530717958647692f
#define PiOver2 1.57079632679489661923
#define PiOver4 0.78539816339744830961
#define Sqrt2   1.41421356237309504880

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
//Packing Operations
///////////////////

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

#include "math.h"
