typedef f32 lane_f32; 
typedef u32 lane_u32; 


struct lane_v2 
{
    lane_f32 x, y;
};

struct lane_v3 
{
    lane_f32 x, y, z;
};

struct lane_v4 
{
    lane_f32 x, y, z, w;
};

internal lane_f32 LaneV3Component(lane_v3 Vector, lane_u32 ComponentIndex) {
    return *( ((lane_f32*)(&Vector)) + ComponentIndex);
}

internal void ConditionalAssign(lane_u32 *Dest, lane_u32 Mask, lane_u32 Source) {
    Mask = Mask ? 0xFFFFFFFF : 0;
    *Dest = ((~Mask & *Dest) | (Mask & Source));
}

internal void ConditionalAssign(lane_f32 *Dest, lane_u32 Mask, lane_f32 Source) {
    ConditionalAssign((lane_u32*)Dest, Mask, *(lane_u32*)&Source);
}


internal lane_f32 Max(lane_f32 A, lane_f32 B) {
    lane_f32 Result = ((A>B) ? A : B);
    return Result;
}

internal lane_f32 Abs(lane_f32 A) {
    lane_f32 Result = (A < 0) ? -A : A;
    return Result;
}

internal lane_f32 Min(lane_f32 A, lane_f32 B) {
    lane_f32 Result = ((A<B) ? A : B);
    return Result;
}

internal b32 MaskIsZeroed(lane_u32 LaneMask) {
    b32 Result = (LaneMask == 0);
    return Result;
}

internal f32 HorizontalAdd(lane_f32 A) {
    f32 Result = A;
    return Result;
}

internal u32 HorizontalAdd(lane_u32 A) {
    u32 Result = A;
    return Result;
}


internal lane_f32 LaneF32FromU32(lane_u32 V) {
    lane_f32 Result = (lane_f32)V;
    return Result;
}
internal lane_f32 LaneF32FromF32(f32 Repl) {
    lane_f32 Result = Repl;
    return Result;
}

internal lane_u32 LaneU32FromF32(lane_f32 Repl) {
    lane_u32 Result = (lane_u32)Repl;
    return Result;
}


internal lane_u32 LaneU32FromU32(u32 Repl)  {
    lane_u32 Result = Repl;
    return Result;
}


internal lane_v3 operator&(lane_u32 A, lane_v3 B) 
{
    lane_v3 Result;
    A = (A ? 0xFFFFFFFF : 0);
    
    u32 x = A & *(u32 *)&B.x;
    u32 y = A & *(u32 *)&B.y;
    u32 z = A & *(u32 *)&B.z;
    
    Result.x = *(f32 *)&x;
    Result.y = *(f32 *)&y;
    Result.z = *(f32 *)&z;
    return Result;
}

internal lane_u32 LaneU32FromU32(u32 R0,u32 R1,u32 R2,u32 R3) {
    lane_u32 Result = R0;
    return Result;
}


internal lane_f32 GatherF32_(void *BasePtr, u32 Stride, lane_u32 Index) {
    lane_f32 Result = *(f32*)((u8*)BasePtr+ Index * Stride);
    return Result;
}


internal lane_f32 GatherU32_(void *BasePtr, u32 Stride, lane_u32 Index) {
    lane_u32 Result = *(u32*)((u8*)BasePtr+ Index * Stride);
    return Result;
}

// inline lane_v3
// Lerp(lane_v3 A, lane_f32 t, lane_v3 B)
// {
//     lane_v3 Result = (1.0f - t)*A + t*B;
    
//     return(Result);
// }

#include <cmath>
internal lane_f32 Log(lane_f32 Input)
{
    return std::log(Input);
}

internal lane_f32 Exp(lane_f32 Input)
{
    return std::exp(Input);
}

internal lane_f32 Cos(lane_f32 Input)
{
    return cos(Input);
}

internal lane_f32 Sine(lane_f32 Input)
{
    return sin(Input);
}

internal lane_f32 Tangent(lane_f32 Input)
{
    return tan(Input);
}


internal lane_f32 Floor(lane_f32 Input)
{
    return floor(Input);
}

lane_u32 AndNot(lane_u32 A, lane_u32 B=0) 
{
    return !A;
}

lane_u32 Not(lane_u32 A) 
{
    return !A;
}

