
struct lane_f32 {
    __m128 V;
    lane_f32 &operator=(f32 A);
};

struct lane_u32 {
    __m128i V;
    lane_u32 &operator=(u32 A);
};


struct lane_v2 {
  lane_f32 x, y;  
};

struct lane_v3 {
  lane_f32 x, y, z;  
};

struct lane_v4 {
  lane_f32 x, y, z, w;  
};

internal lane_u32 operator<<(lane_u32 A, u32 Shift) {
    lane_u32 Result;
    Result.V = _mm_slli_epi32(A.V, Shift);
    return Result;
}

internal lane_u32 operator>>(lane_u32 A, u32 Shift) {
    lane_u32 Result;
    Result.V = _mm_srli_epi32(A.V, Shift);
    return Result;
}

internal lane_u32 operator&(lane_u32 A, lane_u32 B) 
{
    lane_u32 Result;
    Result.V = _mm_and_si128(A.V, B.V);
    return Result;
}

internal lane_f32 operator&(lane_u32 A, lane_f32 B) 
{
    lane_f32 Result;
    Result.V = _mm_and_ps(_mm_castsi128_ps(A.V), B.V);
    return Result;
}

internal lane_u32 AndNot(lane_u32 A, lane_u32 B) 
{
    lane_u32 Result;
    Result.V = _mm_andnot_si128(A.V, B.V);
    return Result;
}

internal lane_u32 operator|(lane_u32 A, lane_u32 B) 
{
    lane_u32 Result;
    Result.V = _mm_or_si128(A.V, B.V);
    return Result;
}

internal lane_u32 operator^(lane_u32 A, lane_u32 B) 
{
    lane_u32 Result;
    Result.V = _mm_xor_si128(A.V, B.V);
    return Result;
}


internal lane_f32 LaneF32FromU32(lane_u32 A) {
    lane_f32 Result;
    Result.V = _mm_cvtepi32_ps(A.V);
    return Result;
}


internal lane_f32 LaneF32FromU32(u32 Repl) {
    lane_f32 Result;
    Result.V = _mm_set1_ps((f32)Repl);
    return Result;
}

internal lane_f32 LaneF32FromF32(f32 Repl) {
    lane_f32 Result;
    Result.V = _mm_set1_ps(Repl);
    return Result;
}

internal lane_f32 LaneF32FromF32(f32 R0,f32 R1,f32 R2,f32 R3) {
    lane_f32 Result;
    Result.V = _mm_setr_ps(R0, R1, R2, R3);
    return Result;
}

internal lane_u32 LaneU32FromU32(u32 Repl) {
    lane_u32 Result;
    Result.V = _mm_set1_epi32(Repl);
    return Result;
}

internal lane_u32 LaneU32FromU32(u32 R0,u32 R1,u32 R2,u32 R3) {
    lane_u32 Result;
    Result.V = _mm_setr_epi32(R0, R1, R2, R3);
    return Result;
}

internal lane_f32 operator/(lane_f32 A, lane_f32 B) {
    lane_f32 Result;
    Result.V = _mm_div_ps(A.V, B.V);
    return Result;
}

internal lane_f32 operator+(lane_f32 A, lane_f32 B) {
    lane_f32 Result;
    Result.V = _mm_add_ps(A.V, B.V);
    return Result;
}

internal lane_u32 operator+(lane_u32 A, lane_u32 B) {
    lane_u32 Result;
    Result.V = _mm_add_epi32(A.V, B.V);
    return Result;
}

internal lane_f32 operator-(lane_f32 A, lane_f32 B) {
    lane_f32 Result;
    Result.V = _mm_sub_ps(A.V, B.V);
    return Result;
}

internal lane_f32 operator*(lane_f32 A, lane_f32 B) {
    lane_f32 Result;
    Result.V = _mm_mul_ps(A.V, B.V);
    return Result;
}

internal lane_u32 operator*(lane_u32 A, lane_u32 B) {
    lane_u32 Result;
    Result.V = _mm_mul_epi32(A.V, B.V);
    return Result;
}

internal lane_f32 SquareRoot(lane_f32 A) {
    lane_f32 Result;
    Result.V = _mm_sqrt_ps(A.V);
    return Result;
}

internal lane_u32 operator<(lane_f32 A, lane_f32 B) {
    lane_u32 Result;

    Result.V = _mm_castps_si128(_mm_cmplt_ps(A.V, B.V));

    return Result;
}

internal lane_u32 operator>(lane_f32 A, lane_f32 B) {
    lane_u32 Result;

    Result.V = _mm_castps_si128(_mm_cmpgt_ps(A.V, B.V));

    return Result;
}

internal lane_u32 operator<=(lane_f32 A, lane_f32 B) {
    lane_u32 Result;

    Result.V = _mm_castps_si128(_mm_cmple_ps(A.V, B.V));

    return Result;
}

internal lane_u32 operator>=(lane_f32 A, lane_f32 B) {
    lane_u32 Result;

    Result.V = _mm_castps_si128(_mm_cmpge_ps(A.V, B.V));

    return Result;
}

internal lane_u32 operator==(lane_f32 A, lane_f32 B) {
    lane_u32 Result;

    Result.V = _mm_castps_si128(_mm_cmpeq_ps(A.V, B.V));

    return Result;
}

internal lane_u32 operator!=(lane_f32 A, lane_f32 B) {
    lane_u32 Result;

    Result.V = _mm_castps_si128(_mm_cmpneq_ps(A.V, B.V));

    return Result;
}

internal lane_u32 operator!=(lane_u32 A, lane_u32 B) {
    lane_u32 Result;

    Result.V = _mm_xor_si128(_mm_cmpeq_epi32(A.V, B.V), _mm_set1_epi32(0xFFFFFFFF));

    return Result;
}


internal lane_u32 operator==(lane_u32 A, lane_u32 B) {
    lane_u32 Result;

    Result.V = _mm_cmpeq_epi32(A.V, B.V);

    return Result;
}


internal void ConditionalAssign(lane_f32 *Dest, lane_u32 Mask, lane_f32 Source) {
    __m128 MaskPS =  _mm_castsi128_ps(Mask.V);
    Dest->V = _mm_or_ps(
                    _mm_andnot_ps(MaskPS, Dest->V),
                    _mm_and_ps(MaskPS, Source.V)
                    );
}

internal lane_f32 Min(lane_f32 A, lane_f32 B) {
    lane_f32 Result;
    Result.V = _mm_min_ps(A.V, B.V);
    return Result;
}

internal lane_f32 Max(lane_f32 A, lane_f32 B) {
    lane_f32 Result;
    Result.V = _mm_max_ps(A.V, B.V);
    return Result;
}

internal lane_f32 GatherF32_(void *BasePtr, u32 Stride, lane_u32 Indices) {
    lane_f32 Result;
    u32* V = (u32*) &Indices.V;
    Result.V = _mm_setr_ps( 
                                    *(f32*)((u8*)BasePtr+ V[0] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[1] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[2] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[3] * Stride)
                                );
    return Result;
}


internal lane_u32 GatherU32_(void *BasePtr, u32 Stride, lane_u32 Indices) {
    lane_u32 Result;
    u32* V = (u32*) &Indices.V;
    Result.V = _mm_setr_epi32( 
                                *(u32*)((u8*)BasePtr+ V[0] * Stride),
                                *(u32*)((u8*)BasePtr+ V[1] * Stride),
                                *(u32*)((u8*)BasePtr+ V[2] * Stride),
                                *(u32*)((u8*)BasePtr+ V[3] * Stride)
                                );
    return Result;
}

internal b32 MaskIsZeroed(lane_u32 A) {
    int Mask = _mm_movemask_epi8(A.V);
    return (Mask==0);
}

internal f32 HorizontalAdd(lane_f32 A) {
    f32 *V = (f32 *)&A.V;
    f32 Result = V[0] +V[1] +V[2] +V[3]; 
    return Result;
}

internal u32 HorizontalAdd(lane_u32 A) {
    u32 *V = (u32 *)&A.V;
    u64 Result = (u64)V[0] +(u64)V[1] +(u64)V[2] +(u64)V[3]; 
    return (u32)Result;
}

internal lane_f32 Abs(lane_f32 A) {
    lane_u32 Mask = A <= LaneF32FromF32(0.0f);
    lane_f32 Result = A;
    ConditionalAssign(&Result, Mask, LaneF32FromF32(-1.0f) * A);
    return Result;
}