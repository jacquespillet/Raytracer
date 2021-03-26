
//CANNOT TEST THAT.

struct lane_f32 {
    __m512 V;
    lane_f32 &operator=(f32 A);
};

struct lane_u32 {
    __m512i V;
    lane_u32 &operator=(u32 A);
};

struct lane_v3 {
  lane_f32 x, y, z;  
};

internal lane_u32 operator<(lane_f32 A, lane_f32 B) {
    lane_u32 Result;

    __mmask16 ComparisonMask = _mm512_cmp_ps_mask (A.V, B.V, _CMP_LT_OQ);
    Result.V = _mm512_mask_set1_epi32 (_mm512_set1_epi32(0), ComparisonMask, 0xFFFFFFFF);

    return Result;
}

internal lane_u32 operator<=(lane_f32 A, lane_f32 B) {
    lane_u32 Result;

    __mmask16 ComparisonMask = _mm512_cmp_ps_mask (A.V, B.V, _CMP_LE_OQ);
    Result.V = _mm512_mask_set1_epi32 (_mm512_set1_epi32(0), ComparisonMask, 0xFFFFFFFF);

    return Result;
}


internal lane_u32 operator>(lane_f32 A, lane_f32 B) {
    lane_u32 Result;

    __mmask16 ComparisonMask = _mm512_cmp_ps_mask (A.V, B.V, _CMP_GT_OQ);
    Result.V = _mm512_mask_set1_epi32 (_mm512_set1_epi32(0), ComparisonMask, 0xFFFFFFFF);

    return Result;
}

internal lane_u32 operator>=(lane_f32 A, lane_f32 B) {
    lane_u32 Result;

    __mmask16 ComparisonMask = _mm512_cmp_ps_mask (A.V, B.V, _CMP_GE_OQ);
    Result.V = _mm512_mask_set1_epi32 (_mm512_set1_epi32(0), ComparisonMask, 0xFFFFFFFF);

    return Result;
}


internal lane_u32 operator==(lane_f32 A, lane_f32 B) {
    lane_u32 Result;

    __mmask16 ComparisonMask = _mm512_cmp_ps_mask (A.V, B.V, _CMP_EQ_OQ);
    Result.V = _mm512_mask_set1_epi32 (_mm512_set1_epi32(0), ComparisonMask, 0xFFFFFFFF);

    return Result;
}


internal lane_u32 operator!=(lane_f32 A, lane_f32 B) {
    lane_u32 Result;

    __mmask16 ComparisonMask = _mm512_cmp_ps_mask (A.V, B.V, _CMP_NEQ_OQ);
    Result.V = _mm512_mask_set1_epi32 (_mm512_set1_epi32(0), ComparisonMask, 0xFFFFFFFF);
    return Result;
}

internal lane_u32 operator!=(lane_u32 A, lane_u32 B) {
    lane_u32 Result;

    __mmask16 ComparisonMask = _mm512_cmp_epi32_mask (A.V, B.V, _CMP_EQ_OQ);
    __m512i EqualityMask = _mm512_mask_set1_epi32 (_mm512_set1_epi32(0), ComparisonMask, 0xFFFFFFFF);
    
    Result.V = _mm512_xor_si512(EqualityMask, _mm512_set1_epi32(0xFFFFFFFF));

    return Result;
}


internal lane_u32 operator^(lane_u32 A, lane_u32 B) 
{
    lane_u32 Result;
    Result.V = _mm512_xor_si512(A.V, B.V);
    return Result;
}

internal lane_u32 operator&(lane_u32 A, lane_u32 B) 
{
    lane_u32 Result;
    Result.V = _mm512_and_si512(A.V, B.V);
    return Result;
}

internal lane_f32 operator&(lane_u32 A, lane_f32 B) 
{
    lane_f32 Result;
    Result.V = _mm512_and_ps(_mm512_castsi512_ps(A.V), B.V);
    return Result;
}

internal lane_u32 operator|(lane_u32 A, lane_u32 B) 
{
    lane_u32 Result;
    Result.V = _mm512_or_si512(A.V, B.V);
    return Result;
}


internal lane_u32 operator<<(lane_u32 A, u32 Shift) {
    lane_u32 Result;
    Result.V = _mm512_slli_epi32(A.V, Shift);
    return Result;
}

internal lane_u32 operator>>(lane_u32 A, u32 Shift) {
    lane_u32 Result;
    Result.V = _mm512_srli_epi32(A.V, Shift);
    return Result;
}

internal lane_u32 AndNot(lane_u32 A, lane_u32 B) 
{
    lane_u32 Result;
    Result.V = _mm512_andnot_si512(A.V, B.V);
    return Result;
}

internal lane_f32 operator/(lane_f32 A, lane_f32 B) {
    lane_f32 Result;
    Result.V = _mm512_div_ps(A.V, B.V);
    return Result;
}

internal lane_f32 operator+(lane_f32 A, lane_f32 B) {
    lane_f32 Result;
    Result.V = _mm512_add_ps(A.V, B.V);
    return Result;
}

internal lane_u32 operator+(lane_u32 A, lane_u32 B) {
    lane_u32 Result;
    Result.V = _mm512_add_epi32(A.V, B.V);
    return Result;
}

internal lane_f32 operator-(lane_f32 A, lane_f32 B) {
    lane_f32 Result;
    Result.V = _mm512_sub_ps(A.V, B.V);
    return Result;
}

internal lane_f32 operator*(lane_f32 A, lane_f32 B) {
    lane_f32 Result;
    Result.V = _mm512_mul_ps(A.V, B.V);
    return Result;
}



internal lane_f32 LaneF32FromU32(lane_u32 A) {
    lane_f32 Result;
    Result.V = _mm512_cvtepi32_ps(A.V);
    return Result;
}


internal lane_f32 LaneF32FromU32(u32 Repl) {
    lane_f32 Result;
    Result.V = _mm512_set1_ps((f32)Repl);
    return Result;
}

internal lane_f32 LaneF32FromF32(f32 Repl) {
    lane_f32 Result;
    Result.V = _mm512_set1_ps(Repl);
    return Result;
}

internal lane_u32 LaneU32FromU32(u32 Repl) {
    lane_u32 Result;
    Result.V = _mm512_set1_epi32(Repl);
    return Result;
}

internal lane_u32 LaneU32FromU32(u32 R0,u32 R1,u32 R2,u32 R3,u32 R4,u32 R5,u32 R6,u32 R7, u32 R8,u32 R9,u32 R10,u32 R11,u32 R12,u32 R13,u32 R14,u32 R15) {
    lane_u32 Result;
    Result.V = _mm512_setr_epi32(R0, R1, R2, R3, R4, R5, R6, R7,R8, R9, R10, R11, R12, R13, R14, R15);
    return Result;
}

internal lane_f32 SquareRoot(lane_f32 A) {
    lane_f32 Result;
    Result.V = _mm512_sqrt_ps(A.V);
    return Result;
}

internal void ConditionalAssign(lane_f32 *Dest, lane_u32 Mask, lane_f32 Source) {
    __m512 MaskPS =  _mm512_castsi512_ps(Mask.V);
    Dest->V = _mm512_or_ps(
                    _mm512_andnot_ps(MaskPS, Dest->V),
                    _mm512_and_ps(MaskPS, Source.V)
                    );
}

internal lane_f32 Min(lane_f32 A, lane_f32 B) {
    lane_f32 Result;
    Result.V = _mm512_min_ps(A.V, B.V);
    return Result;
}

internal lane_f32 Max(lane_f32 A, lane_f32 B) {
    lane_f32 Result;
    Result.V = _mm512_max_ps(A.V, B.V);
    return Result;
}

internal lane_f32 GatherF32_(void *BasePtr, u32 Stride, lane_u32 Indices) {
    lane_f32 Result;
    u32* V = (u32*) &Indices.V;
    Result.V = _mm512_setr_ps( 
                                    *(f32*)((u8*)BasePtr+ V[0] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[1] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[2] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[3] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[4] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[5] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[6] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[7] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[8] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[9] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[10] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[11] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[12] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[13] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[14] * Stride),
                                    *(f32*)((u8*)BasePtr+ V[15] * Stride)
                                );
    return Result;
}


internal b32 MaskIsZeroed(lane_u32 A) {
    __mmask16 Mask = _mm512_movepi32_mask(A.V);
    return (Mask==0);
}

internal f32 HorizontalAdd(lane_f32 A) {
    f32 *V = (f32 *)&A.V;
    f32 Result = V[0] +V[1] +V[2] +V[3]; 
    return Result;
}

internal u32 HorizontalAdd(lane_u32 A) {
    u32 *V = (u32 *)&A.V;
    u64 Result = (u64)V[0] +(u64)V[1] +(u64)V[2] +(u64)V[3] +(u64)V[4] +(u64)V[5] +(u64)V[6] +(u64)V[7];
    return (u32)Result;
}