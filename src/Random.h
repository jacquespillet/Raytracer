
internal lane_u32 XorShift32(random_series *Series) {
	lane_u32 x = Series->State;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
    Series->State = x;
	return x;    
}

internal lane_f32 RandomUnilateral(random_series *Series) {

    lane_f32 Result = LaneF32FromU32(XorShift32(Series)>>1) / (f32)(U32Max>>1);
    return Result;
}

internal lane_f32 RandomBilateral(random_series *Series) {
    lane_f32 Result = (2.0f * RandomUnilateral(Series)) - 1.0f;
    return Result;
}

internal lane_v3 RandomInSphere(random_series *Series) {
    lane_v3 Result;

    lane_u32 IsLengthLessThanOneMask = LaneU32FromU32(0x00000000);

    while(MaskIsZeroed(IsLengthLessThanOneMask)) {
        Result = {RandomBilateral(Series), RandomBilateral(Series), RandomBilateral(Series)};

        IsLengthLessThanOneMask = (Lane_LengthSq(Result) <= 1);
    }
    return Result;
}

internal lane_v3 RandomInDisk(random_series *Series) {
    lane_v3 Result;

    lane_u32 IsLengthLessThanOneMask = LaneU32FromU32(0x00000000);

    while(MaskIsZeroed(IsLengthLessThanOneMask)) {
        Result = {RandomBilateral(Series), RandomBilateral(Series), 0};
        IsLengthLessThanOneMask = (Lane_LengthSq(Result) <= 1);
    }
    return Result;
}

internal f32 RandomUnilateralSlow(){
    return (f32)rand() / (f32)RAND_MAX;
}

internal f32 RandomBilateralSlow(){
    f32 Result = (2.0f * RandomUnilateralSlow()) - 1.0f;
    return Result;
}

internal f32 RandomBilateralSlow(f32 min, f32 max){
    f32 range = max-min;
    f32 Result = (range * RandomUnilateralSlow()) - min;
    return Result;
}


///IMPORTANCE SAMPLING

internal void Jitter(lane_v2 *WideSamples, u32 NumSamples, random_series *Series) {
    v2 *Samples = (v2*) malloc(NumSamples * sizeof(v2));

    u32 Width = (u32)SquareRoot(NumSamples);
    f32 CellSize = 1.0f / (f32)(Width);

    for(u32 i=0; i<Width;i++) {
        for(u32 j=0; j<Width;j++) {
            Samples[i * Width + j].x = i * CellSize + Extract0(RandomUnilateral(Series)) * CellSize;
            Samples[i * Width + j].y = j * CellSize + Extract0(RandomUnilateral(Series)) * CellSize;
        }
    }

        
    if(Width * Width < NumSamples)
    {
        for(u32 i=Width * Width-1; i<NumSamples; i++)
        {
            Samples[i].x = Extract0(RandomUnilateral(Series));
            Samples[i].y = Extract0(RandomUnilateral(Series));
        }
    }

    for(u32 SingleSampleIndex=0; SingleSampleIndex<NumSamples; SingleSampleIndex+= LANE_WIDTH) {
#if(LANE_WIDTH == 1)
        
        (*WideSamples).x = LaneF32FromF32(Samples[SingleSampleIndex+0].x);
        (*WideSamples).y = LaneF32FromF32(Samples[SingleSampleIndex+0].y);
#elif (LANE_WIDTH == 4)
        
        (*WideSamples).x = LaneF32FromF32(Samples[SingleSampleIndex+0].x,
                                        Samples[SingleSampleIndex+1].x,
                                        Samples[SingleSampleIndex+2].x,
                                        Samples[SingleSampleIndex+3].x
                                        );
        (*WideSamples).y = LaneF32FromF32(Samples[SingleSampleIndex+0].y,
                                        Samples[SingleSampleIndex+1].y,
                                        Samples[SingleSampleIndex+2].y,
                                        Samples[SingleSampleIndex+3].y
                                        );
#elif(LANE_WIDTH == 8)
        
            (*WideSamples).x = LaneF32FromF32(Samples[SingleSampleIndex+0].x,
                                            Samples[SingleSampleIndex+1].x,
                                            Samples[SingleSampleIndex+2].x,
                                            Samples[SingleSampleIndex+3].x,
                                            Samples[SingleSampleIndex+4].x,
                                            Samples[SingleSampleIndex+5].x,
                                            Samples[SingleSampleIndex+6].x,
                                            Samples[SingleSampleIndex+7].x
                                            );
            (*WideSamples).y = LaneF32FromF32(Samples[SingleSampleIndex+0].y,
                                            Samples[SingleSampleIndex+1].y,
                                            Samples[SingleSampleIndex+2].y,
                                            Samples[SingleSampleIndex+3].y,
                                            Samples[SingleSampleIndex+4].y,
                                            Samples[SingleSampleIndex+5].y,
                                            Samples[SingleSampleIndex+6].y,
                                            Samples[SingleSampleIndex+7].y
                                            );
#endif
        WideSamples++;
    }

    free(Samples);
}

// internal void NRooks(lane_v2 *WideSamples, u32 NumSamples, random_series *Series)
// {
//     for(u32 i=0; i<NumSamples; i++)
//     {
//         WideSamples[i].x = (((f32)i + RandomUnilateral(Series)) / (f32)NumSamples);
//         WideSamples[i].y = (((f32)i + RandomUnilateral(Series)) / (f32)NumSamples);
//     }

//     for(u32 i=NumSamples-2; i>0; i--)
//     {
//         u32 Target = (u32)(RandomUnilateral(Series) * (f32)i);
//         f32 tmp = WideSamples[i+1].x;
//         WideSamples[i+1].x = WideSamples[Target].x;
//         WideSamples[Target].x = tmp;
//     }
// }

internal void MultiJitter(lane_v2 *WideSamples, u32 NumSamples, random_series *Series) {
    v2 *Samples = (v2*) malloc(NumSamples * sizeof(v2));

    u32 Width = (u32)SquareRoot(NumSamples);
    f32 CellSize = 1.0f / (f32)(NumSamples);

    for(u32 i=0; i<Width;i++) {
        for(u32 j=0; j<Width;j++) {
            Samples[i * Width + j].x = i * Width * CellSize + j * CellSize + Extract0(RandomUnilateral(Series)) * CellSize;
            Samples[i * Width + j].y = j * Width * CellSize + i * CellSize + Extract0(RandomUnilateral(Series)) * CellSize;
        }
    }

    for(u32 i=0; i<Width;i++) {
        for(u32 j=0; j<Width;j++) {
            u32 k = j + u32( Extract0(RandomUnilateral(Series)) * (Width - j - 1));
            f32 t = Samples[i * Width + j].x;
            Samples[i * Width + j].x = Samples[i * Width + k].x;
            Samples[i * Width + k].x = t;

            
            k = j + u32(Extract0(RandomUnilateral(Series)) * (Width - j - 1));
            t = Samples[j * Width + i].y;
            Samples[j * Width + i].y = Samples[k * Width + i].y;
            Samples[k * Width + i].y = t;
        }
    }
    
    if(Width * Width < NumSamples)
    {
        for(u32 i=Width * Width-1; i<NumSamples; i++)
        {
            Samples[i].x = Extract0(RandomUnilateral(Series));
            Samples[i].y = Extract0(RandomUnilateral(Series));
        }
    }

    for(u32 SingleSampleIndex=0; SingleSampleIndex<NumSamples; SingleSampleIndex+= LANE_WIDTH) {
#if(LANE_WIDTH == 1)
        
        (*WideSamples).x = LaneF32FromF32(Samples[SingleSampleIndex+0].x);
        (*WideSamples).y = LaneF32FromF32(Samples[SingleSampleIndex+0].y);
#elif (LANE_WIDTH == 4)
        
        (*WideSamples).x = LaneF32FromF32(Samples[SingleSampleIndex+0].x,
                                        Samples[SingleSampleIndex+1].x,
                                        Samples[SingleSampleIndex+2].x,
                                        Samples[SingleSampleIndex+3].x
                                        );
        (*WideSamples).y = LaneF32FromF32(Samples[SingleSampleIndex+0].y,
                                        Samples[SingleSampleIndex+1].y,
                                        Samples[SingleSampleIndex+2].y,
                                        Samples[SingleSampleIndex+3].y
                                        );
#elif(LANE_WIDTH == 8)
        
            (*WideSamples).x = LaneF32FromF32(Samples[SingleSampleIndex+0].x,
                                            Samples[SingleSampleIndex+1].x,
                                            Samples[SingleSampleIndex+2].x,
                                            Samples[SingleSampleIndex+3].x,
                                            Samples[SingleSampleIndex+4].x,
                                            Samples[SingleSampleIndex+5].x,
                                            Samples[SingleSampleIndex+6].x,
                                            Samples[SingleSampleIndex+7].x
                                            );
            (*WideSamples).y = LaneF32FromF32(Samples[SingleSampleIndex+0].y,
                                            Samples[SingleSampleIndex+1].y,
                                            Samples[SingleSampleIndex+2].y,
                                            Samples[SingleSampleIndex+3].y,
                                            Samples[SingleSampleIndex+4].y,
                                            Samples[SingleSampleIndex+5].y,
                                            Samples[SingleSampleIndex+6].y,
                                            Samples[SingleSampleIndex+7].y
                                            );
#endif
        WideSamples++;
    }

    free(Samples);
}

void ShuffleArray(lane_v2 *Samples, u32 NumSamples, random_series *Entropy)
{
    u32 n = NumSamples;
    for (u32 i = 0; i < NumSamples; i++)
    {
        u32 j = Extract0(RandomUnilateral(Entropy))  * NumSamples - 1;
        lane_v2 t = Samples[j];
        Samples[j] = Samples[i];
        Samples[i] = t;
    }
}