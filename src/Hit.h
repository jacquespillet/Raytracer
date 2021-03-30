struct hit {
    lane_v3 Position;
    lane_v3 Normal;
    lane_u32 MaterialIndex;
    lane_f32 Distance;
    lane_mat4 Transform;
    lane_mat4 InverseTransform;
};

void HitPlane(plane *Plane, hit *Hit, lane_v3 RayOrigin, lane_v3 RayDirection, lane_f32 Tolerance, lane_f32 MinHitDistance) {
    lane_v3 PlaneN = LaneLaneV3FromLaneV3(Plane->N);
    lane_f32 PlaneD = LaneF32FromF32(Plane->d);

    lane_f32 Denom = Inner(PlaneN, RayDirection);
    lane_u32 DenomMask = ((Denom < -Tolerance) | (Denom > Tolerance));
    if(!MaskIsZeroed(DenomMask)) {
        lane_f32 t = (-PlaneD - Inner(PlaneN, RayOrigin)) / Denom;

        lane_u32 tMask = ((t > MinHitDistance) & (t < Hit->Distance));
        if(!MaskIsZeroed(tMask)) { //If the mask does not contain anything, all the rays in the lane missed. 
            lane_u32 HitMask = DenomMask & tMask;
            
            lane_u32 PlaneMatIndex = LaneU32FromU32(Plane->MatIndex);
            ConditionalAssign(&Hit->Distance, HitMask, t);
            ConditionalAssign(&Hit->MaterialIndex, HitMask, PlaneMatIndex);
            ConditionalAssign(&Hit->Position, HitMask, RayOrigin + t * RayDirection);
            ConditionalAssign(&Hit->Normal, HitMask, PlaneN);

            lane_v3 HitPosition = RayOrigin + t * RayDirection;
            

            ConditionalAssign(&Hit->Transform, HitMask, OrthoBasisFromNormal(Hit->Normal));
            ConditionalAssign(&Hit->InverseTransform, HitMask, Inverse(Hit->Transform));
        }
    }
}

b32 HitSphere(sphere *Sphere,hit *Hit, lane_v3 RayOrigin, lane_v3 RayDirection, lane_f32 Tolerance, lane_f32 MinHitDistance) {

    lane_f32 Spherer = LaneF32FromF32( Sphere->r);
    lane_v3 SphereP = TransformPosition(Sphere->Transform, LaneV3(0,0,0));

    lane_v3 SphereRelativeOrigin =  RayOrigin - SphereP;
    lane_f32 a = Inner(RayDirection, RayDirection);
    lane_f32 b = 2.0f * Inner(RayDirection, SphereRelativeOrigin);
    lane_f32 c = Inner(SphereRelativeOrigin, SphereRelativeOrigin) - Sphere->r * Sphere->r;

    lane_f32 RootTerm = SquareRoot(b*b - 4.0f * a * c);
    
    lane_u32 RootMask = (RootTerm > Tolerance);
    if(!MaskIsZeroed(RootMask)) {
        lane_f32 Denom = 2.0f * a;
        
        lane_f32 tp = (-b + RootTerm) / Denom;
        lane_f32 tn = (-b - RootTerm) / Denom;

        lane_f32 t = tp;

        lane_u32 PickMask = ((tn > MinHitDistance) & (tn < tp));
        ConditionalAssign((lane_f32*)&t, PickMask, tn);
        
        lane_u32 tMask = ((t > MinHitDistance) & (t < Hit->Distance));

        if(!MaskIsZeroed(tMask)) {
            lane_u32 HitMask = RootMask & tMask;
            lane_u32 SphereMatIndex = LaneU32FromU32( Sphere->MatIndex);

            ConditionalAssign(&Hit->Distance, HitMask, t);
            ConditionalAssign(&Hit->MaterialIndex, HitMask, SphereMatIndex);

            ConditionalAssign(&Hit->Position, HitMask, RayOrigin + t * RayDirection);
            ConditionalAssign(&Hit->Normal, HitMask, NOZ(Hit->Position - SphereP));

            ConditionalAssign(&Hit->Transform, HitMask, OrthoBasisFromNormal(Hit->Normal));
            ConditionalAssign(&Hit->InverseTransform, HitMask, Inverse(Hit->Transform));
            
            return true;
        }
    }
    return false;
}

lane_u32 HitAABB(lane_v3 RayOrigin, lane_v3 RayDirection, lane_f32 MinDistance, lane_f32 MaxDistance, aabb AABB) {
    lane_f32 tmp;

    // lane_f32 tmin = LaneF32FromF32(0.001f);
    // lane_f32 tmax = LaneF32FromF32(1000.0f);

    lane_u32 LaneMask = LaneU32FromU32(0xffffffff);

    for(int i=0; i<3; i++) {
        lane_f32 InvD = LaneF32FromF32(1.0f) / LaneV3Component(RayDirection, i);
        lane_f32 t0 = (LaneV3Component(AABB.min, i) - LaneV3Component(RayOrigin, i)) * InvD;
        lane_f32 t1 = (LaneV3Component(AABB.max, i) - LaneV3Component(RayOrigin, i)) * InvD;

        lane_u32 InvDMask = InvD < 0.0f;
        ConditionalAssign(&tmp, InvDMask, t0);
        ConditionalAssign(&t0, InvDMask, t1);
        ConditionalAssign(&t1, InvDMask, tmp);
        
        lane_u32 TMinMaks = t0 > MinDistance; 
        ConditionalAssign(&MinDistance, TMinMaks, t0);
        
        lane_u32 TMaxMaks = t1 < MaxDistance; 
        ConditionalAssign(&MaxDistance, TMaxMaks, t1);
        
        lane_u32 MaxIsLessMask = MaxDistance >= MinDistance;
        
        LaneMask |= MaxIsLessMask;
        if(MaskIsZeroed(MaxIsLessMask)) {
            return MaxIsLessMask;
        }
    }

    return LaneMask; 
}

internal b32 HitBVH(lane_v3 RayOrigin, lane_v3 RayDirection, lane_f32 Tolerance, lane_f32 MinHitDistance, bvh *BVH, hit *Hit, u32 *level) {
    lane_u32 HitBoxMask = HitAABB(RayOrigin, RayDirection,LaneF32FromF32(0.001f),LaneF32FromF32(1000.0f),  BVH->AABB);
    if(!MaskIsZeroed( HitBoxMask)) {
        if(BVH->Spheres == nullptr) { //We're not at a root --> We intersect with the child bvhs
            hit HitLeftRecord = *Hit;
            hit HitRightRecord = HitLeftRecord;
            
            b32 HitLeft = HitBVH(RayOrigin, RayDirection, Tolerance, MinHitDistance,  BVH->Left, &HitLeftRecord, level);
            b32 HitRight = HitBVH(RayOrigin, RayDirection, Tolerance, MinHitDistance,  BVH->Right, &HitRightRecord, level); 

            if(HitLeft && HitRight) {
                lane_u32 MaxIsLessMask = HitLeftRecord.Distance < HitRightRecord.Distance;
                
                //When hit both, we need to take the closest from both
                *Hit = HitRightRecord;
                ConditionalAssign(&Hit->Distance, MaxIsLessMask, HitLeftRecord.Distance);
                ConditionalAssign(&Hit->Normal, MaxIsLessMask, HitLeftRecord.Normal);
                ConditionalAssign(&Hit->MaterialIndex, MaxIsLessMask, HitLeftRecord.MaterialIndex);
                ConditionalAssign(&Hit->Position, MaxIsLessMask, HitLeftRecord.Position);
                
                return true;
            } else  
            if(HitRight) {
                *Hit = HitRightRecord;
                return true;
            } 
            else if(HitLeft) {
                *Hit = HitLeftRecord;
                return true;
            }
            else return false;
        } else {
            b32 Result =HitSphere(BVH->Spheres, Hit, RayOrigin, RayDirection, Tolerance, MinHitDistance); 
            return Result;
        }
    } else {
        return false;
    }
}

void HitVolume(volume *Volume, hit *Hit, lane_v3 RayOrigin, lane_v3 RayDirection, lane_f32 Tolerance, lane_f32 MinHitDistance, random_series *Series) {
    lane_u32 db = RandomUnilateral(Series) < LaneF32FromF32(0.00001f);
    db = LaneU32FromU32(0x00000000);

    hit HitRecord1, HitRecord2;
    lane_u32 HitBoundaryMask = HitAABB(RayOrigin, RayDirection, LaneF32FromF32(-F32Max), LaneF32FromF32(F32Max), Volume->AABB);   
}