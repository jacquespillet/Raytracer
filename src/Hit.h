struct hit {
    lane_v3 Position;
    lane_v3 Normal;
    lane_u32 MaterialIndex;
    lane_f32 Distance;
    lane_mat4 Transform;
    lane_mat4 InverseTransform;

    lane_v3 Wo;
};

void HitPlane(plane *Plane, hit *Hit, lane_v3 RayOrigin, lane_v3 RayDirection, lane_f32 Tolerance, lane_f32 MinHitDistance) {
    lane_v3 PlaneN = LaneV3FromV3(Plane->N);
    lane_f32 PlaneD = LaneF32FromF32(Plane->d);

    lane_f32 Denom = Lane_Inner(PlaneN, RayDirection);
    lane_u32 DenomMask = ((Denom < -Tolerance) | (Denom > Tolerance));
    if(!MaskIsZeroed(DenomMask)) {
        lane_f32 t = (-PlaneD - Lane_Inner(PlaneN, RayOrigin)) / Denom;

        lane_u32 tMask = ((t > MinHitDistance) & (t < Hit->Distance));
        if(!MaskIsZeroed(tMask)) { //If the mask does not contain anything, all the rays in the lane missed. 
            lane_u32 HitMask = DenomMask & tMask;
            
            lane_u32 PlaneMatIndex = LaneU32FromU32(Plane->MatIndex);
            ConditionalAssign(&Hit->Distance, HitMask, t);
            ConditionalAssign(&Hit->MaterialIndex, HitMask, PlaneMatIndex);
            ConditionalAssign(&Hit->Position, HitMask, RayOrigin + t * RayDirection);
            ConditionalAssign(&Hit->Normal, HitMask, PlaneN);

            lane_v3 HitPosition = RayOrigin + t * RayDirection;
            

            ConditionalAssign(&Hit->Transform, HitMask, Lane_OrthoBasisFromNormal(Hit->Normal));
            ConditionalAssign(&Hit->InverseTransform, HitMask, Lane_Inverse(Hit->Transform));
        }
    }
}

b32 HitSphere(shape *Shape,hit *Hit, lane_v3 RayOrigin, lane_v3 RayDirection, lane_f32 Tolerance, lane_f32 MinHitDistance) {
    lane_mat4 ShapeTransform = LaneMat4FromMat4(Shape->Transform);

    lane_f32 Spherer = LaneF32FromF32(Shape->Fields.SphereFields.Sphere.r);
    lane_v3 SphereP = Lane_TransformPosition(ShapeTransform, LaneV3(0,0,0));

    lane_v3 SphereRelativeOrigin =  RayOrigin - SphereP;
    lane_f32 a = Lane_Inner(RayDirection, RayDirection);
    lane_f32 b = 2.0f * Lane_Inner(RayDirection, SphereRelativeOrigin);
    lane_f32 c = Lane_Inner(SphereRelativeOrigin, SphereRelativeOrigin) - Spherer * Spherer;

    lane_f32 RootTerm = Lane_SquareRoot(b*b - 4.0f * a * c);
    
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
            lane_u32 SphereMatIndex = LaneU32FromU32( Shape->MatIndex);

            ConditionalAssign(&Hit->Distance, HitMask, t);
            ConditionalAssign(&Hit->MaterialIndex, HitMask, SphereMatIndex);

            ConditionalAssign(&Hit->Position, HitMask, RayOrigin + t * RayDirection);
            ConditionalAssign(&Hit->Normal, HitMask, Lane_NOZ(Hit->Position - SphereP));

            ConditionalAssign(&Hit->Transform, HitMask, Lane_OrthoBasisFromNormal(Hit->Normal));
            ConditionalAssign(&Hit->InverseTransform, HitMask, Lane_Inverse(Hit->Transform));
            
            return true;
        }
    }
    return false;
}

b32 HitTriangle(shape *Shape,hit *Hit, lane_v3 RayOrigin, lane_v3 RayDirection, lane_f32 Tolerance, lane_f32 MinHitDistance) {
    
    lane_v3 v0 = LaneV3FromV3(Shape->Fields.TriangleFields.Triangle.V1);
    lane_v3 v1= LaneV3FromV3(Shape->Fields.TriangleFields.Triangle.V2);
    lane_v3 v2 = LaneV3FromV3(Shape->Fields.TriangleFields.Triangle.V3);

    lane_v3 n0 = LaneV3FromV3(Shape->Fields.TriangleFields.Triangle.N1);
    lane_v3 n1 = LaneV3FromV3(Shape->Fields.TriangleFields.Triangle.N2);
    lane_v3 n2 = LaneV3FromV3(Shape->Fields.TriangleFields.Triangle.N3);
    
    lane_f32 kEpsilon = LaneF32FromF32(1e-8f); 

    lane_v3 v0v1 = v1 - v0; 
    lane_v3 v0v2 = v2 - v0; 
    lane_v3 pvec = Lane_Cross(RayDirection, v0v2);
    lane_f32 det = Lane_Inner(v0v1, pvec);

    lane_u32 MissLane = LaneU32FromU32(0.0f);
    // ray and triangle are parallel if det is close to 0
    // if (std::abs(det) < kEpsilon) return false; 
    MissLane |= Lane_Abs(det) < kEpsilon;

    lane_f32 invDet = LaneF32FromF32(1.0f) / det; 

    lane_v3 tvec = RayOrigin - v0; 
    lane_f32 u =Lane_Inner(tvec, pvec) * invDet; 
    //if (u < 0 || u > 1) return false; 
    MissLane |= (u < LaneF32FromF32(0.0f) | u > LaneF32FromF32(1.0f));

    lane_v3 qvec = Lane_Cross(tvec, v0v1);   
    lane_f32 v =Lane_Inner(RayDirection, qvec) * invDet; 
    //if (v < 0 | u + v > 1) return false; 
    MissLane |= (v < LaneF32FromF32(0.0f) | u + v > LaneF32FromF32(1.0f));

    lane_f32 t =Lane_Inner(v0v2, qvec) * invDet; 

    MissLane |= (t <= LaneF32FromF32(0.0f));

    lane_v3 HitNormal = (1 - u - v) * n0 + u * n1 + v * n2;

    lane_u32 ReverseNormalMask = Lane_Inner(HitNormal, -RayDirection) <= 0.0f;
    ConditionalAssign(&HitNormal, ReverseNormalMask, HitNormal * LaneF32FromF32(-1.0f));


    lane_u32 IntersectMask = AndNot(MissLane, LaneU32FromU32(0xffffffff));

    ConditionalAssign(&Hit->Distance, IntersectMask, t);
    ConditionalAssign(&Hit->MaterialIndex, IntersectMask, LaneU32FromU32(Shape->MatIndex));

    ConditionalAssign(&Hit->Position, IntersectMask, RayOrigin + t * RayDirection);
    ConditionalAssign(&Hit->Normal, IntersectMask, HitNormal);

    ConditionalAssign(&Hit->Transform, IntersectMask, Lane_OrthoBasisFromNormal(Hit->Normal));
    ConditionalAssign(&Hit->InverseTransform, IntersectMask, Lane_Inverse(Hit->Transform));

    return true;
}

lane_u32 HitAABB(lane_v3 RayOrigin, lane_v3 RayDirection, lane_f32 MinDistance, lane_f32 MaxDistance, aabb AABB) {
    lane_f32 tmp;

    // lane_f32 tmin = LaneF32FromF32(0.001f);
    // lane_f32 tmax = LaneF32FromF32(1000.0f);

    lane_u32 LaneMask = LaneU32FromU32(0xffffffff);

    for(int i=0; i<3; i++) {
        lane_f32 InvD = LaneF32FromF32(1.0f) / LaneV3Component(RayDirection, i);
        lane_f32 t0 = (LaneF32FromF32(V3Component(AABB.min, i)) - LaneV3Component(RayOrigin, i)) * InvD;
        lane_f32 t1 = (LaneF32FromF32(V3Component(AABB.max, i)) - LaneV3Component(RayOrigin, i)) * InvD;

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
        if(BVH->Shapes == nullptr) { //We're not at a root --> We intersect with the child bvhs
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
                ConditionalAssign(&Hit->Transform, MaxIsLessMask, HitLeftRecord.Transform);
                ConditionalAssign(&Hit->InverseTransform, MaxIsLessMask, HitLeftRecord.InverseTransform);
                
                
                
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
            if(BVH->Shapes->Type == shape_type::sphereType)
            {
                b32 Result =HitSphere(BVH->Shapes, Hit, RayOrigin, RayDirection, Tolerance, MinHitDistance); 
                return Result;
            } else if(BVH->Shapes->Type == shape_type::triangleType)
            {
                b32 Result =HitTriangle(BVH->Shapes, Hit, RayOrigin, RayDirection, Tolerance, MinHitDistance); 
                return Result;
            }
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