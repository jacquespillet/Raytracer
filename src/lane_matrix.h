//Row major.
//In the memory, first come the individual Rows

struct lane_mat2
{
    lane_f32 Elements[4];
};

struct lane_mat3
{
    lane_f32 Elements[9];
};

struct lane_mat4
{
    lane_f32 Elements[16];
};

internal void ConditionalAssign(lane_mat4 *Dest, lane_u32 Mask, lane_mat4 Source) {
    ConditionalAssign(&Dest->Elements[0], Mask, Source.Elements[0]);
    ConditionalAssign(&Dest->Elements[1], Mask, Source.Elements[1]);
    ConditionalAssign(&Dest->Elements[2], Mask, Source.Elements[2]);
    ConditionalAssign(&Dest->Elements[3], Mask, Source.Elements[3]);

    ConditionalAssign(&Dest->Elements[4], Mask, Source.Elements[4]);
    ConditionalAssign(&Dest->Elements[5], Mask, Source.Elements[5]);
    ConditionalAssign(&Dest->Elements[6], Mask, Source.Elements[6]);
    ConditionalAssign(&Dest->Elements[7], Mask, Source.Elements[7]);
    
    ConditionalAssign(&Dest->Elements[8], Mask, Source.Elements[8]);
    ConditionalAssign(&Dest->Elements[9], Mask, Source.Elements[9]);
    ConditionalAssign(&Dest->Elements[10], Mask, Source.Elements[10]);
    ConditionalAssign(&Dest->Elements[11], Mask, Source.Elements[11]);
    
    ConditionalAssign(&Dest->Elements[12], Mask, Source.Elements[12]);
    ConditionalAssign(&Dest->Elements[13], Mask, Source.Elements[13]);
    ConditionalAssign(&Dest->Elements[14], Mask, Source.Elements[14]);
    ConditionalAssign(&Dest->Elements[15], Mask, Source.Elements[15]);
}

lane_mat2 operator*(lane_mat2 &Matrix, lane_f32 &Value)
{
    lane_mat2 Result = {};
    for(u8 Index=0; Index<4; Index++)
    {
        Result.Elements[Index] = Value * Matrix.Elements[Index];
    }
    return Result;
}

lane_mat2 operator*(lane_f32 &Value, lane_mat2 &Matrix)
{
    lane_mat2 Result = Matrix * Value;
    return Result;
}

lane_mat3 operator*(lane_mat3 &Matrix, lane_f32 &Value)
{
    lane_mat3 Result = {};
    for(u8 Index=0; Index<9; Index++)
    {
        Result.Elements[Index] = Value * Matrix.Elements[Index];
    }
    return Result;
}

lane_mat3 operator*(lane_f32 &Value, lane_mat3 &Matrix)
{
    lane_mat3 Result = Matrix * Value;
    return Result;
}

lane_mat4 operator*(lane_mat4 Matrix, lane_f32 Value)
{
    lane_mat4 Result = {};
    for(u8 Index=0; Index<16; Index++)
    {
        Result.Elements[Index] = Value * Matrix.Elements[Index];
    }
    return Result;
}

lane_mat4 operator*(lane_f32 Value, lane_mat4 Matrix)
{
    lane_mat4 Result = Matrix * Value;
    return Result;
}


lane_f32 Lane_GetMatrixElement(lane_mat2 Matrix, u32 Row, u32 Column)
{
    u32 Index = 2 * Row + Column;
    
    return Matrix.Elements[Index];
}

lane_f32 Lane_GetMatrixElement(lane_mat3 Matrix, u32 Row, u32 Column)
{
    u32 Index = 3 * Row + Column;
    return Matrix.Elements[Index];
}

lane_f32 Lane_GetMatrixElement(lane_mat4 Matrix, u32 Row, u32 Column)
{
    u32 Index = 4 * Row + Column;
    return Matrix.Elements[Index];
}

void Lane_SetMatrixElement(lane_mat2 *Matrix, u32 Row, u32 Column, lane_f32 Value)
{
    u32 Index = 2 * Row + Column;
    Matrix->Elements[Index] = Value;
}

void Lane_SetMatrixElement(lane_mat3 *Matrix, u32 Row, u32 Column, lane_f32 Value)
{
    u32 Index = 3 * Row + Column;
    Matrix->Elements[Index] = Value;
}

void Lane_SetMatrixElement(lane_mat4 *Matrix, u32 Row, u32 Column, lane_f32 Value)
{
    s32 Index = 4 * Row + Column;
    Matrix->Elements[Index] = Value;
}

lane_mat2 LaneMat2F(lane_f32 Value)
{
    lane_mat2 Result ={};
    Lane_SetMatrixElement(&Result, 0, 0, Value);
    Lane_SetMatrixElement(&Result, 1, 1, Value);
    return Result;
}

lane_mat3 LaneMat3F(lane_f32 Value)
{
    lane_mat3 Result ={};
    Lane_SetMatrixElement(&Result, 0, 0, Value);
    Lane_SetMatrixElement(&Result, 1, 1, Value);
    Lane_SetMatrixElement(&Result, 2, 2, Value);
    return Result;
}

lane_mat4 LaneMat4F(lane_f32 Value)
{
    lane_mat4 Result ={};
    Lane_SetMatrixElement(&Result, 0, 0, Value);
    Lane_SetMatrixElement(&Result, 1, 1, Value);
    Lane_SetMatrixElement(&Result, 2, 2, Value);
    Lane_SetMatrixElement(&Result, 3, 3, Value);
    return Result;
}

lane_mat2 Lane_SubMatrix(lane_mat3 Matrix, u8 Row, u8 Column)
{
    lane_mat2 Result = {};
    u8 RowIndices[2] = {};
    u8 ColumnIndices[2] = {};
    u8 runningRow=0;
    u8 runningColumn=0;
    for(u8 Index=0; Index<3; Index++)
    {
        if(Index != Row)
        {
            RowIndices[runningRow++] = Index;
        }
        if(Index != Column)
        {
            ColumnIndices[runningColumn++] = Index;
        }
    }

    for(u8 RowIndex=0; RowIndex < 2; RowIndex++)
    {
        for(u8 ColumnIndex=0; ColumnIndex < 2; ColumnIndex++)
        {
            lane_f32 Value = Lane_GetMatrixElement(Matrix, RowIndices[RowIndex], ColumnIndices[ColumnIndex]);
            Lane_SetMatrixElement(&Result, RowIndex, ColumnIndex, Value);
        }     
    }
    return Result;
}

lane_mat3 Lane_SubMatrix(lane_mat4 Matrix, u8 Row, u8 Column)
{
    lane_mat3 Result = {};
    u8 RowIndices[3] = {};
    u8 ColumnIndices[3] = {};
    u8 runningRow=0;
    u8 runningColumn=0;
    for(u8 Index=0; Index<4; Index++)
    {
        if(Index != Row)
        {
            RowIndices[runningRow++] = Index;
        }
        if(Index != Column)
        {
            ColumnIndices[runningColumn++] = Index;
        }
    }

    for(u8 RowIndex=0; RowIndex < 3; RowIndex++)
    {
        for(u8 ColumnIndex=0; ColumnIndex < 3; ColumnIndex++)
        {
            lane_f32 Value = Lane_GetMatrixElement(Matrix, RowIndices[RowIndex], ColumnIndices[ColumnIndex]);
            Lane_SetMatrixElement(&Result, RowIndex, ColumnIndex, Value);
        }     
    }
    return Result;
}

lane_mat2 operator*(lane_mat2 &A, lane_mat2 &B)

{
    lane_mat2 Result = {};

    for(u8 RowIndex=0; RowIndex<2; RowIndex++)
    {
        for(u8 ColumnIndex=0; ColumnIndex<2; ColumnIndex++)
        {
            lane_f32 Value =   Lane_GetMatrixElement(A, RowIndex, 0) * Lane_GetMatrixElement(B, 0, ColumnIndex) 
                        + Lane_GetMatrixElement(A, RowIndex, 1) * Lane_GetMatrixElement(B, 1, ColumnIndex);

            Lane_SetMatrixElement(&Result, RowIndex, ColumnIndex, Value);
        }
    }
    return Result;
}

lane_mat3 operator*(lane_mat3 &A, lane_mat3 &B)

{
    lane_mat3 Result = {};

    for(u8 RowIndex=0; RowIndex<3; RowIndex++)
    {
        for(u8 ColumnIndex=0; ColumnIndex<3; ColumnIndex++)
        {
            lane_f32 Value =   Lane_GetMatrixElement(A, RowIndex, 0) * Lane_GetMatrixElement(B, 0, ColumnIndex) 
                        + Lane_GetMatrixElement(A, RowIndex, 1) * Lane_GetMatrixElement(B, 1, ColumnIndex) 
                        + Lane_GetMatrixElement(A, RowIndex, 2) * Lane_GetMatrixElement(B, 2, ColumnIndex);
            Lane_SetMatrixElement(&Result, RowIndex, ColumnIndex, Value);
        }
    }

    return Result;
}

lane_mat4 operator*(lane_mat4 &MatrixA, lane_mat4 &MatrixB)

{
    lane_mat4 Result = {};

    for(u8 RowIndex=0; RowIndex<4; RowIndex++)
    {
        for(u8 ColumnIndex=0; ColumnIndex<4; ColumnIndex++)
        {
            lane_f32 A = Lane_GetMatrixElement(MatrixA, RowIndex, 0) * Lane_GetMatrixElement(MatrixB, 0, ColumnIndex);
            lane_f32 B = Lane_GetMatrixElement(MatrixA, RowIndex, 1) * Lane_GetMatrixElement(MatrixB, 1, ColumnIndex);
            lane_f32 C = Lane_GetMatrixElement(MatrixA, RowIndex, 2) * Lane_GetMatrixElement(MatrixB, 2, ColumnIndex);
            lane_f32 D = Lane_GetMatrixElement(MatrixA, RowIndex, 3) * Lane_GetMatrixElement(MatrixB, 3, ColumnIndex);
            lane_f32 Value = A + B + C + D;
            Lane_SetMatrixElement(&Result, RowIndex, ColumnIndex, Value);
        }
    }

    return Result;
}


lane_v2 operator*(lane_mat2 &Matrix, lane_v2 &Vector)
{
    lane_v2 Result = {};

    Result.x = Lane_GetMatrixElement(Matrix, 0, 0) * Vector.x + Lane_GetMatrixElement(Matrix, 0, 1) * Vector.y;
    Result.y = Lane_GetMatrixElement(Matrix, 1, 0) * Vector.x + Lane_GetMatrixElement(Matrix, 1, 1) * Vector.y;

    return Result;
}

lane_v3 operator*(lane_mat3 &Matrix, lane_v3 &Vector)
{
    lane_v3 Result = {};

    Result.x = Lane_GetMatrixElement(Matrix, 0, 0) * Vector.x + Lane_GetMatrixElement(Matrix, 0, 1) * Vector.y + Lane_GetMatrixElement(Matrix, 0, 2) * Vector.z;
    Result.y = Lane_GetMatrixElement(Matrix, 1, 0) * Vector.x + Lane_GetMatrixElement(Matrix, 1, 1) * Vector.y + Lane_GetMatrixElement(Matrix, 1, 2) * Vector.z;
    Result.z = Lane_GetMatrixElement(Matrix, 2, 0) * Vector.x + Lane_GetMatrixElement(Matrix, 2, 1) * Vector.y + Lane_GetMatrixElement(Matrix, 2, 2) * Vector.z;

    return Result;
}

lane_v4 operator*(lane_mat4 &Matrix, lane_v4 &Vector)
{
    lane_v4 Result = {};

    Result.x = Lane_GetMatrixElement(Matrix, 0, 0) * Vector.x + Lane_GetMatrixElement(Matrix, 0, 1) * Vector.y + Lane_GetMatrixElement(Matrix, 0, 2) * Vector.z + Lane_GetMatrixElement(Matrix, 0, 3) * Vector.w;
    Result.y = Lane_GetMatrixElement(Matrix, 1, 0) * Vector.x + Lane_GetMatrixElement(Matrix, 1, 1) * Vector.y + Lane_GetMatrixElement(Matrix, 1, 2) * Vector.z + Lane_GetMatrixElement(Matrix, 1, 3) * Vector.w;
    Result.z = Lane_GetMatrixElement(Matrix, 2, 0) * Vector.x + Lane_GetMatrixElement(Matrix, 2, 1) * Vector.y + Lane_GetMatrixElement(Matrix, 2, 2) * Vector.z + Lane_GetMatrixElement(Matrix, 2, 3) * Vector.w;
    Result.w = Lane_GetMatrixElement(Matrix, 3, 0) * Vector.x + Lane_GetMatrixElement(Matrix, 3, 1) * Vector.y + Lane_GetMatrixElement(Matrix, 3, 2) * Vector.z + Lane_GetMatrixElement(Matrix, 3, 3) * Vector.w;

    return Result;
}

lane_f32 Lane_Determinant(lane_mat2 Input)
{
    lane_f32 Result = LaneF32FromF32(0);

    lane_f32 a = Lane_GetMatrixElement(Input, 0, 0);
    lane_f32 b = Lane_GetMatrixElement(Input, 0, 1);
    lane_f32 c = Lane_GetMatrixElement(Input, 1, 0);
    lane_f32 d = Lane_GetMatrixElement(Input, 1, 1);

    Result = a * d - b * c;
    return Result;
}

lane_f32 Lane_Determinant(lane_mat3 Input)
{
   //Find cofactors
    lane_f32 DetA = Lane_Determinant(Lane_SubMatrix(Input, 0, 0));
    lane_f32 DetB = Lane_Determinant(Lane_SubMatrix(Input, 0, 1));
    lane_f32 DetC = Lane_Determinant(Lane_SubMatrix(Input, 0, 2));

    //Find determinant
    lane_f32 Result = Lane_GetMatrixElement(Input, 0, 0) * DetA - Lane_GetMatrixElement(Input, 0, 1) * DetB + Lane_GetMatrixElement(Input, 0, 2) * DetC;
    return Result;
}

lane_f32 Lane_Determinant(lane_mat4 Input)
{
   //Find cofactors
    lane_f32 DetA = Lane_Determinant(Lane_SubMatrix(Input, 0, 0));
    lane_f32 DetB = Lane_Determinant(Lane_SubMatrix(Input, 0, 1));
    lane_f32 DetC = Lane_Determinant(Lane_SubMatrix(Input, 0, 2));
    lane_f32 DetD = Lane_Determinant(Lane_SubMatrix(Input, 0, 3));

    //Find determinant
    lane_f32 Result = Lane_GetMatrixElement(Input, 0, 0) * DetA - Lane_GetMatrixElement(Input, 0, 1) * DetB + Lane_GetMatrixElement(Input, 0, 2) * DetC - Lane_GetMatrixElement(Input, 0, 3) * DetD;
    return Result;
}

lane_mat2 Lane_Inverse(lane_mat2 Input)
{
    lane_mat2 Result = {};

    lane_f32 OneOverDet = 1.0f / Lane_Determinant(Input);
    
    
    lane_mat2 SwappedMatrix = {};
    
    lane_f32 a = Lane_GetMatrixElement(Input, 0, 0);
    lane_f32 b = Lane_GetMatrixElement(Input, 0, 1);
    lane_f32 c = Lane_GetMatrixElement(Input, 1, 0);
    lane_f32 d = Lane_GetMatrixElement(Input, 1, 1);

    Lane_SetMatrixElement(&SwappedMatrix, 0, 0,  d);
    Lane_SetMatrixElement(&SwappedMatrix, 0, 1, -b);
    Lane_SetMatrixElement(&SwappedMatrix, 1, 0, -c);
    Lane_SetMatrixElement(&SwappedMatrix, 1, 1,  a);

    Result = OneOverDet * SwappedMatrix;
    return Result;
}

lane_mat3 Lane_Inverse(lane_mat3 Input)
{
    //Find cofactors
    lane_f32 DetA = Lane_Determinant(Lane_SubMatrix(Input, 0, 0));
    lane_f32 DetB = Lane_Determinant(Lane_SubMatrix(Input, 0, 1));
    lane_f32 DetC = Lane_Determinant(Lane_SubMatrix(Input, 0, 2));

    //Find determinant
    lane_f32 MatrixDeterminant = Lane_GetMatrixElement(Input, 0, 0) * DetA - Lane_GetMatrixElement(Input, 0, 1) * DetB + Lane_GetMatrixElement(Input, 0, 2) * DetC;
    lane_f32 OneOverDeterminant = 1.0f / MatrixDeterminant;
    //Build cofactor matrix
    lane_mat3 CofactorMatrix =  {};

	lane_mat3 CofactorSigns = 
	{
		{
			 1.0f, -1.0f, 1.0f, 
		    -1.0f,  1.0f,-1.0f, 
			 1.0f, -1.0f, 1.0f, 
		}
	};

    for(u8 RowIndex=0; RowIndex < 3; RowIndex++)
    {
        for(u8 ColumnIndex=0; ColumnIndex < 3; ColumnIndex++)
        {
            lane_f32 Cofactor = LaneF32FromF32(0);
            if(RowIndex==0 && ColumnIndex==0) Cofactor = DetA;
            if(RowIndex==0 && ColumnIndex==1) Cofactor = DetB;
            if(RowIndex==0 && ColumnIndex==2) Cofactor = DetC;
            else Cofactor = Lane_Determinant(Lane_SubMatrix(Input, RowIndex, ColumnIndex));
            
            
            lane_f32 Sign = Lane_GetMatrixElement(CofactorSigns, RowIndex, ColumnIndex);
            Cofactor = Cofactor * Sign;

            Cofactor = Cofactor * OneOverDeterminant;

            //Transposing in place!
            Lane_SetMatrixElement(&CofactorMatrix, ColumnIndex, RowIndex, Cofactor);
        }        
    }

    return CofactorMatrix;
}


lane_mat4 Lane_Inverse(lane_mat4 Input)
{
//First version is optimized code
//Second version is easy to read code
#if 1
    lane_f32 a00 = Lane_GetMatrixElement(Input, 0, 0);
    lane_f32 a01 = Lane_GetMatrixElement(Input, 0, 1);
    lane_f32 a02 = Lane_GetMatrixElement(Input, 0, 2);
    lane_f32 a03 = Lane_GetMatrixElement(Input, 0, 3);
    
    lane_f32 a10 = Lane_GetMatrixElement(Input, 1, 0);
    lane_f32 a11 = Lane_GetMatrixElement(Input, 1, 1);
    lane_f32 a12 = Lane_GetMatrixElement(Input, 1, 2);
    lane_f32 a13 = Lane_GetMatrixElement(Input, 1, 3);
    
    lane_f32 a20 = Lane_GetMatrixElement(Input, 2, 0);
    lane_f32 a21 = Lane_GetMatrixElement(Input, 2, 1);
    lane_f32 a22 = Lane_GetMatrixElement(Input, 2, 2);
    lane_f32 a23 = Lane_GetMatrixElement(Input, 2, 3);
    
    lane_f32 a30 = Lane_GetMatrixElement(Input, 3, 0);
    lane_f32 a31 = Lane_GetMatrixElement(Input, 3, 1);
    lane_f32 a32 = Lane_GetMatrixElement(Input, 3, 2);
    lane_f32 a33 = Lane_GetMatrixElement(Input, 3, 3);

    lane_f32 A2323 = a22 * a33 - a23 * a32 ;
    lane_f32 A1323 = a21 * a33 - a23 * a31 ;
    lane_f32 A1223 = a21 * a32 - a22 * a31 ;
    lane_f32 A0323 = a20 * a33 - a23 * a30 ;
    lane_f32 A0223 = a20 * a32 - a22 * a30 ;
    lane_f32 A0123 = a20 * a31 - a21 * a30 ;
    lane_f32 A2313 = a12 * a33 - a13 * a32 ;
    lane_f32 A1313 = a11 * a33 - a13 * a31 ;
    lane_f32 A1213 = a11 * a32 - a12 * a31 ;
    lane_f32 A2312 = a12 * a23 - a13 * a22 ;
    lane_f32 A1312 = a11 * a23 - a13 * a21 ;
    lane_f32 A1212 = a11 * a22 - a12 * a21 ;
    lane_f32 A0313 = a10 * a33 - a13 * a30 ;
    lane_f32 A0213 = a10 * a32 - a12 * a30 ;
    lane_f32 A0312 = a10 * a23 - a13 * a20 ;
    lane_f32 A0212 = a10 * a22 - a12 * a20 ;
    lane_f32 A0113 = a10 * a31 - a11 * a30 ;
    lane_f32 A0112 = a10 * a21 - a11 * a20 ;

    lane_f32 det = a00 * ( a11 * A2323 - a12 * A1323 + a13 * A1223 ) 
        - a01 * ( a10 * A2323 - a12 * A0323 + a13 * A0223 ) 
        + a02 * ( a10 * A1323 - a11 * A0323 + a13 * A0123 ) 
        - a03 * ( a10 * A1223 - a11 * A0223 + a12 * A0123 ) ;
    det = 1 / det;

   lane_f32 b00 = det *   ( a11 * A2323 - a12 * A1323 + a13 * A1223 );
   lane_f32 b01 = det * - ( a01 * A2323 - a02 * A1323 + a03 * A1223 );
   lane_f32 b02 = det *   ( a01 * A2313 - a02 * A1313 + a03 * A1213 );
   lane_f32 b03 = det * - ( a01 * A2312 - a02 * A1312 + a03 * A1212 );
   lane_f32 b10 = det * - ( a10 * A2323 - a12 * A0323 + a13 * A0223 );
   lane_f32 b11 = det *   ( a00 * A2323 - a02 * A0323 + a03 * A0223 );
   lane_f32 b12 = det * - ( a00 * A2313 - a02 * A0313 + a03 * A0213 );
   lane_f32 b13 = det *   ( a00 * A2312 - a02 * A0312 + a03 * A0212 );
   lane_f32 b20 = det *   ( a10 * A1323 - a11 * A0323 + a13 * A0123 );
   lane_f32 b21 = det * - ( a00 * A1323 - a01 * A0323 + a03 * A0123 );
   lane_f32 b22 = det *   ( a00 * A1313 - a01 * A0313 + a03 * A0113 );
   lane_f32 b23 = det * - ( a00 * A1312 - a01 * A0312 + a03 * A0112 );
   lane_f32 b30 = det * - ( a10 * A1223 - a11 * A0223 + a12 * A0123 );
   lane_f32 b31 = det *   ( a00 * A1223 - a01 * A0223 + a02 * A0123 );
   lane_f32 b32 = det * - ( a00 * A1213 - a01 * A0213 + a02 * A0113 );
   lane_f32 b33 = det *   ( a00 * A1212 - a01 * A0212 + a02 * A0112 );

    lane_mat4 Result = {};
    Lane_SetMatrixElement(&Result, 0, 0, b00);
    Lane_SetMatrixElement(&Result, 0, 1, b01);
    Lane_SetMatrixElement(&Result, 0, 2, b02);
    Lane_SetMatrixElement(&Result, 0, 3, b03);
    
    Lane_SetMatrixElement(&Result, 1, 0, b10);
    Lane_SetMatrixElement(&Result, 1, 1, b11);
    Lane_SetMatrixElement(&Result, 1, 2, b12);
    Lane_SetMatrixElement(&Result, 1, 3, b13);
    
    Lane_SetMatrixElement(&Result, 2, 0, b20);
    Lane_SetMatrixElement(&Result, 2, 1, b21);
    Lane_SetMatrixElement(&Result, 2, 2, b22);
    Lane_SetMatrixElement(&Result, 2, 3, b23);
    
    Lane_SetMatrixElement(&Result, 3, 0, b30);
    Lane_SetMatrixElement(&Result, 3, 1, b31);
    Lane_SetMatrixElement(&Result, 3, 2, b32);
    Lane_SetMatrixElement(&Result, 3, 3, b33);

    return Result;
    
#else
    //Find cofactors
    lane_f32 DetA = Determinant(SubMatrix(Input, 0, 0));
    lane_f32 DetB = Determinant(SubMatrix(Input, 0, 1));
    lane_f32 DetC = Determinant(SubMatrix(Input, 0, 2));
    lane_f32 DetD = Determinant(SubMatrix(Input, 0, 3));

    //Find determinant
    lane_f32 MatrixDeterminant = Lane_GetMatrixElement(Input, 0, 0) * DetA - Lane_GetMatrixElement(Input, 0, 1) * DetB + Lane_GetMatrixElement(Input, 0, 2) * DetC - Lane_GetMatrixElement(Input, 0, 3) * DetD;
    lane_f32 OneOverDeterminant = 1.0f / MatrixDeterminant;
    //Build cofactor matrix
    lane_mat4 CofactorMatrix =  {};

	lane_mat4 CofactorSigns = 
	{
		{
			 1.0f, -1.0f, 1.0f, -1.0f,
		    -1.0f,  1.0f,-1.0f,  1.0f,
			 1.0f, -1.0f, 1.0f, -1.0f,
			-1.0f,  1.0f,-1.0f,  1.0f
		}
	};

    for(u8 RowIndex=0; RowIndex < 4; RowIndex++)
    {
        for(u8 ColumnIndex=0; ColumnIndex < 4; ColumnIndex++)
        {
            lane_f32 Cofactor = 0;
            if(RowIndex==0 && ColumnIndex==0) Cofactor = DetA;
            if(RowIndex==0 && ColumnIndex==1) Cofactor = DetB;
            if(RowIndex==0 && ColumnIndex==2) Cofactor = DetC;
            if(RowIndex==0 && ColumnIndex==3) Cofactor = DetD;
            else Cofactor = Determinant(SubMatrix(Input, RowIndex, ColumnIndex));
            
            
            lane_f32 Sign = Lane_GetMatrixElement(CofactorSigns, RowIndex, ColumnIndex);
            Cofactor *= Sign;

            Cofactor *= OneOverDeterminant;

            //Transposing in place!
            Lane_SetMatrixElement(&CofactorMatrix, ColumnIndex, RowIndex, Cofactor);
        }        
    }
    
    return CofactorMatrix;
#endif
}


lane_mat3 Lane_Translate(lane_v2 Translation) 
{
    lane_mat3 Result = LaneMat3F(LaneF32FromF32(1.0f));
    Lane_SetMatrixElement(&Result, 0, 2, Translation.x);
    Lane_SetMatrixElement(&Result, 1, 2, Translation.y);
    return Result;
}

lane_mat3 Lane_Scale(lane_v2 Size) 
{
    lane_mat3 Result = LaneMat3F(LaneF32FromF32(1.0f));
    Lane_SetMatrixElement(&Result, 0, 0, Size.x);
    Lane_SetMatrixElement(&Result, 1, 1, Size.y);
    return Result;
}

lane_mat3 Lane_Rotate(f32 Angle) 
{
    lane_mat3 Result = LaneMat3F(LaneF32FromF32(1.0f));
    // Lane_SetMatrixElement(&Result, 0, 0, Cosine(Angle));
    // Lane_SetMatrixElement(&Result, 0, 1, Sine(Angle));
    // Lane_SetMatrixElement(&Result, 1, 0, -Sine(Angle));
    // Lane_SetMatrixElement(&Result, 1, 1, Cosine(Angle));
    return Result;
}


lane_mat4 Lane_LookAt(lane_v3 CameraPosition, lane_v3 Center, lane_v3 UpVector)
{
    lane_v3 Z = Lane_NOZ(Center - CameraPosition);
    UpVector = LaneV3(0,1,0);
    
    lane_u32 IsTooCloseMask = Lane_Abs(Lane_Inner(UpVector, Z)) > LaneF32FromF32(0.99f); 
    ConditionalAssign(&UpVector, IsTooCloseMask, LaneV3(0,0,1));
    
    IsTooCloseMask = Lane_Abs(Lane_Inner(UpVector, Z)) > LaneF32FromF32(0.99f); 
    ConditionalAssign(&UpVector, IsTooCloseMask, LaneV3(1,0,1));


    lane_v3 X = Lane_NOZ(Lane_Cross(UpVector, Z));
    lane_v3 Y = Lane_NOZ(Lane_Cross(Z, X));

    lane_mat4 Result = {};

    Lane_SetMatrixElement(&Result, 0, 0, X.x);
    Lane_SetMatrixElement(&Result, 1, 0, X.y);
    Lane_SetMatrixElement(&Result, 2, 0, X.z);
    Lane_SetMatrixElement(&Result, 3, 0, LaneF32FromF32(0));

    Lane_SetMatrixElement(&Result, 0, 1, Y.x);
    Lane_SetMatrixElement(&Result, 1, 1, Y.y);
    Lane_SetMatrixElement(&Result, 2, 1, Y.z);
    Lane_SetMatrixElement(&Result, 3, 1, LaneF32FromF32(0));
    
    Lane_SetMatrixElement(&Result, 0, 2, Z.x);
    Lane_SetMatrixElement(&Result, 1, 2, Z.y);
    Lane_SetMatrixElement(&Result, 2, 2, Z.z);
    Lane_SetMatrixElement(&Result, 3, 2, LaneF32FromF32(0));
    
    Lane_SetMatrixElement(&Result, 0, 3, CameraPosition.x);
    Lane_SetMatrixElement(&Result, 1, 3, CameraPosition.y);
    Lane_SetMatrixElement(&Result, 2, 3, CameraPosition.z);
    Lane_SetMatrixElement(&Result, 3, 3, LaneF32FromF32(1));
    
    // lane_mat4 Result = {
    //     X.x, Y.x, Z.x, CameraPosition.x,
    //     X.y, Y.y, Z.y, CameraPosition.y,
    //     X.z, Y.z, Z.z, CameraPosition.z,
    //     0  , 0  , 0  , 1
    // };

    return Result;
}

lane_v3 Lane_TransformPosition(lane_mat4 Matrix, lane_v3 Vector)
{
    lane_v3 Result = {};
    lane_v4 ResultV4 = Matrix * LaneV4(Vector.x, Vector.y, Vector.z, LaneF32FromF32(1.0f));

    Result.x = ResultV4.x;
    Result.y = ResultV4.y;
    Result.z = ResultV4.z;

    return Result;
}


lane_v3 Lane_TransformDirection(lane_mat4 Matrix, lane_v3 Vector)
{
    lane_v3 Result = {};
    lane_v4 ResultV4 = Matrix * LaneV4(Vector.x, Vector.y, Vector.z, LaneF32FromF32(0.0f));

    Result.x = ResultV4.x;
    Result.y = ResultV4.y;
    Result.z = ResultV4.z;

    return Result;
}

lane_mat4 Lane_Identity()
{
    return LaneMat4F(LaneF32FromF32(1.0f));
}

lane_mat4 Lane_OrthoBasisFromNormal(lane_v3 Normal)
{
    lane_v3 UpVector = LaneV3(0,1,0);
    
    // if(Abs(Inner(UpVector, Normal)) > 0.99)
    // {
    //     UpVector = LaneV3(0,0,1);
    // }
    
    // if(Abs(Inner(UpVector, Normal)) > 0.99)
    // {
    //     UpVector = LaneV3(1,0,0);
    // }

    lane_mat4 Result = Lane_LookAt(LaneV3(0.0f, 0.0f, 0.0f), Normal, UpVector);

    return Result;
}

lane_mat4 Lane_Translate(lane_mat4 Matrix, lane_v3 Position)
{
    lane_mat4 Result = Matrix;
    Lane_SetMatrixElement(&Result, 0, 3, Position.x);
    Lane_SetMatrixElement(&Result, 1, 3, Position.y);
    Lane_SetMatrixElement(&Result, 2, 3, Position.z);
    return Result;
}
