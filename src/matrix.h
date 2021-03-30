//Row major.
//In the memory, first come the individual Rows

struct mat2
{
    f32 Elements[4];
};

struct mat3
{
    f32 Elements[9];
};

struct mat4
{
    f32 Elements[16];
};

mat2 operator*(mat2 &Matrix, f32 &Value)
{
    mat2 Result = {};
    for(u8 Index=0; Index<4; Index++)
    {
        Result.Elements[Index] = Value * Matrix.Elements[Index];
    }
    return Result;
}

mat2 operator*(f32 &Value, mat2 &Matrix)
{
    mat2 Result = Matrix * Value;
    return Result;
}

mat3 operator*(mat3 &Matrix, f32 &Value)
{
    mat3 Result = {};
    for(u8 Index=0; Index<9; Index++)
    {
        Result.Elements[Index] = Value * Matrix.Elements[Index];
    }
    return Result;
}

mat3 operator*(f32 &Value, mat3 &Matrix)
{
    mat3 Result = Matrix * Value;
    return Result;
}

mat4 operator*(mat4 Matrix, f32 Value)
{
    mat4 Result = {};
    for(u8 Index=0; Index<16; Index++)
    {
        Result.Elements[Index] = Value * Matrix.Elements[Index];
    }
    return Result;
}

mat4 operator*(f32 Value, mat4 Matrix)
{
    mat4 Result = Matrix * Value;
    return Result;
}


f32 GetMatrixElement(mat2 Matrix, u32 Row, u32 Column)
{
    u32 Index = 2 * Row + Column;
    
    return Matrix.Elements[Index];
}

f32 GetMatrixElement(mat3 Matrix, u32 Row, u32 Column)
{
    u32 Index = 3 * Row + Column;
    return Matrix.Elements[Index];
}

f32 GetMatrixElement(mat4 Matrix, u32 Row, u32 Column)
{
    u32 Index = 4 * Row + Column;
    return Matrix.Elements[Index];
}

void SetMatrixElement(mat2 *Matrix, u32 Row, u32 Column, f32 Value)
{
    u32 Index = 2 * Row + Column;
    Matrix->Elements[Index] = Value;
}

void SetMatrixElement(mat3 *Matrix, u32 Row, u32 Column, f32 Value)
{
    u32 Index = 3 * Row + Column;
    Matrix->Elements[Index] = Value;
}

void SetMatrixElement(mat4 *Matrix, u32 Row, u32 Column, f32 Value)
{
    s32 Index = 4 * Row + Column;
    Matrix->Elements[Index] = Value;
}

mat2 Mat2F(f32 Value)
{
    mat2 Result ={};
    SetMatrixElement(&Result, 0, 0, Value);
    SetMatrixElement(&Result, 1, 1, Value);
    return Result;
}

mat3 Mat3F(f32 Value)
{
    mat3 Result ={};
    SetMatrixElement(&Result, 0, 0, Value);
    SetMatrixElement(&Result, 1, 1, Value);
    SetMatrixElement(&Result, 2, 2, Value);
    return Result;
}

mat4 Mat4F(f32 Value)
{
    mat4 Result ={};
    SetMatrixElement(&Result, 0, 0, Value);
    SetMatrixElement(&Result, 1, 1, Value);
    SetMatrixElement(&Result, 2, 2, Value);
    SetMatrixElement(&Result, 3, 3, Value);
    return Result;
}

mat2 SubMatrix(mat3 Matrix, u8 Row, u8 Column)
{
    mat2 Result = {};
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
            f32 Value = GetMatrixElement(Matrix, RowIndices[RowIndex], ColumnIndices[ColumnIndex]);
            SetMatrixElement(&Result, RowIndex, ColumnIndex, Value);
        }     
    }
    return Result;
}

mat3 SubMatrix(mat4 Matrix, u8 Row, u8 Column)
{
    mat3 Result = {};
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
            f32 Value = GetMatrixElement(Matrix, RowIndices[RowIndex], ColumnIndices[ColumnIndex]);
            SetMatrixElement(&Result, RowIndex, ColumnIndex, Value);
        }     
    }
    return Result;
}

mat2 operator*(mat2 &A, mat2 &B)

{
    mat2 Result = {};

    for(u8 RowIndex=0; RowIndex<2; RowIndex++)
    {
        for(u8 ColumnIndex=0; ColumnIndex<2; ColumnIndex++)
        {
            f32 Value =   GetMatrixElement(A, RowIndex, 0) * GetMatrixElement(B, 0, ColumnIndex) 
                        + GetMatrixElement(A, RowIndex, 1) * GetMatrixElement(B, 1, ColumnIndex);

            SetMatrixElement(&Result, RowIndex, ColumnIndex, Value);
        }
    }
    return Result;
}

mat3 operator*(mat3 &A, mat3 &B)

{
    mat3 Result = {};

    for(u8 RowIndex=0; RowIndex<3; RowIndex++)
    {
        for(u8 ColumnIndex=0; ColumnIndex<3; ColumnIndex++)
        {
            f32 Value =   GetMatrixElement(A, RowIndex, 0) * GetMatrixElement(B, 0, ColumnIndex) 
                        + GetMatrixElement(A, RowIndex, 1) * GetMatrixElement(B, 1, ColumnIndex) 
                        + GetMatrixElement(A, RowIndex, 2) * GetMatrixElement(B, 2, ColumnIndex);
            SetMatrixElement(&Result, RowIndex, ColumnIndex, Value);
        }
    }

    return Result;
}

mat4 operator*(mat4 &MatrixA, mat4 &MatrixB)

{
    mat4 Result = {};

    for(u8 RowIndex=0; RowIndex<4; RowIndex++)
    {
        for(u8 ColumnIndex=0; ColumnIndex<4; ColumnIndex++)
        {
            f32 A = GetMatrixElement(MatrixA, RowIndex, 0) * GetMatrixElement(MatrixB, 0, ColumnIndex);
            f32 B = GetMatrixElement(MatrixA, RowIndex, 1) * GetMatrixElement(MatrixB, 1, ColumnIndex);
            f32 C = GetMatrixElement(MatrixA, RowIndex, 2) * GetMatrixElement(MatrixB, 2, ColumnIndex);
            f32 D = GetMatrixElement(MatrixA, RowIndex, 3) * GetMatrixElement(MatrixB, 3, ColumnIndex);
            f32 Value = A + B + C + D;
            SetMatrixElement(&Result, RowIndex, ColumnIndex, Value);
        }
    }

    return Result;
}


v2 operator*(mat2 &Matrix, v2 &Vector)
{
    v2 Result = {};

    Result.x = GetMatrixElement(Matrix, 0, 0) * Vector.x + GetMatrixElement(Matrix, 0, 1) * Vector.y;
    Result.y = GetMatrixElement(Matrix, 1, 0) * Vector.x + GetMatrixElement(Matrix, 1, 1) * Vector.y;

    return Result;
}

v3 operator*(mat3 &Matrix, v3 &Vector)
{
    v3 Result = {};

    Result.x = GetMatrixElement(Matrix, 0, 0) * Vector.x + GetMatrixElement(Matrix, 0, 1) * Vector.y + GetMatrixElement(Matrix, 0, 2) * Vector.z;
    Result.y = GetMatrixElement(Matrix, 1, 0) * Vector.x + GetMatrixElement(Matrix, 1, 1) * Vector.y + GetMatrixElement(Matrix, 1, 2) * Vector.z;
    Result.z = GetMatrixElement(Matrix, 2, 0) * Vector.x + GetMatrixElement(Matrix, 2, 1) * Vector.y + GetMatrixElement(Matrix, 2, 2) * Vector.z;

    return Result;
}

v4 operator*(mat4 &Matrix, v4 &Vector)
{
    v4 Result = {};

    Result.x = GetMatrixElement(Matrix, 0, 0) * Vector.x + GetMatrixElement(Matrix, 0, 1) * Vector.y + GetMatrixElement(Matrix, 0, 2) * Vector.z + GetMatrixElement(Matrix, 0, 3) * Vector.w;
    Result.y = GetMatrixElement(Matrix, 1, 0) * Vector.x + GetMatrixElement(Matrix, 1, 1) * Vector.y + GetMatrixElement(Matrix, 1, 2) * Vector.z + GetMatrixElement(Matrix, 1, 3) * Vector.w;
    Result.z = GetMatrixElement(Matrix, 2, 0) * Vector.x + GetMatrixElement(Matrix, 2, 1) * Vector.y + GetMatrixElement(Matrix, 2, 2) * Vector.z + GetMatrixElement(Matrix, 2, 3) * Vector.w;
    Result.w = GetMatrixElement(Matrix, 3, 0) * Vector.x + GetMatrixElement(Matrix, 3, 1) * Vector.y + GetMatrixElement(Matrix, 3, 2) * Vector.z + GetMatrixElement(Matrix, 3, 3) * Vector.w;

    return Result;
}

f32 Determinant(mat2 Input)
{
    f32 Result = (0);

    f32 a = GetMatrixElement(Input, 0, 0);
    f32 b = GetMatrixElement(Input, 0, 1);
    f32 c = GetMatrixElement(Input, 1, 0);
    f32 d = GetMatrixElement(Input, 1, 1);

    Result = a * d - b * c;
    return Result;
}

f32 Determinant(mat3 Input)
{
   //Find cofactors
    f32 DetA = Determinant(SubMatrix(Input, 0, 0));
    f32 DetB = Determinant(SubMatrix(Input, 0, 1));
    f32 DetC = Determinant(SubMatrix(Input, 0, 2));

    //Find determinant
    f32 Result = GetMatrixElement(Input, 0, 0) * DetA - GetMatrixElement(Input, 0, 1) * DetB + GetMatrixElement(Input, 0, 2) * DetC;
    return Result;
}

f32 Determinant(mat4 Input)
{
   //Find cofactors
    f32 DetA = Determinant(SubMatrix(Input, 0, 0));
    f32 DetB = Determinant(SubMatrix(Input, 0, 1));
    f32 DetC = Determinant(SubMatrix(Input, 0, 2));
    f32 DetD = Determinant(SubMatrix(Input, 0, 3));

    //Find determinant
    f32 Result = GetMatrixElement(Input, 0, 0) * DetA - GetMatrixElement(Input, 0, 1) * DetB + GetMatrixElement(Input, 0, 2) * DetC - GetMatrixElement(Input, 0, 3) * DetD;
    return Result;
}

mat2 Inverse(mat2 Input)
{
    mat2 Result = {};

    f32 OneOverDet = 1.0f / Determinant(Input);
    
    
    mat2 SwappedMatrix = {};
    
    f32 a = GetMatrixElement(Input, 0, 0);
    f32 b = GetMatrixElement(Input, 0, 1);
    f32 c = GetMatrixElement(Input, 1, 0);
    f32 d = GetMatrixElement(Input, 1, 1);

    SetMatrixElement(&SwappedMatrix, 0, 0,  d);
    SetMatrixElement(&SwappedMatrix, 0, 1, -b);
    SetMatrixElement(&SwappedMatrix, 1, 0, -c);
    SetMatrixElement(&SwappedMatrix, 1, 1,  a);

    Result = OneOverDet * SwappedMatrix;
    return Result;
}

mat3 Inverse(mat3 Input)
{
    //Find cofactors
    f32 DetA = Determinant(SubMatrix(Input, 0, 0));
    f32 DetB = Determinant(SubMatrix(Input, 0, 1));
    f32 DetC = Determinant(SubMatrix(Input, 0, 2));

    //Find determinant
    f32 MatrixDeterminant = GetMatrixElement(Input, 0, 0) * DetA - GetMatrixElement(Input, 0, 1) * DetB + GetMatrixElement(Input, 0, 2) * DetC;
    f32 OneOverDeterminant = 1.0f / MatrixDeterminant;
    //Build cofactor matrix
    mat3 CofactorMatrix =  {};

	mat3 CofactorSigns = 
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
            f32 Cofactor = (0);
            if(RowIndex==0 && ColumnIndex==0) Cofactor = DetA;
            if(RowIndex==0 && ColumnIndex==1) Cofactor = DetB;
            if(RowIndex==0 && ColumnIndex==2) Cofactor = DetC;
            else Cofactor = Determinant(SubMatrix(Input, RowIndex, ColumnIndex));
            
            
            f32 Sign = GetMatrixElement(CofactorSigns, RowIndex, ColumnIndex);
            Cofactor = Cofactor * Sign;

            Cofactor = Cofactor * OneOverDeterminant;

            //Transposing in place!
            SetMatrixElement(&CofactorMatrix, ColumnIndex, RowIndex, Cofactor);
        }        
    }

    return CofactorMatrix;
}


mat4 Inverse(mat4 Input)
{
//First version is optimized code
//Second version is easy to read code
#if 1
    f32 a00 = GetMatrixElement(Input, 0, 0);
    f32 a01 = GetMatrixElement(Input, 0, 1);
    f32 a02 = GetMatrixElement(Input, 0, 2);
    f32 a03 = GetMatrixElement(Input, 0, 3);
    
    f32 a10 = GetMatrixElement(Input, 1, 0);
    f32 a11 = GetMatrixElement(Input, 1, 1);
    f32 a12 = GetMatrixElement(Input, 1, 2);
    f32 a13 = GetMatrixElement(Input, 1, 3);
    
    f32 a20 = GetMatrixElement(Input, 2, 0);
    f32 a21 = GetMatrixElement(Input, 2, 1);
    f32 a22 = GetMatrixElement(Input, 2, 2);
    f32 a23 = GetMatrixElement(Input, 2, 3);
    
    f32 a30 = GetMatrixElement(Input, 3, 0);
    f32 a31 = GetMatrixElement(Input, 3, 1);
    f32 a32 = GetMatrixElement(Input, 3, 2);
    f32 a33 = GetMatrixElement(Input, 3, 3);

    f32 A2323 = a22 * a33 - a23 * a32 ;
    f32 A1323 = a21 * a33 - a23 * a31 ;
    f32 A1223 = a21 * a32 - a22 * a31 ;
    f32 A0323 = a20 * a33 - a23 * a30 ;
    f32 A0223 = a20 * a32 - a22 * a30 ;
    f32 A0123 = a20 * a31 - a21 * a30 ;
    f32 A2313 = a12 * a33 - a13 * a32 ;
    f32 A1313 = a11 * a33 - a13 * a31 ;
    f32 A1213 = a11 * a32 - a12 * a31 ;
    f32 A2312 = a12 * a23 - a13 * a22 ;
    f32 A1312 = a11 * a23 - a13 * a21 ;
    f32 A1212 = a11 * a22 - a12 * a21 ;
    f32 A0313 = a10 * a33 - a13 * a30 ;
    f32 A0213 = a10 * a32 - a12 * a30 ;
    f32 A0312 = a10 * a23 - a13 * a20 ;
    f32 A0212 = a10 * a22 - a12 * a20 ;
    f32 A0113 = a10 * a31 - a11 * a30 ;
    f32 A0112 = a10 * a21 - a11 * a20 ;

    f32 det = a00 * ( a11 * A2323 - a12 * A1323 + a13 * A1223 ) 
        - a01 * ( a10 * A2323 - a12 * A0323 + a13 * A0223 ) 
        + a02 * ( a10 * A1323 - a11 * A0323 + a13 * A0123 ) 
        - a03 * ( a10 * A1223 - a11 * A0223 + a12 * A0123 ) ;
    det = 1 / det;

   f32 b00 = det *   ( a11 * A2323 - a12 * A1323 + a13 * A1223 );
   f32 b01 = det * - ( a01 * A2323 - a02 * A1323 + a03 * A1223 );
   f32 b02 = det *   ( a01 * A2313 - a02 * A1313 + a03 * A1213 );
   f32 b03 = det * - ( a01 * A2312 - a02 * A1312 + a03 * A1212 );
   f32 b10 = det * - ( a10 * A2323 - a12 * A0323 + a13 * A0223 );
   f32 b11 = det *   ( a00 * A2323 - a02 * A0323 + a03 * A0223 );
   f32 b12 = det * - ( a00 * A2313 - a02 * A0313 + a03 * A0213 );
   f32 b13 = det *   ( a00 * A2312 - a02 * A0312 + a03 * A0212 );
   f32 b20 = det *   ( a10 * A1323 - a11 * A0323 + a13 * A0123 );
   f32 b21 = det * - ( a00 * A1323 - a01 * A0323 + a03 * A0123 );
   f32 b22 = det *   ( a00 * A1313 - a01 * A0313 + a03 * A0113 );
   f32 b23 = det * - ( a00 * A1312 - a01 * A0312 + a03 * A0112 );
   f32 b30 = det * - ( a10 * A1223 - a11 * A0223 + a12 * A0123 );
   f32 b31 = det *   ( a00 * A1223 - a01 * A0223 + a02 * A0123 );
   f32 b32 = det * - ( a00 * A1213 - a01 * A0213 + a02 * A0113 );
   f32 b33 = det *   ( a00 * A1212 - a01 * A0212 + a02 * A0112 );

    mat4 Result = {};
    SetMatrixElement(&Result, 0, 0, b00);
    SetMatrixElement(&Result, 0, 1, b01);
    SetMatrixElement(&Result, 0, 2, b02);
    SetMatrixElement(&Result, 0, 3, b03);
    
    SetMatrixElement(&Result, 1, 0, b10);
    SetMatrixElement(&Result, 1, 1, b11);
    SetMatrixElement(&Result, 1, 2, b12);
    SetMatrixElement(&Result, 1, 3, b13);
    
    SetMatrixElement(&Result, 2, 0, b20);
    SetMatrixElement(&Result, 2, 1, b21);
    SetMatrixElement(&Result, 2, 2, b22);
    SetMatrixElement(&Result, 2, 3, b23);
    
    SetMatrixElement(&Result, 3, 0, b30);
    SetMatrixElement(&Result, 3, 1, b31);
    SetMatrixElement(&Result, 3, 2, b32);
    SetMatrixElement(&Result, 3, 3, b33);

    return Result;
    
#else
    //Find cofactors
    f32 DetA = Determinant(SubMatrix(Input, 0, 0));
    f32 DetB = Determinant(SubMatrix(Input, 0, 1));
    f32 DetC = Determinant(SubMatrix(Input, 0, 2));
    f32 DetD = Determinant(SubMatrix(Input, 0, 3));

    //Find determinant
    f32 MatrixDeterminant = GetMatrixElement(Input, 0, 0) * DetA - GetMatrixElement(Input, 0, 1) * DetB + GetMatrixElement(Input, 0, 2) * DetC - GetMatrixElement(Input, 0, 3) * DetD;
    f32 OneOverDeterminant = 1.0f / MatrixDeterminant;
    //Build cofactor matrix
    mat4 CofactorMatrix =  {};

	mat4 CofactorSigns = 
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
            f32 Cofactor = 0;
            if(RowIndex==0 && ColumnIndex==0) Cofactor = DetA;
            if(RowIndex==0 && ColumnIndex==1) Cofactor = DetB;
            if(RowIndex==0 && ColumnIndex==2) Cofactor = DetC;
            if(RowIndex==0 && ColumnIndex==3) Cofactor = DetD;
            else Cofactor = Determinant(SubMatrix(Input, RowIndex, ColumnIndex));
            
            
            f32 Sign = GetMatrixElement(CofactorSigns, RowIndex, ColumnIndex);
            Cofactor *= Sign;

            Cofactor *= OneOverDeterminant;

            //Transposing in place!
            SetMatrixElement(&CofactorMatrix, ColumnIndex, RowIndex, Cofactor);
        }        
    }
    
    return CofactorMatrix;
#endif
}


mat3 Translate(v2 Translation) 
{
    mat3 Result = Mat3F((1.0f));
    SetMatrixElement(&Result, 0, 2, Translation.x);
    SetMatrixElement(&Result, 1, 2, Translation.y);
    return Result;
}

mat3 Scale(v2 Size) 
{
    mat3 Result = Mat3F((1.0f));
    SetMatrixElement(&Result, 0, 0, Size.x);
    SetMatrixElement(&Result, 1, 1, Size.y);
    return Result;
}

mat3 Rotate(f32 Angle) 
{
    mat3 Result = Mat3F((1.0f));
    // SetMatrixElement(&Result, 0, 0, Cosine(Angle));
    // SetMatrixElement(&Result, 0, 1, Sine(Angle));
    // SetMatrixElement(&Result, 1, 0, -Sine(Angle));
    // SetMatrixElement(&Result, 1, 1, Cosine(Angle));
    return Result;
}


mat4 LookAt(v3 CameraPosition, v3 Center, v3 UpVector)
{
    v3 Z = NOZ(Center - CameraPosition);
    UpVector = LaneV3(0,1,0);
    
    u32 IsTooCloseMask = Abs(Inner(UpVector, Z)) > (0.99f); 
    ConditionalAssign(&UpVector, IsTooCloseMask, LaneV3(0,0,1));
    
    IsTooCloseMask = Abs(Inner(UpVector, Z)) > (0.99f); 
    ConditionalAssign(&UpVector, IsTooCloseMask, LaneV3(1,0,1));


    v3 X = NOZ(Cross(UpVector, Z));
    v3 Y = NOZ(Cross(Z, X));

    mat4 Result = {};

    SetMatrixElement(&Result, 0, 0, X.x);
    SetMatrixElement(&Result, 1, 0, X.y);
    SetMatrixElement(&Result, 2, 0, X.z);
    SetMatrixElement(&Result, 3, 0, (0));

    SetMatrixElement(&Result, 0, 1, Y.x);
    SetMatrixElement(&Result, 1, 1, Y.y);
    SetMatrixElement(&Result, 2, 1, Y.z);
    SetMatrixElement(&Result, 3, 1, (0));
    
    SetMatrixElement(&Result, 0, 2, Z.x);
    SetMatrixElement(&Result, 1, 2, Z.y);
    SetMatrixElement(&Result, 2, 2, Z.z);
    SetMatrixElement(&Result, 3, 2, (0));
    
    SetMatrixElement(&Result, 0, 3, CameraPosition.x);
    SetMatrixElement(&Result, 1, 3, CameraPosition.y);
    SetMatrixElement(&Result, 2, 3, CameraPosition.z);
    SetMatrixElement(&Result, 3, 3, (1));
    
    // mat4 Result = {
    //     X.x, Y.x, Z.x, CameraPosition.x,
    //     X.y, Y.y, Z.y, CameraPosition.y,
    //     X.z, Y.z, Z.z, CameraPosition.z,
    //     0  , 0  , 0  , 1
    // };

    return Result;
}

v3 TransformPosition(mat4 Matrix, v3 Vector)
{
    v3 Result = {};
    v4 ResultV4 = Matrix * V4(Vector.x, Vector.y, Vector.z, (1.0f));

    Result.x = ResultV4.x;
    Result.y = ResultV4.y;
    Result.z = ResultV4.z;

    return Result;
}

v3 TransformDirection(mat4 Matrix, v3 Vector)
{
    v3 Result = {};
    v4 ResultV4 = Matrix * V4(Vector.x, Vector.y, Vector.z, (0.0f));

    Result.x = ResultV4.x;
    Result.y = ResultV4.y;
    Result.z = ResultV4.z;

    return Result;
}

mat4 Identity()
{
    return Mat4F((1.0f));
}

mat4 OrthoBasisFromNormal(v3 Normal)
{
    v3 UpVector = LaneV3(0,1,0);
    
    // if(Abs(Inner(UpVector, Normal)) > 0.99)
    // {
    //     UpVector = LaneV3(0,0,1);
    // }
    
    // if(Abs(Inner(UpVector, Normal)) > 0.99)
    // {
    //     UpVector = LaneV3(1,0,0);
    // }

    mat4 Result = LookAt(LaneV3(0.0f, 0.0f, 0.0f), Normal, UpVector);

    return Result;
}

mat4 Translate(mat4 Matrix, v3 Position)
{
    mat4 Result = Matrix;
    SetMatrixElement(&Result, 0, 3, Position.x);
    SetMatrixElement(&Result, 1, 3, Position.y);
    SetMatrixElement(&Result, 2, 3, Position.z);
    return Result;
}
