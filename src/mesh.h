#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags


shape* MeshFromFile(char *FileName, u32 *TriangleNumber, u32 MaterialIndex, mat4 Transform)
{
    Assimp::Importer import;
    const aiScene *scene = import.ReadFile(FileName, aiProcess_Triangulate | aiProcess_MakeLeftHanded  | aiProcess_ValidateDataStructure| aiProcess_FlipUVs); 

    v3 *Vertices = nullptr;
    v3 *Normals = nullptr;
    u32 *Indices = nullptr;

    u32 CurrentVerticesCount=0;
    u32 CurrentIndicesCount=0;

    if(!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) 
    {
        printf("ModelLoader:LoaderModel: ERROR::\n"); 
        return nullptr;
    }


    for(unsigned int j = 0; j < scene->mNumMeshes; j++) {
        aiMesh *aiMesh = scene->mMeshes[j]; 
        printf("ModelLoader:LoaderModel: Vertices number : %d", aiMesh->mNumVertices);
        
        Vertices = (v3*)realloc(Vertices, (CurrentVerticesCount + aiMesh->mNumVertices) * sizeof(v3));
        Normals = (v3*)realloc(Normals, (CurrentVerticesCount + aiMesh->mNumVertices) * sizeof(v3));
        for(unsigned int i = 0; i < aiMesh->mNumVertices; i++)
        {
            Vertices[i] = {
                aiMesh->mVertices[i].x,
                aiMesh->mVertices[i].y,
                aiMesh->mVertices[i].z
            };

            Normals[i] = {
                aiMesh->mNormals[i].x,
                aiMesh->mNormals[i].y,
                aiMesh->mNormals[i].z
            };
            
            CurrentVerticesCount ++;
        }
        
        
        Indices = (u32*)realloc(Indices, (CurrentIndicesCount + (3 * aiMesh->mNumFaces)) * sizeof(u32));
        for(unsigned int i = 0; i < aiMesh->mNumFaces; i++)
        {
            aiFace face = aiMesh->mFaces[i];
            Indices[CurrentIndicesCount++] = face.mIndices[0];
            Indices[CurrentIndicesCount++] = face.mIndices[1];
            Indices[CurrentIndicesCount++] = face.mIndices[2];
        }
    }

    u32 TotalNumber = CurrentIndicesCount / 3;

    shape *OutputTriangles = (shape*) malloc(TotalNumber * sizeof(shape));
    u32 TrianglesAdded=0;
    for(u32 i=0; i<CurrentIndicesCount; i+=3)
    {
        u32 Index0 = Indices[i];
        u32 Index1 = Indices[i + 1];
        u32 Index2 = Indices[i + 2];

        v3 Vertex0 = Vertices[Index0];
        v3 Vertex1 = Vertices[Index1];
        v3 Vertex2 = Vertices[Index2];

        v3 Normal0 = Normals[Index0];
        v3 Normal1 = Normals[Index1];
        v3 Normal2 = Normals[Index2];

        *(OutputTriangles + TrianglesAdded) = Triangle(Vertex0, Vertex1, Vertex2, Normal0, Normal1, Normal2, Transform, MaterialIndex);
        TrianglesAdded++;
    }

    *TriangleNumber = TrianglesAdded;


    free(Indices);
    free(Vertices);
    free(Normals);


    return OutputTriangles;
}