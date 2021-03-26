
///////////////////
//IMAGE
///////////////////
#pragma pack(push, 1)
struct bitmap_header
{
    u16 FileType;
    u32 FileSize;
    u16 Reserved1;
    u16 Reserved2;
    u32 BitmapOffset;
    u32 Size;
    s32 Width;
    s32 Height;
    u16 Planes;
    u16 BitsPerPixel;
    u32 Compression;
    u32 SizeOfBitmap;
    s32 HorizontalResolution;
    s32 VerticalResolution;
    u32 ColorsUsed;
    u32 ColorsImportant;
};
#pragma pack(pop)

struct image_32 
{
    u32 Width;
    u32 Height;
    u32 *Pixels;
};


internal u32 GetTotalPixelSize(image_32 Image) 
{
    return sizeof(u32) * Image.Width * Image.Height;
}

internal void WriteImage(image_32 Image, char *OutputFileName) 
{
    u32 OutputPixelSize = GetTotalPixelSize(Image);

    bitmap_header Header = {};
    Header.FileType = 0x4D42;
    Header.FileSize = sizeof(Header) + OutputPixelSize;
    Header.BitmapOffset = sizeof(Header);
    Header.Size = sizeof(Header) - 14;
    Header.Width = Image.Width;
    Header.Height = Image.Height;
    Header.Planes = 1;
    Header.BitsPerPixel = 32;
    Header.Compression = 0;
    Header.SizeOfBitmap = OutputPixelSize;
    Header.HorizontalResolution = 0;
    Header.VerticalResolution = 0;
    Header.ColorsUsed = 0;
    Header.ColorsImportant = 0;

    FILE *OutFile = fopen(OutputFileName, "wb+");
    if(OutFile) 
    {
        fwrite(&Header, sizeof(Header), 1, OutFile);
        fwrite(Image.Pixels, OutputPixelSize, 1, OutFile);
        fclose(OutFile);
    } else {
        fprintf(stderr, "ERROR : UNable to write output file %s. \n", OutputFileName);
    }
}

internal image_32 AllocateImage(u32 Width, u32 Height) {
    image_32 Image = {};
    Image.Width = Width;
    Image.Height = Height;
    u32 OutputPixelSize = GetTotalPixelSize(Image);
    Image.Pixels = (u32*)malloc(OutputPixelSize);
    return Image;
}

internal f32 LinearToSRGB(f32 L) {
    f32 S;
    if(L < 0.0f) L=0.0f;
    if(L > 1.0f) L=1.0f;

    S = L * 12.92f;
    if(L > 0.0031308) {
        S = 1.055f * Pow(L, 1.0f / 2.4f) - 0.055f;
    }

    return S;
}

internal u32 *GetPixelPointer(image_32 Image, u32 X, u32 Y) {
    u32 *Result = Image.Pixels + X + Y * Image.Width;
    return Result;
}
