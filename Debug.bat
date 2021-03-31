@echo off

REM SRC FILES
set  srcFiles=
REM --------------------

@REM set compilerFlags= /D UNICODE /D _UNICODE /DEBUG -D_CRT_SECURE_NO_WARNINGS -MP -MT -nologo -EHa- -Od -Oi -W4 -Z7 -EHsc -wd4201 -wd4310 -wd4100  /I ../src  /I ..\include\
@REM set linkerFlags=  -opt:ref Gdi32.lib Shell32.lib  %d3d12% %dxgi%  %dxguid% %assimpLib%


REM ASSIMP
set assimpLibDir= ..\ext\assimp\lib
set assimpLib=..\ext\assimp\lib\assimp-vc140-mt.lib
set assimpIncludes=..\ext\assimp\include
set assimpBin= ..\ext\assimp\bin\assimp-vc140-mt.dll
REM --------------------


set compilerFlags=-arch:AVX2 -Od -MTd -nologo -fp:fast -fp:except- -Gm- -GR- -EHa- -Zo -Oi -WX -W0 -wd4201 -wd4100 -wd4189 -wd4505 -wd4127 -FC -Z7
set compilerFlags=-DCOMPILER_MSVC -D_CRT_SECURE_NO_WARNINGS -DDEBUG=1 -DRAY_WIN32=1 %compilerFlags% /I %assimpIncludes%
set linkerFlags= -incremental:no -opt:ref user32.lib gdi32.lib winmm.lib opengl32.lib %assimpLib%

IF NOT EXIST .\build mkdir .\build
pushd .\build
"C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29910\bin\Hostx64\x64\cl.exe" %compilerFlags% ..\src\Main.cpp /link %linkerFlags%

copy %assimpBin% .\ >nul

popd
