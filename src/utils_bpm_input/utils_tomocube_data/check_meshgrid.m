x = linspace(0,1,15);
y = linspace(0,2,20);
z = linspace(-1,1,10);

[Xm,Ym,Zm] = meshgrid(x,y,z);
[Xn,Yn,Zn] = ndgrid(x,y,z);

isequal(Xm,pagetranspose(Xn))
isequal(Xn,pagetranspose(Xm))
isequal(Ym,pagetranspose(Yn))
isequal(Yn,pagetranspose(Ym))
isequal(Zm,pagetranspose(Zn))
isequal(Zn,pagetranspose(Zm))