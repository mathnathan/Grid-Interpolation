%read a fortran unformatted binary 3-dimensional array of size n x m x l

n=540;
m=420;
l=80;

fname1='/panfs/storage.local/coaps/home/ndc08/code/hycom2ncom/data/f90data/SMinterp.dat';
fname2='/panfs/storage.local/coaps/home/ndc08/code/hycom2ncom/data/f90data/NCinterp.dat';
fname3='/panfs/storage.local/coaps/home/ndc08/code/hycom2ncom/data/f90data/DIFinterp.dat';
fname4='/panfs/storage.local/coaps/home/ndc08/code/hycom2ncom/data/f90data/longitude.dat';
fname5='/panfs/storage.local/coaps/home/ndc08/code/hycom2ncom/data/f90data/latitude.dat';
fname6='/panfs/storage.local/coaps/home/ndc08/code/hycom2ncom/data/f90data/depth.dat';

%open a file for reading...big endian (we try to always save data in big endian format)
file_smt=fopen(fname1,'r','b');
file_nct=fopen(fname2,'r','b');
file_dif=fopen(fname3,'r','b');
file_lon=fopen(fname4,'r','b');
file_lat=fopen(fname5,'r','b');
file_z=fopen(fname6,'r','b');

%fortran binary records start and end with a 32-bit record length in bytes..we just read it as a dummy variable since we know the array size
reclen=fread(file_smt,1,'int32'); 
reclen=fread(file_nct,1,'int32'); 
reclen=fread(file_dif,1,'int32'); 
reclen=fread(file_lon,1,'int32'); 
reclen=fread(file_lat,1,'int32'); 
reclen=fread(file_z,1,'int32'); 

%read a record of length n*m*l
smt=fread(file_smt,n*m*l,'float'); 
nct=fread(file_nct,n*m*l,'float'); 
dif=fread(file_dif,n*m*l,'float'); 
lon=fread(file_lon,n,'float'); 
lat=fread(file_lat,m,'float'); 
z=fread(file_z,n*m*l,'float'); 

%read the record length again (we don't really need to do this unless we would be reading more records from the fortran binary file)
reclen=fread(file_smt,1,'int32'); 
reclen=fread(file_nct,1,'int32'); 
reclen=fread(file_dif,1,'int32'); 
reclen=fread(file_lon,1,'int32'); 
reclen=fread(file_lat,1,'int32'); 
reclen=fread(file_z,1,'int32'); 
fclose(file_smt);
fclose(file_nct);
fclose(file_dif);
fclose(file_lon);
fclose(file_lat);
fclose(file_z);

%reshape the vector we just read into an nxmxl array
smt=reshape(smt,n,m,l); 
nct=reshape(nct,n,m,l); 
dif=reshape(dif,n,m,l); 
z=reshape(z,n,m,l);

%matlab switches rows and columns, so we'll perumte the first and second dimensions
smt=permute(smt,[2 1 3]); 
nct=permute(nct,[2 1 3]); 
dif=permute(dif,[2 1 3]); 

%plot the surface temperature z=1
f = figure('visible', 'off');
clf;

subplot(2,1,1);
pcolor(smt(:,:,1)); shading flat; colorbar; 
title('Steve''s xy-slice z=1');xlabel('Longitude');ylabel('Latitude');

subplot(2,1,2);
pcolor(nct(:,:,1)); shading flat; colorbar; 
title('Nathan''s xy-slice z=1');xlabel('Longitude');ylabel('Latitude');

print( f, '-djpeg', '../data/graphs/xyslice_z=1' )
clf;

pcolor(dif(:,:,1)); shading flat; colorbar; 
title('Steve''s - Nathan''s xy-slice z=1');xlabel('Longitude');ylabel('Latitude');

print( f, '-djpeg', '../data/graphs/xyslice_z=1_diff' )
clf;

%plot the surface temperature z=40
subplot(2,1,1);
pcolor(smt(:,:,40)); shading flat; colorbar; 
title('Steve''s xy-slice z=40');xlabel('Longitude');ylabel('Latitude');

subplot(2,1,2);
pcolor(nct(:,:,40)); shading flat; colorbar; 
title('Nathan''s xy-slice z=40');xlabel('Longitude');ylabel('Latitude');

print( f, '-djpeg', '../data/graphs/xyslice_z=40' )
clf;

pcolor(dif(:,:,40)); shading flat; colorbar; 
title('Steve''s - Nathan''s xy-slice z=40');xlabel('Longitude');ylabel('Latitude');

print( f, '-djpeg', '../data/graphs/xyslice_z=40_diff' )
clf;

%plot the surface temperature z=80
subplot(2,1,1);
pcolor(smt(:,:,80)); shading flat; colorbar; 
title('Steve''s xy-slice z=80');xlabel('Longitude');ylabel('Latitude');

subplot(2,1,2);
pcolor(nct(:,:,80)); shading flat; colorbar; 
title('Nathan''s xy-slice z=80');xlabel('Longitude');ylabel('Latitude');

print( f, '-djpeg', '../data/graphs/xyslice_z=80' )
clf;

pcolor(dif(:,:,80)); shading flat; colorbar; 
title('Steve''s - Nathan''s xy-slice z=80');xlabel('Longitude');ylabel('Latitude');

print( f, '-djpeg', '../data/graphs/xyslice_z=80_diff' )
clf;

%plot an x-z slice along row 50
%keyboard
[LN,dmm]=meshgrid(lon,dimension(z,3));

subplot(2,1,1);
A=squeeze(smt(50,:,:))';
pcolor(lon,squeeze(z(1,50,:)),A); shading interp; colorbar; 
title('Steve''s xz-slice y=50');xlabel('Longitude');ylabel('Depth');

subplot(2,1,2);
pcolor(rot90(squeeze(nct(:,50,:)))); shading interp; colorbar; 
title('Nathan''s xz-slice y=50');xlabel('Longitude');ylabel('Depth');

print( f, '-djpeg', '../data/graphs/xzslice_y=50' )
clf;

pcolor(rot90(squeeze(dif(:,50,:)))); shading interp; colorbar; 
title('Steve''s - Nathan''s xz-slice y=50');xlabel('Longitude');ylabel('Depth');

print( f, '-djpeg', '../data/graphs/xzslice_y=50_diff' )
clf;

%plot an x-y slice along column 50
subplot(2,1,1);
pcolor(rot90(squeeze(smt(50,:,:)))); shading interp; colorbar; 
title('Steve''s yz-slice x=50');xlabel('Latitude');ylabel('Depth');

subplot(2,1,2);
pcolor(rot90(squeeze(nct(50,:,:)))); shading interp; colorbar; 
title('Nathan''s yz-slice x=50');xlabel('Latitude');ylabel('Depth');

print( f, '-djpeg', '../data/graphs/yzslice_x=50' )
clf;

pcolor(rot90(squeeze(dif(50,:,:)))); shading interp; colorbar; 
title('Steve''s - Nathan''s yz-slice x=50');xlabel('Latitude');ylabel('Depth');

print( f, '-djpeg', '../data/graphs/yzslice_x=50_diff' )
clf;

disp('FINISHED');
