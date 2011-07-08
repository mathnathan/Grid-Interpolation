%read a fortran unformatted binary 3-dimensional array of size n x m x l

n=540;
m=420;

fname1='/panfs/storage.local/coaps/home/ndc08/code/hycom2ncom/data/f90data/SMinterp_2d.dat';
fname2='/panfs/storage.local/coaps/home/ndc08/code/hycom2ncom/data/f90data/NCinterp_2d.dat';
fname3='/panfs/storage.local/coaps/home/ndc08/code/hycom2ncom/data/f90data/DIFinterp_2d.dat';
fname4='/panfs/storage.local/coaps/home/ndc08/code/hycom2ncom/data/f90data/longitude_2d.dat';
fname5='/panfs/storage.local/coaps/home/ndc08/code/hycom2ncom/data/f90data/latitude_2d.dat';

%open a file for reading...big endian (we try to always save data in big endian format)
file_smt=fopen(fname1,'r','b');
file_nct=fopen(fname2,'r','b');
file_dif=fopen(fname3,'r','b');
file_lon=fopen(fname4,'r','b');
file_lat=fopen(fname5,'r','b');

%fortran binary records start and end with a 32-bit record length in bytes..we just read it as a dummy variable since we know the array size
reclen=fread(file_smt,1,'int32'); 
reclen=fread(file_nct,1,'int32'); 
reclen=fread(file_dif,1,'int32'); 
reclen=fread(file_lon,1,'int32'); 
reclen=fread(file_lat,1,'int32'); 

%read a record of length n*m*l
smt=fread(file_smt,n*m,'float'); 
nct=fread(file_nct,n*m,'float'); 
dif=fread(file_dif,n*m,'float'); 
lon=fread(file_lon,n,'float'); 
lat=fread(file_lat,m,'float'); 

%read the record length again (we don't really need to do this unless we would be reading more records from the fortran binary file)
reclen=fread(file_smt,1,'int32'); 
reclen=fread(file_nct,1,'int32'); 
reclen=fread(file_dif,1,'int32'); 
reclen=fread(file_lon,1,'int32'); 
reclen=fread(file_lat,1,'int32'); 
fclose(file_smt);
fclose(file_nct);
fclose(file_dif);
fclose(file_lon);
fclose(file_lat);

%reshape the vector we just read into an nxmxl array
smt=reshape(smt,n,m); 
nct=reshape(nct,n,m); 
dif=reshape(dif,n,m); 

%matlab switches rows and columns, so we'll perumte the first and second dimensions
smt=permute(smt,[2 1 3]); 
nct=permute(nct,[2 1 3]); 
dif=permute(dif,[2 1 3]); 

%plot the surface temperature z=1
f = figure('visible', 'off');
clf;

subplot(2,1,1);
pcolor(smt(:,:)); shading flat; colorbar; 
title('Steve''s xy-slice');xlabel('Longitude');ylabel('Latitude');

subplot(2,1,2);
pcolor(nct(:,:)); shading flat; colorbar; 
title('Nathan''s 2d interpolation');xlabel('Longitude');ylabel('Latitude');

print( f, '-djpeg', '../data/graphs/xyslice_2d' )
clf;

pcolor(dif(:,:)); shading flat; colorbar; 
title('Steve''s xy-slice - Nathan''s 2d interpolation');xlabel('Longitude');ylabel('Latitude');

print( f, '-djpeg', '../data/graphs/xyzlice_2d_diff' )
clf;

disp('FINISHED')
