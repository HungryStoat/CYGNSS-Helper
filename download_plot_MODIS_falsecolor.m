function [] = download_plot_MODIS_falsecolor(latitude_limits,longitude_limits,date_i_want,username,password,out_path)
%DOWNLOAD_PLOT_MODIS_FALSECOLOR This function will download MODIS/Terra
%Surface Reflectance data and plot it as a false color image. The surface
%reflectance files will then be deleted because oh boy they are kind of
%large files! This is because this function downloads 500 m resolution data, so if
%you want to plot the whole world...maybe don't. This function downloads
%version 6.0 of the MODIS data, which are imminently going to be replaced
%with version 6.1. If the code breaks, this is likely why. You'll have to
%update the paths to the files once this happens yourself. You will need at
%least the mapping toolbox for this to work; I think that's all though.

%INPUTS
% latitude_limits: For the region you want to download and plot, what are
% the latitude limits of the bounding box, in degrees? For example, if you want to plot
% the majority of the United States, you might input [25 50]. Do not include
% North/South designations, just use negative numbers for southern
% latitudes.
% longitude_limits: Same as latitude limits, but for longitude. Use the -180
% to 180 degrees convention, not 0 - 360.
% date_i_want: What date do you want to plot in 'yyyymmdd' format. For
% example, if you want to plot data from Aug 27, 2022, you would input
% '20220827'
% 'username': You need to create a free Earthdata login here:
% https://urs.earthdata.nasa.gov/, in single quotes
% 'password': Password for your Earthdata account, in single quotes
% out_path: Where do you want everything saved on your computer? For
% example, if I wanted things saved straight to my desktop (which I do not
% recommend) I would type: '/Users/username/Desktop/'
% Include the single quotes and all appropriate forward slashes, please.

%OUTPUT
% There is no output other than a saved figure in the out_path directory.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Citing the MODIS/Terra data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vermote, E., Wolfe, R. (2015). MOD09GA MODIS/Terra Surface Reflectance Daily L2G Global 1kmand 500m SIN Grid V006 [Data set].
% NASA EOSDIS Land Processes DAAC. Accessed 2022-09-22 from https://doi.org/10.5067/MODIS/MOD09GA.006

% Determine which MODIS tiles you'll want to download. Do I
% remember how I made this file many moons ago? Not really. I think it was
% derived from here: https://modis-land.gsfc.nasa.gov/MODLAND_grid.html
load('modis_tiles.mat');
h_and_v = [];
for ll = 1:length(tilehelp)
    [xg,yg] = meshgrid(tilehelp(ll,3):1:tilehelp(ll,4),tilehelp(ll,5):1:tilehelp(ll,6));
    inqd = ingeoquad(yg(:),xg(:),latitude_limits,longitude_limits);
    if(~isempty(find(inqd==1)))
        h_and_v = vertcat(h_and_v,[tilehelp(ll,2) tilehelp(ll,1)]);
    end
end

[numtiles,~] = size(h_and_v);

% Path name for where MODIS data live, etc.
modpath = 'https://e4ftl01.cr.usgs.gov/MOLT/MOD09GA.006/';
options = weboptions('Username',username,'Password',password,'Timeout',40);

% Separate the date
myyr = date_i_want(1:4);
mymo = date_i_want(5:6);
mydy = date_i_want(7:8);

% Make a directory to save all files to, most will get deleted later
mkdir([out_path date_i_want '/']);

%% Loop through the tiles you need
data = webread([modpath myyr '.' mymo '.' mydy '/'],options);
for tt = 1:numtiles
    % Name of the tile
    tname = ['h' sprintf('%02d',h_and_v(tt,1)) 'v' sprintf('%02d',h_and_v(tt,2))];
    kx = strfind(data,tname);
    for uu = 1:length(kx)
        myfilename = data(kx(uu)-17:kx(uu)+30);
        rf = strfind(myfilename,'</a');
        if(~isempty(rf))
            myfilename = myfilename(1:end-3);
            break;
        end
    end
    % This is the path for where the hdf file lives online
    url = [modpath myyr '.' mymo '.' mydy '/' myfilename];

    % Save .hdf reflectance file
    fname = websave([out_path date_i_want '/' myfilename],url,options);

    % Assign a grid of x and y values, to be converted into lat/lon later
    info = hdfinfo(fname,'eos');
    xcorr = info.Grid(2).UpperLeft(1);
    ycorr = info.Grid(2).UpperLeft(2);
    xcorr2 = info.Grid(2).LowerRight(1);
    ycorr2 = info.Grid(2).LowerRight(2);

    xvec = linspace(xcorr,xcorr2,2400);
    yvec = linspace(ycorr,ycorr2,2400);
    [xmesh,ymesh] = meshgrid(xvec,yvec);

    % Got the projection parameters below from here: https://modis-land.gsfc.nasa.gov/GCTP.html
    % The code below will convert x,y into lat/lon. Another way to do this would
    % be to cannibalize code found here: http://hdfeos.org/zoo/NSIDC/MYD10A1F_A2020131_h18v03_061_2020335131245_hdf.m
    radius = 6371007.181;

    mstruct = defaultm('sinusoid');
    mstruct.geoid = [1 0];
    mstruct.origin = [0 0];
    mstruct.falsenorthing = 0; mstruct.falseeasting = 0;
    mstruct = defaultm(mstruct);
    [lats,lons]=projinv(mstruct,xmesh./radius,ymesh./radius);

    % Determine if there are reflectance points in your region of interest
    inqd = ingeoquad(lats(1:5:end),lons(1:5:end),latitude_limits,longitude_limits);
    hj = find(inqd==1);
    % If so, save the bands for the false color image
    if(~isempty(hj))
        % Read bands 7, 2, and 1, for false color image
        data7 = hdfread(fname,'sur_refl_b07_1');
        data2 = hdfread(fname,'sur_refl_b02_1');
        data1 = hdfread(fname,'sur_refl_b01_1');

        % Include 0.0001 scaling factor in reflectance data
        data7 = double(data7).*0.0001;
        data2 = double(data2).*0.0001;
        data1 = double(data1).*0.0001;

        % False color image data
        [rs,cs] = size(data7);
        image = nan(rs,cs,3);
        image(:,:,1) = data7;
        image(:,:,2) = data2;
        image(:,:,3) = data1;
        
        % Make single to make files smaller, doesn't really matter
        lats = single(lats); lons = single(lons);

        % Save false color information for an individual tile, these will
        % be deleted later
        save([out_path date_i_want '/' myfilename '_forplot.mat'],'image','lats','lons');
    end
    % Remove large reflectance .hdf file
    delete([out_path date_i_want '/' myfilename]);
end

%% Now that you've saved false color information for individual tiles, plot all the tiles in one figure
% You will also be deleting the false color information for that tile
files = dir([[out_path date_i_want '/'],'*_forplot.mat']);
figure; set(gcf,'visible','off');
ax = worldmap(latitude_limits,longitude_limits);
for ff = 1:length(files)
    fname = [out_path date_i_want '/' files(ff).name];
    load(fname);
    geoshow(lats,lons,image);
    delete(fname); % Delete false color information for the tile. Comment out if you want to keep it.
end
exportgraphics(gcf,[out_path date_i_want '/FalseColor_' date_i_want '.png'],'Resolution',1000); % Print the figure at 1000 dpi.
close; % Close the figure;
end