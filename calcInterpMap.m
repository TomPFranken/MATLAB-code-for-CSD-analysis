function map = calcInterpMap(M,N_INTERP_X,N_INTERP_Y,method)

Msiz = size(M);

if ~exist('N_INTERP_X','var')
    N_INTERP_X = 50;
end

if ~exist('N_INTERP_Y','var')
    N_INTERP_Y = 50;
end

if ~exist('method','var')
    method = 'cubic';
end

map.Xlim = linspace(1,Msiz(2),N_INTERP_X);
map.Ylim = linspace(1,Msiz(1),N_INTERP_Y);

[XI,YI] = meshgrid(map.Xlim, map.Ylim);

map.ZI = interp2(M,XI,YI,method);

% map.ZI = imgaussfilt(map.ZI,4); %SD=4 seems good for 192 channels from Neuropixels

