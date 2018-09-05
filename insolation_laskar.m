% Currently uses Laskar 2010(a-d) for eccentricity and Laskar 2004 for
% obliquity and precession.
% 
% IN:
% lat: latitude (degrees). if multiple latitudes and times are requested, 
%   the output insolation will be an array with columns corresponding to
%   latitudes and rows corresponding to times
% t: time vector in kyr. Laskar's solutions will be interpolated to match 
%   the spacing of the time vector. Must be negative, 0 is present day.
% 'e2010': (default 'a') Which Laskar simulation to use for eccentricity, 
%   can be 'a','b','c','d'
% 
% any name-value pair parameters that can be given to insolation can be
% passed to this function as well, except those corresponding to
% instantaneous insolation, i.e. 'mode' is always 'daily'. 
%
% OUT: 
% I: insolation
%
% Adrian Tasistro-Hart 05.09.2018

function I = insolation_laskar(lat,t,varargin)

%% parse and validate
parser = inputParser;
parser.KeepUnmatched = true;
parser.PartialMatching = false;

addRequired(parser,'lat',@isnumeric)
addRequired(parser,'t',@isnumeric)
addParameter(parser,'e2010','a',@ischar)

parse(parser,lat,t,varargin{:})

lat = parser.Results.lat;
t = parser.Results.t;
e2010 = parser.Results.e2010;

% make columns
lat = lat(:);
t = t(:);

nlat = length(lat);
nt = length(t);

% validate e2010
e2010 = validatestring(e2010,{'a','b','c','d'});

%% load data and prepare for insolation.m
la04 = load('INSOLN.LA2004.BTL.250.txt');
la10 = load(sprintf('La2010%s_ecc3L.dat',e2010));

idx = la04(:,1) >= min(t) & la04(:,1) <= max(t);
e = la10(idx,2);
o = rad2deg(la04(idx,3));
w = rad2deg(la04(idx,4));

% interpolate
wrap = @(V) mod(V, 2*pi);
e = interp1(la10(idx,1),e,t,'linear');
o = interp1(la04(idx,1),o,t,'linear');
w = interp1(la04(idx,1),unwrap(deg2rad(w)),t,'linear');
w = rad2deg(wrap(w));

I = zeros(nt,nlat);
for ii = 1:nlat
    % note that Berger's equations utilizes the longitude of 
    I(:,ii) = insolation(lat(ii),e,o,w+180,'mode','daily',...
        parser.Unmatched);
end

end