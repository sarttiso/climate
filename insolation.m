% Insolation Computer
% Generate insolation for any point on the globe. e, o, and w can be
% vectors but must all be of the same size
% from Huybers (PhD thesis)
%
% IN:
% lat: latitude for which to compute insolation, degrees. Can compute
%   insolation at a single latitude for multiple orbital parameters, or can 
%   compute insolation at multiple latitudes for fixed orbital parameters
% e: Eccentricity, vector or scalar
% o: Obliquity (degrees), vector or scalar
% w: Argument of the perigree/perihelion (degrees), vector or scalar
% 'I0': solar constant (default 1368) W/m2
% 'mode': 'instant' or 'daily', determines whether or not to average
%   daily insolation at a latitude or to calculate instantaneous insolation
%   at a latitutde and longitude. 
% 'slon': Solar longitude (degrees). For NH: 0-spring 90-summer 180-fall     
%   270-winter
% 'hour': Hour angle (degrees) (specify only for 'instant' mode), must be
%   scalar
% 'lon': longitude at which to compute insolation (specify only for 
%   'instant' mode), must be scalar
% 
% OUT:
% I: insolation
%
% Adrian Tasistro-Hart, 05.09.2018

function I = insolation(lat,e,o,w,varargin)
    
%% parse and validate

parser = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(parser,'lat',@isnumeric)
addRequired(parser,'e',@isnumeric)
addRequired(parser,'o',@isnumeric);
addRequired(parser,'w',@isnumeric);
addParameter(parser,'I0',1368,validScalarPosNum)
addParameter(parser,'mode','daily',@ischar);
addParameter(parser,'slon',0,@isnumeric);
addParameter(parser,'hour',0,@isnumeric);
addParameter(parser,'lon',0,@isnumeric);

parse(parser,lat,e,o,w,varargin{:})

lat = parser.Results.lat;
e = parser.Results.e;
o = parser.Results.o;
w = parser.Results.w;
I0 = parser.Results.I0;
mode = parser.Results.mode;
slon = parser.Results.slon;
hour = parser.Results.hour;
lon = parser.Results.lon;

% validate number of latitudes vs number of orbital parameters
nlat = length(lat);
ne = length(e);
no = length(o);
nw = length(w);
norbit = max([ne; no; nw]);
% varying orbital parameters 
if nlat == 1
    % validate eccentricity
    if ne == 1
        e = e*ones(norbit,1);
    else
        assert(ne == norbit,'inappropriate number of eccentricity values')
    end
    % validate obliquity
    if no == 1
        o = o*ones(norbit,1);
    else
        assert(no == norbit,'inappropriate number of obliquity values')
    end
    % validate precession
    if nw == 1
        w = w*ones(norbit,1);
    else
        assert(nw == norbit,'inappropriate number of precession values')
    end
% varying latitudes
elseif nlat > 1
    assert(ne == 1 && no == 1 && nw == 1, ...
        ['for multiple latitudes, specify only single values for each'...
        ' orbital parameter'])
    e = e*ones(nlat,1);
    o = o*ones(nlat,1);
    w = w*ones(nlat,1);
else
    error('must specify at least one latitude')
end

% validate mode
mode = validatestring(mode,{'daily','instant'});
% if we need hour, lon parameters, make sure they're the correct size
if strcmp(mode,'instant')
    hour = hour*ones(max([norbit,nlat]),1);
    lon = lon*ones(max([norbit,nlat]),1);
end

%% convert values and compute insolation

% convert to radians
lat = deg2rad(lat);
lon = deg2rad(lon);
o = deg2rad(o);
w = deg2rad(w);
slon = deg2rad(slon);
hour = deg2rad(hour);

% eccentricity component
R = ((1+e.*cos(slon-w))./(1-e.^2)).^2;
% solar declination
delta = asin(sin(o)*sin(slon));

switch mode
    case 'instant'
    I = I0 * R .* ...
        (sin(lat).*sin(delta) + cos(lat).*cos(delta).*cos(lon-hour));
    case 'daily'
    % hour angle at sunset
    H = real(acos(-tan(lat)*tan(delta)));
    I = I0/pi * R .* ...
        (H.*sin(lat).*sin(delta) + cos(lat).*cos(delta).*sin(H));
end

end