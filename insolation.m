% Insolation Computer
% Generate insolation for any point on the globe. e, o, and w can be
% vectors but must all be of the same size
% I0: Solar constant
% e: Eccentricity
% o: Obliquity (degrees)
% w: Argument of the perigree/perihelion (degrees)
% mode: 'instant' or 'daily', determines whether or not to average
%   daily insolation at a latitude or to calculate instantaneous insolation
%   at a latitutde and longitude. Instant requires specifing lon and n
% slon: Solar longitude (degrees). For NH: 0-spring 90-summer 180-fall     
%   270-winter
% n: Hour angle (degrees)
% from Huybers

function I = insolation(lat,e,o,w,mode,slon,n,lon,I0)
    
    % default values
    defval('I0',1368)       % Solar flux, 1368 W/m2
    defval('e',0.01672)     % Modern eccentricity
    defval('o',23.5)        % Modern obliquity º
    defval('w',90)          % Southern hemisphere summer solstice at peri.
    defval('slon',0)        % Earth at vernal equinox 
    defval('mode','daily')  % Daily average insolation
    defval('n',0)           % Prime meridian
    defval('lon',0);        % Prime meridian
    
    % convert to radians
    lat = deg2rad(lat);
    lon = deg2rad(lon);
    o = deg2rad(o);
    w = deg2rad(w);
    slon = deg2rad(slon);
    n = deg2rad(n);
    
    % eccentricity component
    R = ((1+e.*cos(slon-w))./(1-e.^2)).^2;
    % solar declination
    delta = asin(sin(o)*sin(slon));
    
    if(strcmp(mode,'instant'))
        I = I0 * R .* ...
            (sin(lat).*sin(delta) + cos(lat).*cos(delta).*cos(lon-n));
    elseif(strcmp(mode,'daily'))
        % hour angle at sunset
        H = real(acos(-tan(lat)*tan(delta)));
        I = I0/pi * R .* ...
            (H.*sin(lat).*sin(delta) + cos(lat).*cos(delta).*sin(H));
    else
        error('Mode was specified incorrectly.')
    end

end