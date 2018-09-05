% Currently uses Laskar 2010a for eccentricity and Laskar 2004 for
% obliquity and precession.
% lat: latitude (degrees)
% t: time vector in kyr, equally spaced. Laskar's solutions will be
%   interpolated to match the spacing of the time vector. Must be negative,
%   0 is present day.
% slon: Solar longitude, see insolation
% I0: Solar constant
% e2010: Which Laskar simulation to use for eccentricity (a,b,c,d)


function I = insolation_orbital(lat,t,slon,I0,e2010)
    defval('e2010','a')
    la04 = load('/Users/adrianraph/Documents/research/_published_data/orbital_solutions/INSOLN.LA2004.BTL.250.txt');
    la10 = load(sprintf('/Users/adrianraph/Documents/research/_published_data/orbital_solutions/La2010%s_ecc3L.dat',...
        e2010));
    
    defval('I0',1368);
    
    idx = la04(:,1) >= min(t) & la04(:,1) <= max(t);
    e = la10(idx,2);
    o = rad2deg(la04(idx,3));
    w = rad2deg(la04(idx,4));
    
    % interpolate
    wrap = @(V) mod(V, 2*pi);
    e = interp1(la10(idx,1),e,t,'spline');
    o = interp1(la04(idx,1),o,t,'spline');
    w = interp1(la04(idx,1),unwrap(deg2rad(w)),t,'spline');
    w = rad2deg(wrap(w));
    
    % note that Berger's equations utilizes the longitude of 
    I = insolation(lat,e,o,w+180,'daily',slon,[],[],I0);
end