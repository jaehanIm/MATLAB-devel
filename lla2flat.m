%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lla2flat
%
% Im, Jaehan, NearthLab, 200301
% Last modified date: 200305
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% work reference : https://kr.mathworks.com/help/aerotbx/ug/lla2flat.html

function FLAT = lla2flat(lla, llo, psio, href, ~)

d2r = pi/180;
psio = psio * d2r;
len = size(lla,1);
FLAT = zeros(len,3);
o_lat = llo(1); o_lon = llo(2); % reference location of ned [2]
o_lat = o_lat*d2r; o_lon = o_lon*d2r;

for i = 1:len
    % Prepocessing
    lat = lla(i,1); lon = lla(i,2); alt = lla(i,3);
    lat = lat*d2r; lon = lon*d2r;
    
    d_lat = lat - o_lat;
    d_lon = lon - o_lon;
    
    % WGS84 parameters
    a = 6378137;
    f = 0.003352810664747;
    
    % Geodetic lat lon to North East
    RN = a/sqrt(1-(2*f-f^2)*sin(o_lat)^2);
    RM = RN * (1-(2*f-f^2))/(1-(2*f-f^2)*sin(o_lat)^2);
    NED = zeros(3,1);
    NED(1) = d_lat/atan(1/RM);
    NED(2) = d_lon/atan(1/RN/cos(o_lat));
    NED(3) = -alt -href;

    % NED to Flat
    T = [ cos(psio)  , sin(psio) , 0  ;
          -sin(psio) , cos(psio) , 0  ;
          0          , 0         , 1 ];

    FLAT(i,:) = T * NED;
end

end