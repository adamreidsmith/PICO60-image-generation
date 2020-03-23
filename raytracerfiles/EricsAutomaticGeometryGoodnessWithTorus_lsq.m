% function residuals = EricsAutomaticGeometryGoodnessWithTorus_lsq(fitparams, allparams, whichparamstofit, ...
%     pixels, izphi, pixel_err, zcoordtype)
%
% This function is a stab at a chisq calculator for the COUPP60 (or other
% planar-viewport geometry) bubble chamber.  It is structured so that it
% can be used by an anonymous function with any subset of the (many)
% geometry parameters fixed vs floating.  It takes some shortcuts (compared
% to the typical RayTracer) that the user should keep in mind.  These
% include:
%
% - The rays followed from each camera are propagated by hand through the
% viewport before ever being sent to RayTracer -- if the fit routine moves
% a camera to the wrong side of its window, there may be problems.
%
% - After entering the hydraulic fluid, all rays are propagated through a
% 6-surface geometry by RayTracer.  The four surfaces are the inside and
% outside of the jar hemisphere/ellipsoid, jar torus knuckle, and jar
% cylinder.  There is no buffer fluid in this geometry, so you cannot look
% at fiducial marks above the liquid level
%
% The chisq calculation is done by following pixels to a given surface and
% comparing their intersection point at that surface to the expected
% spatial position.  Examples:
%
% A fiducial mark at (phi,z) on the front side of the jar should get i=1 in
% the (i,z,phi) triplet.  This tells the code to look at the first
% intersection point for the pixel that hits that fiducial mark (outside of
% the jar) and compare it to the true position
%
% If z or phi is NaN, then the coordinate is unrestrained in that
% dimension.  I.e., z is NaN if the aziumuth of a point is known but the
% height is not.
%
% A point on meniscus should get z=Inf, phi=NaN and i=2 or 3.  z=Inf
% is a code indicating that z is at the liquidlevel.
%
% A point at the bottom of the jar should get phi=NaN, z=-Inf and i=3 or 4.
%
% A fiducial mark on the back side of the jar should get i=4
%
% Rays that have a special relation to a surface, that is they are tangent
% to it or intersect it at a critical angle, are indicated by z=phi=NaN and
% i<0.  The value of i indicates the specific relation: i=-1 for tangent to
% the outer surface of the jar, i=-3 for tangent to the inner surface of
% the jar, i=-4 for meeting the inner surface of the jar at the critical
% angle for total internal reflection.
%
% One can also simply require that a ray intersect n surfaces by setting
% both phi=NaN and z=NaN.  e.g., a ray that should undergo total-internal
% reflection off the inner-wall of the jar should have
% [i,z,phi]=[3,NaN,NaN]
%
% pixel_err is in length-units, same units as used in allparams (typically
% cm).  For points where z=phi=NaN and i>=0, pixel_err has units of
% #-of-surfaces intersected.  For points where z=phi=NaN and i=-1 or i=-3,
% pixel_err has units of length^2 (interpreted as the distance between
% intersection points for the ray and the jar surface).  For z=phi=NaN and
% i=-4. pixel_err is unitless (residual is delta(sin(theta)))
%
% Finally, you can pick out pairs of (pixel,interface) coordinates and
% require that they be spatially close to one another.   To do this, list
% the two pixels consecutively in the 'pixels' list.  For the first pixel,
% make i<0, z=+1, phi=Inf -- for the second, i<0, z=-1, phi=Inf.  The
% absolute value of i should still be the index of the scattering surface.
%
%
% 2013-08-26, CED

function residuals = EricsAutomaticGeometryGoodnessWithTorus_lsq(fitparams, allparams, whichparamstofit, ...
    pixels, izphi, pixel_err, zcoordtype)

%% check inputs

residuals = [];

max_scatters = 4;

if nargin<5
    disp('You need 5 arguments when calling EricsAutomaticGeometryGoodness');
    return
end

if nargin<6 || isempty(pixel_err)
    pixel_err = ones(size(pixels,1),1);
end

if nargin<7 || isempty(zcoordtype)
    zcoordtype = 'arc';
end

pixel_err = pixel_err(:);
fitparams = fitparams(:);
allparams = allparams(:);
whichparamstofit = whichparamstofit(:);

if size(izphi,2)~=3 || size(pixels,2)~=3 || size(izphi,1)~=size(pixels,1) || size(izphi,1)~=length(pixel_err)
    disp('Improper pixel or izphi or pixel_err input to EricsAutomaticGeometryGoodness');
    return
end

if ~islogical(whichparamstofit) || length(whichparamstofit)~=length(allparams) || length(fitparams)~=sum(whichparamstofit)
    disp('Improper fitparams or allparams or whichparamstofit input to EricsAutomaticGeometryGoodness');
    return
end

if any(abs(izphi(:,1))>max_scatters)
    disp('bad i''s in izphi, in input to EricsAutomaticGeometryGoodness');
    return
end

%% write down params
cam = struct('x',[],'y',[],'z',[],'pitch',[],'yaw',[],'roll',[], ...
    'f',[],'f1',[],'f2',[],'i0',[],'j0',[],'i_pitch',[],'j_pitch',[],'theta',[],'phi',[],'bf',[],'lens_type',[], ...
    'win_d',[],'win_t',[],'win_phi',[],'win_pitch',[],'win_yaw',[]);

fn = fieldnames(cam);
n_camparams = length(fn);

n_otherparams = 16;

n_cam = (length(allparams)-n_otherparams)/n_camparams;

if n_cam ~= round(n_cam)
    disp('Improper allparams input -- length gives non-integer number of cameras');
    return
end

cam(n_cam).x = [];



allparams(whichparamstofit) = fitparams;


for i_cam = 1:n_cam
    for i_f = 1:n_camparams
        cam(i_cam).(fn{i_f}) = allparams((i_cam-1)*n_camparams+i_f);
    end
end

param_offset = n_cam*n_camparams;

jar_phi = allparams(param_offset + 1);
liquidlevel = allparams(param_offset + 2);

jar_wall = allparams(param_offset + 3);
jar_OD = allparams(param_offset + 4);
jar_ID = jar_OD - 2*jar_wall;

jar_axwall = allparams(param_offset + 5);
jar_Oaxrad = allparams(param_offset + 6);
jar_Iaxrad = jar_Oaxrad - jar_axwall;

jar_knucklewall = allparams(param_offset + 7);
jar_knuckleOrad = allparams(param_offset + 8);
jar_knuckleIrad = jar_knuckleOrad - jar_knucklewall;

jar_xpitch = allparams(param_offset + 9);
jar_ypitch = allparams(param_offset + 10);
%%% add torus geospecs here

n_air = allparams(param_offset + 11);
n_window = allparams(param_offset + 12);
n_hydraulic = allparams(param_offset + 13);
n_jar = allparams(param_offset + 14);
n_target = allparams(param_offset + 15);

z_offset = allparams(param_offset + 16);

%% create rays from cameras and propogate through window

raylist = repmat([0 1 0 0 0 1 1 0 0 0],size(pixels,1),1);
ray_startingpoints = repmat([0 0 0],size(pixels,1),1);

for i_cam = 1:n_cam
    this_ix = find(pixels(:,1)==i_cam);
    ccd_ijk = (pixels(this_ix,[2,3,1]) - repmat([cam(i_cam).i0, cam(i_cam).j0, 0],length(this_ix),1) ) .* ...
        repmat([-cam(i_cam).i_pitch, cam(i_cam).j_pitch, 0],length(this_ix),1);
    ccd_tiltmat_A = [ ...
        cos(cam(i_cam).phi), sin(cam(i_cam).phi), 0 ; ...
        -sin(cam(i_cam).phi), cos(cam(i_cam).phi), 0 ; ...
        0, 0, 1];
    ccd_tiltmat_B = [ ...
        cos(cam(i_cam).theta), 0, sin(cam(i_cam).theta) ; ...
        0, 1, 0 ; ...
        -sin(cam(i_cam).theta), 0, cos(cam(i_cam).theta) ] ;
    ccd_ijk = ccd_ijk * ccd_tiltmat_A' * ccd_tiltmat_B' * ccd_tiltmat_A;
    
    ccd_ijk(:,1:2) = ccd_ijk(:,1:2) .* repmat(cam(i_cam).bf ./ (cam(i_cam).bf - ccd_ijk(:,3)),1,2);
        
    ccd_d2 = sum(ccd_ijk(:,1:2).^2,2);
    barrel_d = [cam(i_cam).f1, cam(i_cam).f2];
    
    effective_f = cam(i_cam).f .* (1 + ...
        sum( repmat(barrel_d,length(ccd_d2),1) .* ...
        (repmat(cam(i_cam).f.^-2 .* ccd_d2,1,length(barrel_d)).^repmat(1:length(barrel_d),length(ccd_d2),1)), 2) );

    switch cam(i_cam).lens_type
        case 1 %'theta'
            theta = sqrt(ccd_d2)./effective_f;
        case 2 %'sin'
            theta = asin(sqrt(ccd_d2)./effective_f);
        case 3 %'tan'
            theta = atan(sqrt(ccd_d2)./effective_f);
        otherwise
            theta = atan(sqrt(ccd_d2)./effective_f);
    end
    
    raylist(this_ix,2) = cos(theta);
    raylist(this_ix,[1 3]) = -(ccd_ijk(:,1:2)./sqrt(repmat(ccd_d2,1,2))).*repmat(sin(theta),1,2);
    
    M1 = [cos(cam(i_cam).yaw) -sin(cam(i_cam).yaw) 0 ; sin(cam(i_cam).yaw) cos(cam(i_cam).yaw) 0 ; 0 0 1];
    M2 = [1 0 0 ; 0 cos(cam(i_cam).pitch) -sin(cam(i_cam).pitch) ; 0 sin(cam(i_cam).pitch) cos(cam(i_cam).pitch)];
    M3 = [cos(cam(i_cam).roll) 0 sin(cam(i_cam).roll) ; 0 1 0 ; -sin(cam(i_cam).roll) 0 cos(cam(i_cam).roll)];

    M = M1*M2*M3;

    raylist(this_ix,1:3) = (M * ((raylist(this_ix,1:3))'))';
    ray_startingpoints(this_ix,1) = cam(i_cam).x;
    ray_startingpoints(this_ix,2) = cam(i_cam).y;
    ray_startingpoints(this_ix,3) = cam(i_cam).z;
    
    ray_startingpoints(this_ix, 1:2) = [ray_startingpoints(this_ix, 1), ray_startingpoints(this_ix, 2)-cam(i_cam).win_d] * ...
        [cos(cam(i_cam).win_yaw), sin(cam(i_cam).win_yaw) ; -sin(cam(i_cam).win_yaw), cos(cam(i_cam).win_yaw)];
    
    window_airsidepoint = [cam(i_cam).win_d*sin(cam(i_cam).win_phi), ...
        -cam(i_cam).win_d*cos(cam(i_cam).win_phi), ...
        cam(i_cam).z];
        
    window_normal = [-cos(cam(i_cam).win_pitch)*sin(cam(i_cam).win_yaw) , ...
        cos(cam(i_cam).win_pitch)*cos(cam(i_cam).win_yaw) , ...
        sin(cam(i_cam).win_pitch)];
    
    window_hydraulicsidepoint = window_airsidepoint + cam(i_cam).win_t * window_normal;
    
    [ray_startingpoints(this_ix,:), normals, lt, or] = RayToPlane(ray_startingpoints(this_ix,:), raylist(this_ix,1:3), window_airsidepoint, window_normal);
    if any(isinf(lt) | isnan(lt) | (lt <= 0) | (or ~= -1))
        fprintf(2,'Watch out, geometry is out of bounds for air-side of window for cam %d\n',i_cam);
    end
    raylist(this_ix,:) = RefractionReflectionAtInterface(raylist(this_ix,:), normals, n_air, n_window);
    
    [ray_startingpoints(this_ix,:), normals, lt, or] = RayToPlane(ray_startingpoints(this_ix,:), raylist(this_ix,1:3), window_hydraulicsidepoint, window_normal);
    if any(isinf(lt) | isnan(lt) | (lt <= 0) | (or ~= -1))
        fprintf(2,'Watch out, geometry is out of bounds for hydraulic-side of window for cam %d\n',i_cam);
    end
    raylist(this_ix,:) = RefractionReflectionAtInterface(raylist(this_ix,:), normals, n_window, n_hydraulic);
    
end

%% Build surface list (ignore buffer fluid in this version)

cyl_axis = [sin(jar_xpitch), sin(jar_ypitch), sqrt(1-sin(jar_xpitch)^2-sin(jar_ypitch)^2)];

r1 = [jar_OD, jar_ID]*.5;
r2 = [jar_knuckleOrad, jar_knuckleIrad];
r3 = [jar_Oaxrad, jar_Iaxrad];
s = r3.*(r1-r2)./(r3-r2);
z = -r2 .* sqrt(1 - (s./r3).^2);
d = r3 .* z .* ((1./r3)-(1./r2));

% quartz_hemi_inside_Q = [(.5*jar_ID)^-2 0 0 ; 0 (.5*jar_ID)^-2 0 ; 0 0 (jar_Iaxrad)^-2 ];
% quartz_hemi_outside_Q = [(.5*jar_OD)^-2 0 0 ; 0 (.5*jar_OD)^-2 0 ; 0 0 (jar_Oaxrad)^-2 ];

if jar_xpitch~=0 || jar_ypitch~=0
    jar_pitch = acos(cyl_axis(3));
    jar_yaw = atan2(jar_ypitch,jar_xpitch);

    hemi_rotmat = [cos(jar_pitch), 0, -sin(jar_pitch) ; 0, 1, 0 ; sin(jar_pitch), 0, cos(jar_pitch)] * ...
        [cos(jar_yaw), sin(jar_yaw), 0 ; -sin(jar_yaw), cos(jar_yaw), 0 ; 0, 0, 1];

%     quartz_hemi_inside_Q = (hemi_rotmat')*quartz_hemi_inside_Q*hemi_rotmat;
%     quartz_hemi_outside_Q = (hemi_rotmat')*quartz_hemi_outside_Q*hemi_rotmat;
else
    hemi_rotmat = [1,0,0;0,1,0;0,0,1];
end

surface_list = struct( ...
    'description', {}, ...
    'intersect_function', {}, ...
    'inbounds_function', {}, ...
    'n_outside', {}, ...
    'n_inside', {}, ...
    'surface_type', {}, ...
    'absorption', {});

surface_list(end+1).description = 'inside surface of quartz cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], cyl_axis, r1(2));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:)*cyl_axis(1) + p(:,2,:)*cyl_axis(2) + p(:,3,:)*cyl_axis(3)) >= 0, size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz cylinder';
surface_list(end).intersect_function = @(sp,indir)RayToCylinder(sp,indir, ...
    [0 0 0], cyl_axis, r1(1));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:)*cyl_axis(1) + p(:,2,:)*cyl_axis(2) + p(:,3,:)*cyl_axis(3)) >= 0, size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of quartz dome';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    cyl_axis*d(2), r3(2));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:)*cyl_axis(1) + p(:,2,:)*cyl_axis(2) + p(:,3,:)*cyl_axis(3)) < z(2), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz dome';
surface_list(end).intersect_function = @(sp,indir)RayToSphere(sp,indir, ...
    cyl_axis*d(1), r3(1));
surface_list(end).inbounds_function = @(p)(reshape( (p(:,1,:)*cyl_axis(1) + p(:,2,:)*cyl_axis(2) + p(:,3,:)*cyl_axis(3)) < z(1), size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'inside surface of quartz knuckle';
surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
    [0 0 0], cyl_axis, r1(2)-r2(2), r2(2));
surface_list(end).inbounds_function = @(p)(reshape( ...
    ((p(:,1,:)*cyl_axis(1) + p(:,2,:)*cyl_axis(2) + p(:,3,:)*cyl_axis(3)) < 0) & ...
    ((p(:,1,:)*cyl_axis(1) + p(:,2,:)*cyl_axis(2) + p(:,3,:)*cyl_axis(3)) >= z(2)) & ...
    ((p(:,1,:).^2+p(:,2,:).^2+p(:,3,:).^2 - (p(:,1,:)*cyl_axis(1) + p(:,2,:)*cyl_axis(2) + p(:,3,:)*cyl_axis(3)).^2)>((r1(2)-r2(2))^2)), size(p,1), [] ));
surface_list(end).n_outside = n_jar;
surface_list(end).n_inside = n_target;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;

surface_list(end+1).description = 'outside surface of quartz knuckle';
surface_list(end).intersect_function = @(sp,indir)RayToTorus(sp,indir, ...
    [0 0 0], cyl_axis, r1(1)-r2(1), r2(1));
surface_list(end).inbounds_function = @(p)(reshape( ...
    ((p(:,1,:)*cyl_axis(1) + p(:,2,:)*cyl_axis(2) + p(:,3,:)*cyl_axis(3)) < 0) & ...
    ((p(:,1,:)*cyl_axis(1) + p(:,2,:)*cyl_axis(2) + p(:,3,:)*cyl_axis(3)) >= z(1)) & ...
    ((p(:,1,:).^2+p(:,2,:).^2+p(:,3,:).^2 - (p(:,1,:)*cyl_axis(1) + p(:,2,:)*cyl_axis(2) + p(:,3,:)*cyl_axis(3)).^2)>((r1(1)-r2(1))^2)), size(p,1), [] ));
surface_list(end).n_outside = n_hydraulic;
surface_list(end).n_inside = n_jar;
surface_list(end).surface_type = 'normal';
surface_list(end).absorption = 0;


%% propagate rays through vessel

raytracer_output = RayTracer(ray_startingpoints, ...
    raylist, surface_list, max_scatters, 1e-5, [0 100]);

%% handle rays that disappear
intersection_matrix = zeros(size(raylist,1),3,max_scatters);
intersection_matrix(:) = NaN;

for n = 1:length(raytracer_output)
    intersection_matrix(raytracer_output(n).ray_index, :, n) = raytracer_output(n).intersection_point;
end

num_scatters = sum(squeeze(all(~isnan(intersection_matrix),2)),2);

surface_targets = izphi(:,1)>0 & isinf(izphi(:,2)) & isnan(izphi(:,3)) & izphi(:,2)>0;
axis_targets = izphi(:,1)>0 & isinf(izphi(:,2)) & isnan(izphi(:,3)) & izphi(:,2)<0;
edge_targets = izphi(:,1)>=0 & isnan(izphi(:,2)) & isnan(izphi(:,3));
edge_targets_smooth = izphi(:,1)<0 & isnan(izphi(:,2)) & isnan(izphi(:,3));
twopixel_targets_A = izphi(:,1)<0 & izphi(:,2)>0 & isinf(izphi(:,3));
twopixel_targets_B = izphi(:,1)<0 & izphi(:,2)<0 & isinf(izphi(:,3));
wall_targets = ~(surface_targets | axis_targets | edge_targets | edge_targets_smooth | twopixel_targets_A | twopixel_targets_B);

problematic_rays = ~(edge_targets | edge_targets_smooth) & num_scatters<abs(izphi(:,1));

if any(problematic_rays)
    if 1
        fprintf(2,'Uh oh, we found %d rays that didn''t make it far enough in the geometry\n',sum(problematic_rays));
    end
    intersection_matrix(problematic_rays,:,:) = 1e3;
end

%% calculate chisq
target_positions = zeros(size(raylist,1), 3);
sim_positions = zeros(size(target_positions));

% z_offset = 0;
switch zcoordtype
    case 'arc'
        if (.5*jar_OD)>jar_Oaxrad
            hemi_e = sqrt(1 - (jar_Oaxrad*2/jar_OD)^2);
            hemi_a = .5*jar_OD;
            [~, hemi_quartercirc] = ellipke(hemi_e);
            z_offset = hemi_quartercirc*hemi_a;
        else
            hemi_e = sqrt(1-(jar_OD/(2*jar_Oaxrad))^2);
            hemi_a = jar_Oaxrad;
            [~, hemi_quartercirc] = ellipke(hemi_e);
            z_offset = hemi_quartercirc*hemi_a;
        end
    case 'elev'
        z_offset = jar_Oaxrad - d(1);
    case 'float'
%     otherwise
%         z_offset = 0;
end
        

izphi(:,3) = izphi(:,3)/(.5*jar_OD);
izphi(:,2) = izphi(:,2) - z_offset;

for n=1:max_scatters
    sim_positions(abs(izphi(:,1))==n,:) = intersection_matrix(abs(izphi(:,1))==n,:,n);
end
sim_positions(edge_targets,3) = num_scatters(edge_targets);


sim_positions(edge_targets_smooth,3) = 0;

tan_outerwall = edge_targets_smooth & izphi(:,1)==-1;
% crit_outerwall = edge_targets_smooth & ixphi(:,2)==-2;
tan_innerwall = edge_targets_smooth & izphi(:,1)==-3;
crit_innerwall = edge_targets_smooth & izphi(:,1)==-4;

if any(tan_outerwall)
    [ip_cyl, ~, dt_cyl] = surface_list(2).intersect_function(ray_startingpoints(tan_outerwall,:), ...
        raylist(tan_outerwall,1:3));
    ib_cyl = surface_list(2).inbounds_function(real(ip_cyl));
    [ip_hemi, ~, dt_hemi] = surface_list(4).intersect_function(ray_startingpoints(tan_outerwall,:), ...
        raylist(tan_outerwall,1:3));
    ib_hemi = surface_list(4).inbounds_function(real(ip_hemi));
    [~, ~, dt_tor] = surface_list(6).intersect_function(ray_startingpoints(tan_outerwall,:), ...
        raylist(tan_outerwall,1:3));
    
    res_cyl = real(diff(dt_cyl,1,2).^2);
    res_hemi = real(diff(dt_hemi,1,2).^2);
    res_tor_complex = [ ...
        dt_tor(:,2)-dt_tor(:,1), ...
        dt_tor(:,3)-dt_tor(:,1), ...
        dt_tor(:,4)-dt_tor(:,1), ...
        dt_tor(:,3)-dt_tor(:,2), ...
        dt_tor(:,4)-dt_tor(:,2), ...
        dt_tor(:,4)-dt_tor(:,3), ...
        ].^2;
    res_tor_real = real(res_tor_complex);
    res_tor_good = abs(imag(res_tor_complex)) < 1e-10;
    res_tor_real(~res_tor_good) = inf;
    [~, which_res_tor] = min(abs(res_tor_real), [], 2);
    which_res_tor_cut = repmat(1:6, length(which_res_tor), 1) == repmat(which_res_tor, 1, 6);
    res_tor_real = res_tor_real';
    which_res_tor_cut = which_res_tor_cut';
    res_tor = res_tor_real(which_res_tor_cut);
    
    
    cyl_cut = any(ib_cyl, 2);
    hemi_cut = ~cyl_cut & any(ib_hemi, 2);
    tor_cut = ~(cyl_cut | hemi_cut);

    res_all = res_cyl;    
    res_all(hemi_cut) = res_hemi(hemi_cut);
    res_all(tor_cut) = res_tor(tor_cut);
    
    sim_positions(tan_outerwall,3) = res_all;
    
end

if any(tan_innerwall)
    tan_innerwall_good = tan_innerwall & num_scatters>=2;
    tan_innerwall_bad = tan_innerwall & ~tan_innerwall_good;
    if any(tan_innerwall_good)
        tan_innerwall_ix = find(tan_innerwall_good);
        [sub_tf, sub_ix] = ismember(tan_innerwall_ix,raytracer_output(1).ray_index);
        if ~all(sub_tf)
            disp('We have a problem in tan_innerwall');
        end
        [ip_cyl, ~, dt_cyl] = surface_list(1).intersect_function(raytracer_output(1).intersection_point(sub_ix,:), ...
            raytracer_output(1).refracted_ray(sub_ix,1:3));
        ib_cyl = surface_list(1).inbounds_function(real(ip_cyl));
        [ip_hemi, ~, dt_hemi] = surface_list(3).intersect_function(raytracer_output(1).intersection_point(sub_ix,:), ...
            raytracer_output(1).refracted_ray(sub_ix,1:3));
        ib_hemi = surface_list(3).inbounds_function(real(ip_hemi));
        [~, ~, dt_tor] = surface_list(5).intersect_function(raytracer_output(1).intersection_point(sub_ix,:), ...
            raytracer_output(1).refracted_ray(sub_ix,1:3));
        res_cyl = real(diff(dt_cyl,1,2).^2);
        res_hemi = real(diff(dt_hemi,1,2).^2);
        res_tor_complex = [ ...
            dt_tor(:,2)-dt_tor(:,1), ...
            dt_tor(:,3)-dt_tor(:,1), ...
            dt_tor(:,4)-dt_tor(:,1), ...
            dt_tor(:,3)-dt_tor(:,2), ...
            dt_tor(:,4)-dt_tor(:,2), ...
            dt_tor(:,4)-dt_tor(:,3), ...
            ].^2;
        res_tor_real = real(res_tor_complex);
        res_tor_good = abs(imag(res_tor_complex)) < 1e-10;
        res_tor_real(~res_tor_good) = inf;
        [~, which_res_tor] = min(abs(res_tor_real), [], 2);
        which_res_tor_cut = repmat(1:6, length(which_res_tor), 1) == repmat(which_res_tor, 1, 6);
        res_tor_real = res_tor_real';
        which_res_tor_cut = which_res_tor_cut';
        res_tor = res_tor_real(which_res_tor_cut);

        cyl_cut = any(ib_cyl, 2);
        hemi_cut = ~cyl_cut & any(ib_hemi, 2);
        tor_cut = ~(cyl_cut | hemi_cut);

        res_all = res_cyl;    
        res_all(hemi_cut) = res_hemi(hemi_cut);
        res_all(tor_cut) = res_tor(tor_cut);

        sim_positions(tan_innerwall_good,3) = res_all;
    end
    if any(tan_innerwall_bad)
        sim_positions(tan_innerwall_bad,3) = -100;
        fprintf(2,'Whoops, %d tan_innerwall rays didn''t get far enough in the geometry\n',sum(tan_innerwall_bad));
    end
end

if any(crit_innerwall)
    crit_innerwall_good = crit_innerwall & num_scatters>=3;
    crit_innerwall_bad = crit_innerwall & ~crit_innerwall_good;
    crit_sin_incident_angle = n_target / n_jar;
    if any(crit_innerwall_good)
        crit_innerwall_ix = find(crit_innerwall_good);
        [sub_tf, sub_ix] = ismember(crit_innerwall_ix,raytracer_output(2).ray_index);
        if ~all(sub_tf)
            fprintf(2,'We have a problem in crit_innerwall\n');
        end
        
        sin_incident_angle = sqrt(1-sum(raytracer_output(2).incoming_ray(sub_ix,1:3).*raytracer_output(2).surface_normal(sub_ix,:),2).^2);
        sim_positions(crit_innerwall_good,3) = crit_sin_incident_angle - sin_incident_angle;
        
    end
    if any(crit_innerwall_bad)
        sim_positions(crit_innerwall_bad,3) = crit_sin_incident_angle - 1;
        fprintf(2,'Whoops, %d crit_innerwall rays didn''t get far enough in the geometry\n',sum(crit_innerwall_bad));
    end
end


% remember, izphi is in jar axis, but target_positions and sim_positions
% are in physical axes (relevant if jar is rotated)

% first liquid surface targets
target_positions(surface_targets,1:2) = NaN;
target_positions(surface_targets,3) = liquidlevel;

% then axial targets
target_positions(axis_targets,:) = ...
    repmat([0 0 d(1)-jar_Oaxrad]*hemi_rotmat,sum(axis_targets),1);

% now discrete edge targets
target_positions(edge_targets,1:2) = NaN;
target_positions(edge_targets,3) = izphi(edge_targets,1);

% now smooth edge targets
target_positions(edge_targets_smooth,1:2) = NaN;

% now two-pixel targets
target_positions(twopixel_targets_A,:) = sim_positions(twopixel_targets_B,:);
target_positions(twopixel_targets_B,:) = NaN;

% now the 'normal' targets, which incidentally we assume are on the outside
% of the jar...
cyl_targets = wall_targets & (izphi(:,2)>=0 | isnan(izphi(:,2)));
sph_targets = wall_targets & izphi(:,2) < z(1);
tor_targets = wall_targets & ~(cyl_targets | sph_targets);

s_targets = .5*jar_OD + zeros(size(wall_targets));
s_targets(sph_targets) = sqrt(jar_Oaxrad.^2 - (d(1)-izphi(sph_targets,2)).^2);
s_targets(tor_targets) = r1(1) - r2(1) + sqrt(r2(1).^2-izphi(tor_targets,2).^2);

target_positions(wall_targets,:) = [s_targets(wall_targets).*cos(izphi(wall_targets,3)+jar_phi-(pi/2)), ...
    s_targets(wall_targets).*sin(izphi(wall_targets,3)+jar_phi-(pi/2)), ...
    izphi(wall_targets,2)];

sim_positions(wall_targets, :) = sim_positions(wall_targets, :) * (hemi_rotmat');

residuals = sim_positions - target_positions;

residuals = residuals ./ repmat(pixel_err,1,3);

residuals = residuals';

residuals = residuals(~isnan(residuals));

% residuals(isnan(residuals)) = 0;

% chisq = sum(sum(residuals.^2,2)./pixel_err);
 