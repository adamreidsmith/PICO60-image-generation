
function [ret_array, ray_startingpoints, pixels, generator_params] = GetRaysAndPixels(allparams, resolutions)
% Adapted from EricsQuickLookupTableMakerEtcWithTorus.m

%% check inputs

if iscell(allparams)
    allparams = cell2mat(allparams);
end

if iscell(resolutions)
    resolutions = cell2mat(resolutions);
end

max_scatters = 4;

if nargin<2
    disp('You need 2 arguments when calling GetPixelCoordsFrom3DPoint');
    return
end

if size(resolutions,2)~=2
    disp('Resolutions input is weird in GetPixelCoordsFrom3DPoint');
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

resolutions = repmat(resolutions, n_cam, 1);

if n_cam ~= round(n_cam)
    disp('Improper allparams input -- length gives non-integer number of cameras');
    return
end

if n_cam ~= size(resolutions,1)
    disp('Improper resolutions input -- length disagrees with allparams!');
    return
end

cam(n_cam).x = [];


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

n_air = allparams(param_offset + 11);
n_window = allparams(param_offset + 12);
n_hydraulic = allparams(param_offset + 13);
n_jar = allparams(param_offset + 14);
n_target = allparams(param_offset + 15);

z_offset = allparams(param_offset + 16);

%% make pixel list
pixels_camends = cumsum(prod(resolutions,2));
pixels_camstarts = [1;(pixels_camends(1:end-1)+1)];

pixels = zeros(pixels_camends(end),3);
for i_cam = 1:n_cam
    pixels(pixels_camstarts(i_cam):pixels_camends(i_cam),1) = i_cam;
    pixels(pixels_camstarts(i_cam):pixels_camends(i_cam),2) = reshape(repmat((1:resolutions(i_cam,1))',1,resolutions(i_cam,2)),[],1);
    pixels(pixels_camstarts(i_cam):pixels_camends(i_cam),3) = reshape(repmat(1:resolutions(i_cam,2),resolutions(i_cam,1),1),[],1);
end

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
        fprintf('Watch out, geometry is out of bounds for air-side of window for cam %d',i_cam);
    end
    raylist(this_ix,:) = RefractionReflectionAtInterface(raylist(this_ix,:), normals, n_air, n_window);
    
    [ray_startingpoints(this_ix,:), normals, lt, or] = RayToPlane(ray_startingpoints(this_ix,:), raylist(this_ix,1:3), window_hydraulicsidepoint, window_normal);
    if any(isinf(lt) | isnan(lt) | (lt <= 0) | (or ~= -1))
        fprintf('Watch out, geometry is out of bounds for hydraulic-side of window for cam %d',i_cam);
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

if jar_xpitch~=0 || jar_ypitch~=0
    jar_pitch = acos(cyl_axis(3));
    jar_yaw = atan2(jar_ypitch,jar_xpitch);

    hemi_rotmat = [cos(jar_pitch), 0, -sin(jar_pitch) ; 0, 1, 0 ; sin(jar_pitch), 0, cos(jar_pitch)] * ...
        [cos(jar_yaw), sin(jar_yaw), 0 ; -sin(jar_yaw), cos(jar_yaw), 0 ; 0, 0, 1];
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

%% Propagate rays through vessel

raytracer_output = RayTracer(ray_startingpoints, ...
    raylist, surface_list, max_scatters, 1e-5, [0 100]);

%%

total_num_scat_pts = 0;
for i=1:length(raytracer_output)
    total_num_scat_pts = total_num_scat_pts + size(raytracer_output(i).ray_index,1);
end

ret_array = zeros(total_num_scat_pts,4);

n = 1;
for i=1:length(raytracer_output)
    % Scattering points in this iteration
    scatterpoints = raytracer_output(i).intersection_point;
    % Indices of the scattering points
    indices = raytracer_output(i).ray_index;
    
    ret_array(n:n+size(scatterpoints)-1, 1:3) = scatterpoints;
    ret_array(n:n+size(scatterpoints)-1, 4) = indices;
    
    n = n + size(scatterpoints,1);
end

generator_params = [r1 r2 r3 liquidlevel cyl_axis];

