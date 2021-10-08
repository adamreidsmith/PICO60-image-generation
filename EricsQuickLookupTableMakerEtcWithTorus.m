% function chisq = EricsAutomaticGeometryGoodness(fitparams, allparams, whichparamstofit, pixels, izphi, pixel_err)
%
%
% 2013-08-26, CED

function EricsQuickLookupTableMakerEtcWithTorus(allparams, resolutions, output_filename, imagefiles)

%% check inputs

max_scatters = 4;

if nargin<4 || ~iscell(imagefiles)
    imagefiles = {};
end

if nargin<3 || isempty(output_filename)
    output_filename = 'COUPP60XYZlookup.mat';
end

if nargin<2
    disp('You need 2 arguments when calling EricsQuickLookupTableMakerEtc');
    return
end


if size(resolutions,2)~=2
    disp('Resolutions input is weird in EricsQuickLookupTableMakerEtc');
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
%%% add torus geospecs here

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
        disp(sprintf('Watch out, geometry is out of bounds for air-side of window for cam %d',i_cam));
    end
    raylist(this_ix,:) = RefractionReflectionAtInterface(raylist(this_ix,:), normals, n_air, n_window);
    
    [ray_startingpoints(this_ix,:), normals, lt, or] = RayToPlane(ray_startingpoints(this_ix,:), raylist(this_ix,1:3), window_hydraulicsidepoint, window_normal);
    if any(isinf(lt) | isnan(lt) | (lt <= 0) | (or ~= -1))
        disp(sprintf('Watch out, geometry is out of bounds for hydraulic-side of window for cam %d',i_cam));
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

%%

if adams_arg
    
    %Get random points uniformly distributed in the jar volume up to z=85
    npts = 10000;
    pts = zeros(npts,3);
    for i = 1:npts
        pt = [-99 -99 -99];
        while pt == [-99 -99 -99]
            p = [31*rand() - 15.5, 31*rand() - 15.5, 95*rand() - 10];
            if p(1)^2 + p(2)^2 < r1(2)^2 && p(3) >= 0
               pt = p;
            elseif norm(p - cyl_axis*d(2)) < r3(2) && dot(cyl_axis, p) < z(2)
                pt = p;
            elseif dot(cyl_axis, p) < 0 && dot(cyl_axis, p) >= z(2)
                if p(1)^2 + p(2)^2 < r3(2)^2 - (abs(z(2)) + cyl_axis(3)*d(2))^2
                    pt = p;
                elseif (norm(p)^2 - (r1(2) - r2(2))^2 - r2(2)^2)^2 - 4*(r1(2) - r2(2))^2*(r2(2)^2-p(3)^2) < 0
                    pt = p;
                end
            end
        end
        pts(i,:) = pt;
    end
    
    
    numrays = size(raylist,1);

    % Scatter points as a variable gives the location of the ray at each point
    % it interacts with a surface
    scatter_points = zeros([max_scatters+2, numrays, 3]);
    scatter_points(:) = NaN;
    scatter_points(1,:,:) = ray_startingpoints;

    for n=1:length(raytracer_output)
        scatter_points(n+1,abs(raytracer_output(n).ray_index),:) = raytracer_output(n).intersection_point;
    end

    % The following plots the rays in the yz-plane
    plotrays = reshape(scatter_points,[],3);

    close all
    clf;
    set(gca,'fontsize',16);
    hold on
    axis equal

    SectionPlotter(surface_list, [0 0 0], [1 0 0], 10000);
    plot(plotrays(1:(numrays*(max_scatters+2)),2),plotrays(1:(numrays*(max_scatters+2)),3),'-or');
    scatter(pts(:,2), pts(:,3), 'g')
    ax = gca;
    ax.YLim = [-10,90];
    return

end
















%% make lookup table
scatter_points = zeros([size(raylist,1), 3, 2+max_scatters]) + inf;
directions = zeros([size(raylist,1), 3, 2+max_scatters]) + inf;
scatter_here = false(size(raylist,1), 2+max_scatters);


scatter_points(:,:,1) = ray_startingpoints;
directions(:,:,1) = raylist(:,1:3);
scatter_here(:,1) = true;

for s=1:length(raytracer_output)
    refracted_cut = raytracer_output(s).refracted_ray(:,7) > 0;
    scatter_points(raytracer_output(s).ray_index(refracted_cut),:,s+1) =  ...
        raytracer_output(s).intersection_point(refracted_cut,:);
    directions(raytracer_output(s).ray_index(refracted_cut),:,s+1) = ...
        raytracer_output(s).refracted_ray(refracted_cut,1:3);
    scatter_here(raytracer_output(s).ray_index(refracted_cut), s+1) = true;

    reflected_cut = ~refracted_cut & raytracer_output(s).reflected_ray(:,7) > 0;
    scatter_points(raytracer_output(s).ray_index(reflected_cut),:,s+1) =  ...
        raytracer_output(s).intersection_point(reflected_cut,:);
    directions(raytracer_output(s).ray_index(reflected_cut),:,s+1) = ...
        raytracer_output(s).reflected_ray(reflected_cut,1:3);
    scatter_here(raytracer_output(s).ray_index(reflected_cut), s+1) = true;
end

maxlengths = diff(scatter_points, 1, 3);
maxlengths = squeeze(sqrt(sum(maxlengths.^2,2)));

num_segments = sum(scatter_here, 2);

for c=1:n_cam
    pixel_list = struct('starting_point', [], 'direction', [], 'maxlength', []);
    pixel_list(prod(resolutions(c,:))) = pixel_list;

    for p=1:prod(resolutions(c,:))
        pixel_list(p).starting_point = squeeze( ...
            scatter_points(p+pixels_camstarts(c)-1,:,1:num_segments(p+pixels_camstarts(c)-1)))';
        if size(pixel_list(p).starting_point,2)==1
            pixel_list(p).starting_point = pixel_list(p).starting_point';
        end
        pixel_list(p).direction = squeeze( ...
            directions(p+pixels_camstarts(c)-1,:,1:num_segments(p+pixels_camstarts(c)-1)))';
        if size(pixel_list(p).direction,2)==1
            pixel_list(p).direction = pixel_list(p).direction';
        end
        pixel_list(p).maxlength = ...
            maxlengths( p+pixels_camstarts(c)-1, 1:num_segments(p+pixels_camstarts(c)-1) );
    end

    pixel_list = reshape(pixel_list, resolutions(c,:));

    switch c
        case 1
            cam0_pixels = pixel_list;
        case 2
            cam1_pixels = pixel_list;
        case 3
            cam2_pixels = pixel_list;
        case 4
            cam3_pixels = pixel_list;
    end
end


jar_iaxrad = jar_Iaxrad;
jar_irad = .5*jar_ID;
target_segment  = 3;

save(output_filename, 'cam0_pixels', 'cam1_pixels', 'cam2_pixels', 'cam3_pixels', 'target_segment', 'jar_irad', 'jar_iaxrad', '-mat');

%% make a few more useful output files, first an auxilliary .mat file (small) with some nice functions
jar_irad_fun = @(zpos)( ...
    (.5 .* jar_ID .* (zpos>=0)) + ...
    ((r1(2)-r2(2) + sqrt(r2(2).^2 - zpos.^2)) .* ((zpos<0) & (zpos>=z(2)))) + ...
    (sqrt(r3(2).^2 - (zpos-d(2)).^2) .* (zpos<z(2))));

dwall_horiz = @(rpos, zpos)(jar_irad_fun(zpos) - rpos);

dwall_fitfun = @(rpos, zpos, zwall)( ...
    sqrt((zpos-zwall).^2 + (rpos-jar_irad_fun(zwall)).^2) );

options = optimset('MaxFunEvals',1e3,'MaxIter',500,'GradObj','off','Hessian','off','Display','off');

dwall_rz = @(rpos, zpos)( ...
    dwall_fitfun(rpos, zpos, lsqnonlin(@(zwall)dwall_fitfun(rpos, zpos, zwall), zpos, d(2)-r3(2), max(0,zpos), options)) .* ...
    sign(jar_irad_fun(zpos)-rpos) );

z_fun = @(xyz_pos)(sum(xyz_pos .* repmat(cyl_axis, size(xyz_pos,1), 1), 2));
r_fun = @(xyz_pos)(sqrt(sum(xyz_pos.^2, 2) - z_fun(xyz_pos).^2));

dwall_xyz = @(xyz_pos)(dwall_rz(r_fun(xyz_pos), z_fun(xyz_pos)));

save([output_filename(1:end-4) '_aux.mat'], 'jar_irad_fun', 'hemi_rotmat', 'dwall_rz', 'dwall_horiz', 'dwall_xyz', '-mat');

%% Now a recon-style binary file
fid = fopen([output_filename(1:end-4) '_targetonly.bin'], 'wb');

fwrite(fid, hex2dec('01020304'), 'uint32');

cam_targetmatrix = zeros([4, resolutions(1,:), 3, 2]) + nan;

recon_quartzrays = false;

for ipix=1:resolutions(1,1)
    for jpix=1:resolutions(1,2)
        if size(cam0_pixels(ipix,jpix).starting_point, 1) == 5
            cam_targetmatrix(1,ipix,jpix,:,1) = cam0_pixels(ipix,jpix).starting_point(3,:);
            cam_targetmatrix(1,ipix,jpix,:,2) = cam0_pixels(ipix,jpix).starting_point(3,:) + cam0_pixels(ipix,jpix).direction(3,:) * cam0_pixels(ipix,jpix).maxlength(3);
        elseif recon_quartzrays && size(cam0_pixels(ipix,jpix).starting_point, 1) == 4
            cam_targetmatrix(1,ipix,jpix,:,1) = cam0_pixels(ipix,jpix).starting_point(2,:);
            cam_targetmatrix(1,ipix,jpix,:,2) = cam0_pixels(ipix,jpix).starting_point(3,:) + cam0_pixels(ipix,jpix).direction(3,:) * cam0_pixels(ipix,jpix).maxlength(3);
        elseif recon_quartzrays && size(cam0_pixels(ipix,jpix).starting_point, 1) == 3
            cam_targetmatrix(1,ipix,jpix,:,1) = cam0_pixels(ipix,jpix).starting_point(2,:);
            cam_targetmatrix(1,ipix,jpix,:,2) = cam0_pixels(ipix,jpix).starting_point(2,:) + cam0_pixels(ipix,jpix).direction(2,:) * cam0_pixels(ipix,jpix).maxlength(2);
        end
        
        if size(cam1_pixels(ipix,jpix).starting_point, 1) == 5
            cam_targetmatrix(2,ipix,jpix,:,1) = cam1_pixels(ipix,jpix).starting_point(3,:);
            cam_targetmatrix(2,ipix,jpix,:,2) = cam1_pixels(ipix,jpix).starting_point(3,:) + cam1_pixels(ipix,jpix).direction(3,:) * cam1_pixels(ipix,jpix).maxlength(3);
        elseif recon_quartzrays && size(cam1_pixels(ipix,jpix).starting_point, 1) == 4
            cam_targetmatrix(2,ipix,jpix,:,1) = cam1_pixels(ipix,jpix).starting_point(2,:);
            cam_targetmatrix(2,ipix,jpix,:,2) = cam1_pixels(ipix,jpix).starting_point(3,:) + cam1_pixels(ipix,jpix).direction(3,:) * cam1_pixels(ipix,jpix).maxlength(3);
        elseif recon_quartzrays && size(cam1_pixels(ipix,jpix).starting_point, 1) == 3
            cam_targetmatrix(2,ipix,jpix,:,1) = cam1_pixels(ipix,jpix).starting_point(2,:);
            cam_targetmatrix(2,ipix,jpix,:,2) = cam1_pixels(ipix,jpix).starting_point(2,:) + cam1_pixels(ipix,jpix).direction(2,:) * cam1_pixels(ipix,jpix).maxlength(2);
        end
        
        if size(cam2_pixels(ipix,jpix).starting_point, 1) == 5
            cam_targetmatrix(3,ipix,jpix,:,1) = cam2_pixels(ipix,jpix).starting_point(3,:);
            cam_targetmatrix(3,ipix,jpix,:,2) = cam2_pixels(ipix,jpix).starting_point(3,:) + cam2_pixels(ipix,jpix).direction(3,:) * cam2_pixels(ipix,jpix).maxlength(3);
        elseif recon_quartzrays && size(cam2_pixels(ipix,jpix).starting_point, 1) == 4
            cam_targetmatrix(3,ipix,jpix,:,1) = cam2_pixels(ipix,jpix).starting_point(2,:);
            cam_targetmatrix(3,ipix,jpix,:,2) = cam2_pixels(ipix,jpix).starting_point(3,:) + cam2_pixels(ipix,jpix).direction(3,:) * cam2_pixels(ipix,jpix).maxlength(3);
        elseif recon_quartzrays && size(cam2_pixels(ipix,jpix).starting_point, 1) == 3
            cam_targetmatrix(3,ipix,jpix,:,1) = cam2_pixels(ipix,jpix).starting_point(2,:);
            cam_targetmatrix(3,ipix,jpix,:,2) = cam2_pixels(ipix,jpix).starting_point(2,:) + cam2_pixels(ipix,jpix).direction(2,:) * cam2_pixels(ipix,jpix).maxlength(2);
        end
        
        if size(cam3_pixels(ipix,jpix).starting_point, 1) == 5
            cam_targetmatrix(4,ipix,jpix,:,1) = cam3_pixels(ipix,jpix).starting_point(3,:);
            cam_targetmatrix(4,ipix,jpix,:,2) = cam3_pixels(ipix,jpix).starting_point(3,:) + cam3_pixels(ipix,jpix).direction(3,:) * cam3_pixels(ipix,jpix).maxlength(3);
        elseif recon_quartzrays && size(cam3_pixels(ipix,jpix).starting_point, 1) == 4
            cam_targetmatrix(4,ipix,jpix,:,1) = cam3_pixels(ipix,jpix).starting_point(2,:);
            cam_targetmatrix(4,ipix,jpix,:,2) = cam3_pixels(ipix,jpix).starting_point(3,:) + cam3_pixels(ipix,jpix).direction(3,:) * cam3_pixels(ipix,jpix).maxlength(3);
        elseif recon_quartzrays && size(cam3_pixels(ipix,jpix).starting_point, 1) == 3
            cam_targetmatrix(4,ipix,jpix,:,1) = cam3_pixels(ipix,jpix).starting_point(2,:);
            cam_targetmatrix(4,ipix,jpix,:,2) = cam3_pixels(ipix,jpix).starting_point(2,:) + cam3_pixels(ipix,jpix).direction(2,:) * cam3_pixels(ipix,jpix).maxlength(2);
        end
    end
end
headerstring = sprintf('lookuptable;double;4,%d,%d,3,2;',resolutions(1,1), resolutions(1,2));

fwrite(fid, length(headerstring), 'uint16');
fwrite(fid, headerstring, 'char');
fwrite(fid, 0, 'int32');
fwrite(fid, cam_targetmatrix(:), 'double');
fclose(fid);
%% make pictures

if isempty(imagefiles)
    return
end

for i_cam = 1:n_cam

    figure(100+i_cam);
    clf
    pdata = imread(imagefiles{i_cam});
    if size(pdata,3)==1
        pdata = repmat(pdata,[1,1,3]);
    end
    image(pdata);
    hold on;
    
    pvalues = reshape(num_segments(pixels_camstarts(i_cam):pixels_camends(i_cam)),resolutions(i_cam,:));
    
    vert_boundaries = diff(pvalues,1,1)~=0;
    horz_boundaries = diff(pvalues,1,2)~=0;
    
    xd = (2:resolutions(i_cam,1)) - .5;
    ydl = (1:resolutions(i_cam,2)) - .5;
    ydu = (1:resolutions(i_cam,2)) + .5;

    yd = (2:resolutions(i_cam,2)) - .5;
    xdl = (1:resolutions(i_cam,1)) - .5;
    xdu = (1:resolutions(i_cam,1)) + .5;
    
    xd = reshape(repmat(xd(:),resolutions(i_cam,2),1), size(vert_boundaries));
    ydl = reshape(repmat(ydl(:)',resolutions(i_cam,1)-1,1), size(vert_boundaries));
    ydu = reshape(repmat(ydu(:)',resolutions(i_cam,1)-1,1), size(vert_boundaries));
    
    yd = reshape(repmat(yd(:)',resolutions(i_cam,1),1), size(horz_boundaries));
    xdl = reshape(repmat(xdl(:),resolutions(i_cam,2)-1,1), size(horz_boundaries));
    xdu = reshape(repmat(xdu(:),resolutions(i_cam,2)-1,1), size(horz_boundaries));
    
    xdata = [ reshape([reshape(xd(vert_boundaries),1,[]) ; reshape(xd(vert_boundaries),1,[]) ; repmat(NaN,1,sum(vert_boundaries(:)))] , [], 1) ;...
        reshape([reshape(xdl(horz_boundaries),1,[]) ; reshape(xdu(horz_boundaries),1,[]) ; repmat(NaN,1,sum(horz_boundaries(:)))], [], 1)];

    ydata = [ reshape([reshape(ydl(vert_boundaries),1,[]) ; reshape(ydu(vert_boundaries),1,[]) ; repmat(NaN,1,sum(vert_boundaries(:)))], [], 1) ; ...
        reshape([reshape(yd(horz_boundaries),1,[]) ; reshape(yd(horz_boundaries),1,[]) ; repmat(NaN,1,sum(horz_boundaries(:)))], [], 1)];
    
    
    plot(xdata, ydata, 'r');

end
