% function [rays pixels] = GenerateRaysFromCamera(pinhole_position, resolution, ...
%     pixel_pitch, pixel_center, focal_length, pitch, yaw, roll, radial_distortion)
%
% 12/23/09, CED

function [ray_direction, pixels] = GenerateRaysFromCameraCommented(resolution, pixel_pitch, ...
    pixel_center, focal_length, pitch, yaw, roll, radial_distortion, lens_type)

%% set default
% All this is doing is entering any arguments that have been left out
if nargin<9 || isempty(lens_type)
    lens_type = 'tan';
end

if nargin<8 || isempty(radial_distortion)
    radial_distortion = [];
end

if nargin<7 || isempty(roll)
    roll = 0;
end

if nargin<6 || isempty(yaw)
    yaw = 0;
end

if nargin<5 || isempty(pitch)
    pitch = 0;
end

% Checks to make sure everything that is required is present
if nargin<4 || numel(focal_length)~=1 || numel(yaw)~=1 || ...
        numel(pitch)~=1 || numel(roll)~=1 || numel(radial_distortion)>2 || ...
        numel(pixel_center)~=2 || isempty(pixel_pitch) || numel(pixel_pitch)>2 || ...
        numel(resolution)~=2
    disp('improper input to GenerateRaysFromCamera');
    return
end

% Sets up pixel pitch as two values instead of one
if numel(pixel_pitch)==1
    pixel_pitch = [0 0] + pixel_pitch;
end

%% generate rays in camera-frame
% in camera frame, forward is +y, +x is +i, and +z is -j

% This sets up the image as it will be seen on a camera
% The z is upside down as it is seen on the camera film
% it gives a resolution(1) x resolution(2) array of pixels
i_pix = repmat((1:resolution(1))',1,resolution(2));
j_pix = repmat((1:resolution(2)),resolution(1),1);
pixels = [i_pix(:) j_pix(:)];

% pixel_location is a Xx3 array where X is the number of pixels in the
% image, i.e. resolution(1)xresolution(2)
% Column 1 centers the i_pix such that 0 is at the x-value of pixel_center and
% fills the row with all of the new i_pix values, these values are the
% actual distances away from the x-center each pixel is
% Column 2 has resolution(1) number of entries, all of which are -focal_length
% Column 3 does the same as Row 1 for the z-values
pixel_location = [(pixel_center(1)-i_pix(:)).*pixel_pitch(1), ...
    zeros(numel(i_pix),1)-focal_length, ...
    -(pixel_center(2)-j_pix(:)).*pixel_pitch(2)];

% pixel_d2 is a Xx1 column vector where X is the total number of pixels
% Each row gives the square distance of the pixel to the center of the
% image
pixel_d2 = sum(pixel_location(:,[1 3]).^2,2);

% The next set of variables are all used in calculation the effective focal
% length
% x is an Xx1 column vector, X is the total number of pixels, all values
% are the radial distortion value
x = repmat(radial_distortion(:)',length(pixel_d2),1);
% y is an Xx1 column vector, X is the total number of pixels, all values
% are 1 when radial distortion is present, 0 if not
y = repmat(1:length(radial_distortion),length(pixel_d2),1);
% z is a Xx1 column vector, X is the total number of pixels, each value is
% a pixel distance over the focal length, all squared
z = repmat(focal_length^-2 * pixel_d2,1,length(radial_distortion));
% This is either still z, or all 1's in replace of z
w = (z.^y);
% effective_f is the focal length for each of the individual pixels when radial distortion occurs, it is
% an Xx1 column vector where X is the number of pixels.
effective_f = focal_length * (1 + (sum((x .* w),2)));

% Depending on the lens type this switch chooses a different formula to
% calculate the angle at which a parallel ray from the pixel will exit the
% lens, this is an Xx1 column vector, where X is the number of pixels
switch lens_type
    case 'theta'
        theta = sqrt(pixel_d2)./effective_f;
    case 'sin'
        theta = asin(sqrt(pixel_d2)./effective_f);
    case 'tan'
        theta = atan(sqrt(pixel_d2)./effective_f);
    otherwise
        theta = atan(sqrt(pixel_d2)./effective_f);
end

% phi is again an Xx1 column vector, where X is the number of pixels, this
% represents the three dimensional angle since the image is produced on a
% square not a line
phi = atan2(-pixel_location(:,3),-pixel_location(:,1));

% ray_direction gives the direction of the parallel ray from each pixel as
% the x, y, z components of a unit vector
ray_direction = [sin(theta).*cos(phi), cos(theta), sin(theta).*sin(phi)];
%ray_direction = zeros(size(ray_direction));
%ray_direction(:,2)=ones(size(ray_direction(:,2)));

%% rotate camera
% This adjusts the x, y, z components of the direction vector based on the
% yaw, pitch, and roll of the camera which should be input as radians,
% note: roll is around y-axis, pitch is around x-axis, yaw is around
% z-axis
M1 = [cos(yaw) -sin(yaw) 0 ; sin(yaw) cos(yaw) 0 ; 0 0 1];
M2 = [1 0 0 ; 0 cos(pitch) -sin(pitch) ; 0 sin(pitch) cos(pitch)];
M3 = [cos(roll) 0 sin(roll) ; 0 1 0 ; -sin(roll) 0 cos(roll)];

M = M1*M2*M3;

ray_direction = (M * (ray_direction'))';