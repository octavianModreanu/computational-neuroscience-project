%% Stimulus and visual-field grid
% Loads the bar stimulus used for pRF mapping and sets up a coordinate
% grid so every stimulus pixel has a position in the visual field.
% All positions are in degrees of visual angle, measured from fixation
% (the point the participant looks at, placed at the centre of the image).

% The stimulus is a movie of a bar at different positions and orientations.
% Each column is one frame: a 150x150 binary image (1 = bar, 0 = background)
% flattened into a single column of 22500 values.

load('stimulus.mat');                  % stimulus: 22500 x 304 (pixels x frames)


p.n_px       = 150;                    % image is 150 x 150 pixels
p.fov_radius = 10;                     % half-width of the visual field in degrees,
% so the image spans -10 to +10 deg.
p.frame_dt   = 2;                      % seconds per frame = TR of the scan

R = p.fov_radius;


% Coordinate grid: X(i,j) and Y(i,j) give the visual-field position of
% pixel (i,j) in degrees. Used to draw each pRF as a 2D Gaussian image.
%   x runs left to right  (-R at column 1, +R at column 150)
%   y runs top to bottom  (+R at row 1,    -R at row 150)
[X, Y] = meshgrid(linspace(-R, R, p.n_px),linspace( R, -R, p.n_px));     

A        = stimulus;                             
n_frames = size(A, 2);                           % 304
frame_t  = (0:n_frames-1) * p.frame_dt;          % 0, 2, 4, ... 606 s
T        = n_frames * p.frame_dt;                % total run length: 608 s