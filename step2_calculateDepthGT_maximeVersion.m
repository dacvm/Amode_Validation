%% SUMMARY
% This code is to calculate the ground truth of the depth of each
% individual ultrasound transducer. You can choose variable beamType to be 
% 'line' or 'tube'
% - This script required the MAT file
%   'transducers2_probebased.mat' to be generated from
%   'step1_generateVirtualSetup.m'
% - This script requires geom3d, which can be installed from here
%   https://nl.mathworks.com/matlabcentral/fileexchange/24484-geom3d
% 
% IMPORTANT NOTE: 
% Please specify the correct folder name or file name for these two
% variables, path_measurement and filename_bone.
% for path_measurement, you can choose 'ct1' or 'ct2'
% for filename_bone, you can choose 'CT1_tibia_smooth-remesh.stl' or 'CT2_tibia_raw.stl'


clc; clear; clf;

% path to our data
path_data        = 'data_organized';
% path_measurement = strcat(path_data, filesep, 'CT', filesep, 'ct1');
path_measurement = strcat(path_data, filesep, 'CT', filesep, 'ct2');

% read bone stl data ------------------------------------------------------
% filename_bone = 'CT1_tibia_smooth-remesh.stl';
filename_bone = 'CT2_tibia_raw.stl';
path_bone     = strcat(path_measurement, filesep, 'bone');
fullpath_bone = strcat(path_bone, filesep, filename_bone);
bone_stl      = stlread(fullpath_bone);

% read transducer mat file ------------------------------------------------
filename_transducermat = 'transducers2_probebased.mat';
fullpath_transducermat = strcat(path_measurement, filesep, filename_transducermat);
load(fullpath_transducermat);

% read probe data ---------------------------------------------------------
% path_probes        = 'probes\probe_v2';
path_probes        = 'probes';
fullpath_probes    = strcat(path_measurement, filesep, path_probes);
filelist           = dir(strcat(fullpath_probes, filesep, '*.stl'));
filenames          = strcat(fullpath_probes, filesep, vertcat(filelist.name));
probes             = {};
probefiducial_idxs = [20, 21, 24, 25, 26, 27, 28, 29, 30];
for i=1:size(filenames,1)
    probes{probefiducial_idxs(i)} = stlread(filenames(i,:));
end
probes = probes';

% program mode
display  = true;
beamType = 'tube';

%% Main Loop

% 3 vertices that makes the faces. This variable is required by 
% TriangleRayIntersection() function.
vert1 = bone_stl.Points(bone_stl.ConnectivityList(:,1),:);
vert2 = bone_stl.Points(bone_stl.ConnectivityList(:,2),:);
vert3 = bone_stl.Points(bone_stl.ConnectivityList(:,3),:);

% origin points, if we are not using single points
radius = 3;
step   = 1;
x      = -radius : step : radius;
[X, Y] = meshgrid(x, x);
xy     = [X(:), Y(:)];
xy_circle = xy( sqrt(sum(xy.^2, 2)) <= radius, :);

% prepare a struct to store our usefull information about the depth for
% each probes
% depth = struct('type', '', 'mean', [], 'max', [], 'min', [], 'points', []);

for transducer_idx = probefiducial_idxs

    % take the current fiducials
    transducer = transducers(transducer_idx);

    % there are several type of the beam that we use here:
    % 1) line: A hypothetical line that span from the origin point of
    %          transducers to the bone mesh.
    % 2) tube: A hypotethical tube that span from the entire surface of the
    %          transducer to the bone mesh. For simplicity, we define the
    %          origin point as a circle that lie on the transducer surface.
    %          So basically, it is a line but they are multiple.

    % if the user specify line
    if(strcmp(beamType, 'line'))

        % origin point is assumed to be the center of the transducer
        origins    = transducer.origin';
        % by convention, positive direction is away from the bone
        direction = -transducer.basisvector(:,3);

        % prepare some variables to store the intersection point and
        % intersection distances. actually, since in the "line" mode, we
        % don't need to initialize the variable, you can directly assign
        % the value to a variable when you find the intersection. but, i
        % just wanted to be consitent with "tube" mode. There, since i use
        % parallel for, it is required to initialize the variables first.
        lineplane_intersectpoints = [];
        lineplane_intersectdists  = [];

        % create the 3d line from the markerstick in the direction of the bone mesh
        % using x0, y0, z0 for the point and dx, dy, dz for the direction. The direction
        % is negative with respect to the plane basisvector, and we are only interested 
        % in the Z direction. The line is represented in parametric form : [x0 y0 z0 dx dy dz]
        line = [origins(1), origins(2), origins(3), direction(1), direction(2), direction(3)];  

        % find the intersection points between the line and the bone mesh, using the 3d
        % line, the vertices (= stl points) and faces (= stl connectivity
        % list).
        tmp_inters = intersectLineMesh3d(line, bone_stl.Points, bone_stl.ConnectivityList); % intersection between line and the bone mesh
        
        % find the intersection and ground truth depth
        % when the line and bone mesh do have an intersection 
        if ~isempty(tmp_inters) 
            % in case of multiple intersections (i.e. one intersection when the 
            % line enters the mesh and one when the line leaves the mesh),
            % calculate the distances for both intersections and take the smallest
            % distance as the distance we are interested in
            distancetointers = zeros(1,length(tmp_inters(:,1)));
            for i = 1:length(tmp_inters(:,1))
                % calculate the distance between the intersection point and transducer tip (origin),
                % for both intersections
                distancetointers(i) = distancePoints3d(tmp_inters(i,:),origins); % distancePoints3d(P1, P2) returns distance between points P1 and P2
            end
            % get distance ground truth data, using the smallest distance, because
            % intersections further away do not represent the bone surface 
            % that faces the ultrasound transducer 
            [lineplane_intersectdists, ind] = min(distancetointers); 
            % collect intersection closest to transducer origin in INTERS
            lineplane_intersectpoints = (tmp_inters(ind,:));
            
        % when the line and mesh do not have an intersection (probe is not facing
        % the tibia)
        else
           lineplane_intersectdists  = [];
           % define intersection point as an empty array
           lineplane_intersectpoints = [];
       end    

        % specify the type of the beam
        transducer.depth.type = 'line';
    
    % if the user specify tube
    elseif(strcmp(beamType, 'tube'))

        % Instead of one point, origin is represented as multiple points on
        % the probe beam surface. We do this by making a circle point grid
        % in global frame then transform it to transducer rigid body.
        origins    = homogeneous2cartesian( transducer.rigidbody * cartesian2homogeneous( [xy_circle, zeros(size(xy_circle, 1), 1)] ) );
        % by convention, positive direction is away from the bone
        direction = -transducer.basisvector(:,3);

        % prepare some variables, this is required if are using parfor
        lineplane_intersectpoints = zeros(size(origins));
        lineplane_intersectdists  = zeros(size(origins, 1), 1);

        tic;
        parfor i=1:length(origins)
            % for each point from the circle point grid, the the
            % intersection line to the bone mesh, represented in parametric form : [x0 y0 z0 dx dy dz]
            line = [origins(i,1), origins(i,2), origins(i,3), direction(1), direction(2), direction(3)];
            
            % find the intersection points between the line and the bone mesh, using the 3d
            % line, the vertices (= stl points) and faces (= stl connectivity
            % list).
            tmp_inters = intersectLineMesh3d(line, bone_stl.Points, bone_stl.ConnectivityList); % intersection between line and the bone mesh

            % find the intersection and ground truth depth
            % when the line and bone mesh do have an intersection 
            if ~isempty(tmp_inters)

                % in case of multiple intersections (i.e. one intersection when the 
                % line enters the mesh and one when the line leaves the mesh),
                % calculate the distances for both intersections and take the smallest
                % distance as the distance we are interested in                
                distancetointers = zeros(1,length(tmp_inters(:,1)));
                for j = 1:size(tmp_inters,1)
                    % calculate the distance between the intersection point and transducer tip (origin),
                    % for both intersections
                    distancetointers(j) = distancePoints3d(tmp_inters(j,:), origins(i,:)); % distancePoints3d(P1, P2) returns distance between points P1 and P2
                end

                % get distance ground truth data, using the smallest distance, because
                % intersections further away do not represent the bone surface 
                % that faces the ultrasound transducer 
                [lineplane_intersectdist, ind] = min(distancetointers); 
                % collect intersection closest to transducer origin in INTERS
                lineplane_intersectpoint = (tmp_inters(ind,:));

                % store the data to our variable that is declared outside
                % the parfor loop
                lineplane_intersectpoints(i,:) = lineplane_intersectpoint;
                lineplane_intersectdists(i)    = lineplane_intersectdist;
            end

        end
        toc;

        % If a point doesn't have any intersection, it will return nothing,
        % and since we avoid storing nothing to a variable (which will give 
        % us an error), the initial value (zero) is preserved. We can 
        % actually ignore them but it is annoying when we want to display 
        % the result. There will be points at zero coordinate, so let's
        % delete the delete rows with zero
        lineplane_intersectpoints(~any(lineplane_intersectpoints,2), :) = [];
        lineplane_intersectdists(lineplane_intersectdists==0) = [];

        % specify the type of the beam
        transducer.depth.type = 'tube';

    end

    % store some usefull information about the depth, and only store if 
    % there is an intersection points detected
    if(~isempty(lineplane_intersectdists))
        transducer.depth.mean   = mean(lineplane_intersectdists);
        transducer.depth.max    = max(lineplane_intersectdists);
        transducer.depth.min    = min(lineplane_intersectdists);
        transducer.depth.points = lineplane_intersectpoints;
    end
    % update the value of transducer
    transducers(transducer_idx) = transducer;

    % display
    if(display)

        % display markers
        colorcode = ['r', 'g', 'b', 'm'];
        for i=1:size(transducer.markerstick.markers,1)
            plot3( transducer.markerstick.markers(i).center(1), ...
                   transducer.markerstick.markers(i).center(2), ...
                   transducer.markerstick.markers(i).center(3), ...
                   'o', 'Color', colorcode(i), 'MarkerSize', 5, 'MarkerFaceColor', colorcode(i), 'Tag', 'plot_center');
            hold on;
        end
        grid on; axis equal;
        xlabel('X');
        ylabel('Y');
        zlabel('Z');

        % display markers rigid body
        quiverscale = 10;
        plot3(transducer.markerstick.origin(1), transducer.markerstick.origin(2), transducer.markerstick.origin(3),'ok', 'MarkerFaceColor', 'black');
        quiver3( transducer.markerstick.origin(1), ...
                 transducer.markerstick.origin(2), ...
                 transducer.markerstick.origin(3), ...
                 transducer.markerstick.basisvector(1,1)*quiverscale, ...
                 transducer.markerstick.basisvector(2,1)*quiverscale, ...
                 transducer.markerstick.basisvector(3,1)*quiverscale, ...
                 'r', 'linewidth',2);
        quiver3( transducer.markerstick.origin(1), ...
                 transducer.markerstick.origin(2), ...
                 transducer.markerstick.origin(3), ...
                 transducer.markerstick.basisvector(1,2)*quiverscale, ...
                 transducer.markerstick.basisvector(2,2)*quiverscale, ...
                 transducer.markerstick.basisvector(3,2)*quiverscale, ...
                 'g', 'linewidth',2);
        quiver3( transducer.markerstick.origin(1), ...
                 transducer.markerstick.origin(2), ...
                 transducer.markerstick.origin(3), ...
                 transducer.markerstick.basisvector(1,3)*quiverscale, ...
                 transducer.markerstick.basisvector(2,3)*quiverscale, ...
                 transducer.markerstick.basisvector(3,3)*quiverscale, ...
                 'b', 'linewidth',2);

        % display probes
        trimesh(probes{transducer_idx}, 'FaceColor', '#34495e', 'FaceAlpha', 0.5, 'LineStyle', 'none');

        % diplay ultrasound "beam"
        if(strcmp(beamType, 'line'))

            % transducer rigid body
            plot3(origins(1), origins(2), origins(3),'ok', 'MarkerFaceColor', 'k');
            % ultrasound artificial beam
            if(~isempty(lineplane_intersectpoints))
                beamline = [origins; lineplane_intersectpoints];
                plot3(beamline(:,1), beamline(:,2), beamline(:,3), '-r', 'LineWidth', 2);
                plot3(lineplane_intersectpoints(1), lineplane_intersectpoints(2), lineplane_intersectpoints(3), 'or', 'MarkerFaceColor', 'r');
            end

        elseif(strcmp(beamType, 'tube'))

            % transducer rigid body
            plot3(origins(:,1), origins(:,2), origins(:,3),'.k', 'MarkerSize', 4);
            % ultrasound artificial beam
            if(~isempty(lineplane_intersectpoints))
                % beamline = [origins; lineplane_intersectpoints];
                % plot3(beamline(:,1), beamline(:,2), beamline(:,3), '-r', 'LineWidth', 2);
                plot3(lineplane_intersectpoints(:,1), lineplane_intersectpoints(:,2), lineplane_intersectpoints(:,3), '.r', 'MarkerSize', 4);
            end

         % end beam
        end

    % end display
    end
% end for loop
end

trimesh(bone_stl, 'FaceColor', '#fad390', 'FaceAlpha', 0.5, 'LineStyle', 'none');
lighting gouraud;
camlight('headlight');

% save('data/transducers2_probebased_linedepth.mat', 'transducers');