%% SUMMARY
% This code is to show the how the transducers are placed relative to the
% reconstructed bone from CT scan.
% 
% IMPORTANT NOTE: 
% Please specify the correct folder name or file name for these two
% variables, path_measurement and filename_bone.
% for path_measurement, you can choose 'ct1' or 'ct2'
% for filename_bone, you can choose 'CT1_tibia_smooth-remesh.stl' or 'CT2_tibia_raw.stl'

clc; clear; clf;

addpath('functions\sphereFit\');
addpath('functions\affine_fit\');
addpath('functions\TriangleRayIntersection\');

% specify the path measurement
path_data        = 'data_organized';
% path_measurement = strcat(path_data, filesep, 'CT', filesep, 'ct1');
path_measurement = strcat(path_data, filesep, 'CT', filesep, 'ct2');

% read fiducials stl data -------------------------------------------------
path_fiducials = strcat(path_measurement, filesep, 'fiducials');
filelist       = dir(strcat(path_fiducials, filesep, '*.stl'));
filenames      = strcat(path_fiducials, filesep, vertcat(filelist.name));

fiducialgroup_stls = {};
probefiducial_idxs = [20, 21, 24, 25, 26, 27, 28, 29, 30];
for i=1:size(filenames,1)
    fiducialgroup_stls{probefiducial_idxs(i)} = stlread(filenames(i,:));
end
fiducialgroup_stls = fiducialgroup_stls';

% read probe data ---------------------------------------------------------
% path_probes     = 'probes\probes_v2';
path_probes     = 'probes';
fullpath_probes = strcat(path_measurement, filesep, path_probes);
filelist        = dir(strcat(fullpath_probes, filesep, '*.stl'));
filenames       = strcat(fullpath_probes, filesep, vertcat(filelist.name));

probes = {};
for i=1:size(filenames,1)
    probes{probefiducial_idxs(i)} = stlread(filenames(i,:));
end
probes = probes';

% read bone stl data ------------------------------------------------------
% filename_bone = 'CT1_tibia_smooth-remesh.stl';
filename_bone = 'CT2_tibia_raw.stl';
path_bone     = strcat(path_measurement, filesep, 'bone');
fullpath_bone = strcat(path_bone, filesep, filename_bone);
bone_stl      = stlread(fullpath_bone);

% read marker-transducer calibration data ---------------------------------
filename_transducercalib = 'transducercalibration.csv';
path_transducercalib     = strcat(path_data, filesep, 'CT', filesep, filename_transducercalib); 
calibration_markerstick  = readmatrix(path_transducercalib);

% program modes
display            = true;
probeRigidBodymode = 'probe';

%% Main Loop

% initialize the structures, they will store our useful information:
% 1) marker     : for individual fiducials
% 2) markerstick: for group of fiducials
% 3) depth      : for the detected depth for each individual transducer
% 4) transducer : as the name, the transducer, which also includes every
%                 information mentioned above
marker      = struct('pointcloud', [], 'center', [], 'radius', 0);
markerstick = struct('number', 0, 'origin', [], 'basisvector', [], 'rigidbody', [], 'error', 0, 'markers', struct(marker) );
depth       = struct('type', '', 'mean', 0, 'max', 0, 'min', 0, 'points', []);
transducer  = struct('number', 0, 'origin', [], 'basisvector', [], 'rigidbody', [], 'depth', struct(depth), 'markerstick', struct(markerstick));

% this variable is for collecting all the transducers, but i need to
% initialize it first
markersticks  = markerstick;
transducers   = transducer;

% loop for all the fiducials/markers
for probefiducial_idx = probefiducial_idxs

    % take the current fiducials
    fiducialgroup_current = fiducialgroup_stls{probefiducial_idx};
    probe_current         = probes{probefiducial_idx};

    %=====================================================================%
    % Setting up the rigid body of the markers                            %
    %=====================================================================%
    
    % segment the fiducial group to get the individual
    labels      = pcsegdist( pointCloud(fiducialgroup_current.Points), 0.5);
    label_index = unique(labels);
    tmp_markers   = [];
    for index = 1:length(label_index)
        tmp_marker.pointcloud   = fiducialgroup_current.Points(labels==index, :);
        [tmp_marker.center, tmp_marker.radius] = sphereFit(tmp_marker.pointcloud);
        tmp_markers = [tmp_markers; tmp_marker];
    end

    % rearrange marker
    marker_points = vertcat(tmp_markers.center);
    marker_order  = rearrangeMarker(marker_points);
    markers       = tmp_markers(marker_order);

    % Estimate the plane (the normal)
    marker_points = vertcat(markers.center);
    [plane_center, plane_basisvector, error] = estimateMarkerRigidBody(marker_points);

    % collect the data    
    markerstick.number         = probefiducial_idx;
    markerstick.origin         = plane_center;
    markerstick.basisvector    = plane_basisvector;
    markerstick.rigidbody      = [plane_basisvector, plane_center'; 0 0 0 1];
    markerstick.error          = error;
    markerstick.markers        = markers;
    markersticks(probefiducial_idx) = markerstick;    

    %=====================================================================%
    % Setting up the rigid body of the transducer                         %
    %=====================================================================%
    
    % This block of code is using fiducials + prior knowledge of holder
    % structures (presented in calibration matrix) to calculate the
    % rigidbody of the transducers
    if (strcmp(probeRigidBodymode, 'marker'))

        % obtain the transformations
        T_global_markerstick    = markerstick.rigidbody;
        T_markerstick_usttip    = reshape(calibration_markerstick(probefiducial_idx,:), [4,4])';
        T_global_usttip         = T_global_markerstick * T_markerstick_usttip;
    
        % get the informations
        transducer.number       = probefiducial_idx;
        transducer.origin       = T_global_usttip(1:3, 4);
        transducer.basisvector  = T_global_usttip(1:3, 1:3);
        transducer.rigidbody    = T_global_usttip;
        transducer.markerstick  = markerstick;
        transducers(probefiducial_idx) = transducer;

    % This block of code is using the segmented 3D mesh of the probe to
    % calculate its rigid body transformation
    elseif (strcmp(probeRigidBodymode, 'probe'))

        % estimate rigid body transformation of the probe
        probe_current_rigidbody = estimateProbeRigidBody(probe_current, markerstick, 'pca');

        % get the informations
        transducer.number       = probefiducial_idx;
        transducer.origin       = probe_current_rigidbody(1:3, 4);
        transducer.basisvector  = probe_current_rigidbody(1:3, 1:3);
        transducer.rigidbody    = probe_current_rigidbody;
        transducer.markerstick  = markerstick;
        transducers(probefiducial_idx) = transducer;

    end

    % display
    if(display)

        % markers
        colorcode = ['r', 'g', 'b', 'm'];
        for i=1:size(markers,1)
            % before arrangement
            % plot3( tmp_markers(i).pointcloud(:,1), ...
            %        tmp_markers(i).pointcloud(:,2), ...
            %        tmp_markers(i).pointcloud(:,3), ...
            %        '.', 'Color', colorcode(i), 'MarkerSize', 0.01, 'Tag', 'plot_pointcloud');
            % after arrangement
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

        % markers rigid body
        quiverscale = 10;
        plot3(plane_center(1), plane_center(2), plane_center(3),'ok', 'MarkerFaceColor', 'black');
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

        % transducer rigid body
        plot3(transducer.origin(1), transducer.origin(2), transducer.origin(3),'ok', 'MarkerFaceColor', 'black');
        quiver3( transducer.origin(1), ...
                 transducer.origin(2), ...
                 transducer.origin(3), ...
                 transducer.basisvector(1,1)*quiverscale, ...
                 transducer.basisvector(2,1)*quiverscale, ...
                 transducer.basisvector(3,1)*quiverscale, ...
                 'r', 'linewidth',2);
        quiver3( transducer.origin(1), ...
                 transducer.origin(2), ...
                 transducer.origin(3), ...
                 transducer.basisvector(1,2)*quiverscale, ...
                 transducer.basisvector(2,2)*quiverscale, ...
                 transducer.basisvector(3,2)*quiverscale, ...
                 'g', 'linewidth',2);
        quiver3( transducer.origin(1), ...
                 transducer.origin(2), ...
                 transducer.origin(3), ...
                 transducer.basisvector(1,3)*quiverscale, ...
                 transducer.basisvector(2,3)*quiverscale, ...
                 transducer.basisvector(3,3)*quiverscale, ...
                 'b', 'linewidth',2);

        % probes
        trimesh(probe_current, 'FaceColor', '#34495e', 'FaceAlpha', 0.5, 'LineStyle', 'none');

    end

end

trimesh(bone_stl, 'FaceColor', '#fad390', 'FaceAlpha', 0.5, 'LineStyle', 'none');
lighting gouraud;
camlight('headlight')

% save('data/transducers2_probebased.mat', 'transducers');

















