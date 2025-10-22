% This script is only for generating probe_toIntersectionArea_angle

clc; clear;

% let's add function that we need
addpath('functions\affine_fit\');

% path to our data
path_data        = 'data';
path_measurement = strcat(path_data, filesep, 'ct2');

% load ground truth data (tube based)
filename_transducermat = 'transducers2_probebased_tubedepth.mat';
fullpath_transducermat = strcat(path_measurement, filesep, filename_transducermat);
load(fullpath_transducermat);
transducers_TD = transducers;

% load ground truth data (line based)
filename_transducermat = 'transducers2_probebased_linedepth.mat';
fullpath_transducermat = strcat(path_measurement, filesep, filename_transducermat);
load(fullpath_transducermat);
transducers_LD = transducers;

% read bone stl data ------------------------------------------------------
% filename_bone = 'Tibia_smooth-remesh.stl';
filename_bone = 'CT2_tibia_raw.stl';
path_bone     = strcat(path_measurement, filesep, 'bone');
fullpath_bone = strcat(path_bone, filesep, filename_bone);
bone_stl      = stlread(fullpath_bone);

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

% flag for display
display = true;

%%

% create figure
figure1 = figure('Name', 'Intersection Angle');

% variable to store (as the name suggests) an angle from the probe to the
% intersection area. I hypotheszied there is a correlation between the
% angle and the width of the peak in the ultrasound signal.
probe_toIntersectionArea_angle = zeros(size(transducers_LD))';
probe_toIntersectionArea_sidewayBeamLength = zeros(size(transducers_LD))';

% loop for selected transducers
for transducer_idx = probefiducial_idxs

    % get both transducers data from line method and tube method
    currenttransducer_LD = transducers_LD(transducer_idx);
    currenttransducer_TD = transducers_TD(transducer_idx);
    
    % get the middle point of the intersection
    if(~isempty(currenttransducer_LD.depth.points))
        bone_toTransducer_lineIntersection = currenttransducer_LD.depth.points;
        bone_toTransducer_depth = currenttransducer_LD.depth.mean;

    elseif((~isempty(currenttransducer_TD.depth.points)))
        bone_toTransducer_lineIntersection = mean(currenttransducer_TD.depth.points, 1);
        bone_toTransducer_depth = currenttransducer_TD.depth.mean;

    % if intersection area is not found, then fuck this probe
    else
        continue;
    end

    % get the intersection area
    intersectionArea_idx    = knnsearch(bone_stl.Points, bone_toTransducer_lineIntersection, 'K', 300);
    intersectionArea_points = bone_stl.Points(intersectionArea_idx, :);

    % fit the plane, and obtain only the plane normal
    [intersectionArea_normal, ~, ~] = affine_fit(intersectionArea_points);
    % correct the normal so that it is in the direction of probe normal
    sign_corrector          = sign(dot(currenttransducer_LD.basisvector(:,3), intersectionArea_normal'));
    intersectionArea_normal = sign_corrector*intersectionArea_normal;

    % compute angle
    u         = intersectionArea_normal;
    v         = currenttransducer_LD.basisvector(:,3);
    costheta  = dot(u,v)/(norm(u)*norm(v));
    theta_deg = acosd(costheta);

    % compute sideway beam distance
    sidewayBeam_length = tan(deg2rad(theta_deg))*bone_toTransducer_depth;

    % store the angle
    probe_toIntersectionArea_angle(transducer_idx) = theta_deg;
    probe_toIntersectionArea_sidewayBeamLength(transducer_idx) = sidewayBeam_length;

    % turn this on if you want to display
    if (display)

        % display the intersection area
        plot3(intersectionArea_points(:,1), intersectionArea_points(:,2), intersectionArea_points(:, 3), '.c');
        hold on; grid on; axis equal;

        % display probes
        trimesh(probes{transducer_idx}, 'FaceColor', '#34495e', 'FaceAlpha', 0.5, 'LineStyle', 'none');

        % display ultrasound artificial beam
        if(~isempty(currenttransducer_TD.depth.points))
            plot3(currenttransducer_TD.depth.points(:,1), currenttransducer_TD.depth.points(:,2), currenttransducer_TD.depth.points(:,3), '.r', 'MarkerSize', 4);
        end

        % display the vectors
        quiverscale = 10;
        plot3(bone_toTransducer_lineIntersection(1), bone_toTransducer_lineIntersection(2), bone_toTransducer_lineIntersection(3),'ok', 'MarkerFaceColor', 'black');
        quiver3( bone_toTransducer_lineIntersection(1), ...
                 bone_toTransducer_lineIntersection(2), ...
                 bone_toTransducer_lineIntersection(3), ...
                 intersectionArea_normal(1)*quiverscale, ...
                 intersectionArea_normal(2)*quiverscale, ...
                 intersectionArea_normal(3)*quiverscale, ...
                 'm', 'linewidth',2);
        quiver3( bone_toTransducer_lineIntersection(1), ...
                 bone_toTransducer_lineIntersection(2), ...
                 bone_toTransducer_lineIntersection(3), ...
                 currenttransducer_LD.basisvector(1,3)*quiverscale, ...
                 currenttransducer_LD.basisvector(2,3)*quiverscale, ...
                 currenttransducer_LD.basisvector(3,3)*quiverscale, ...
                 'b', 'linewidth',2);

    % end display loop
    end

% end transducer loop
end


% display the bone
trimesh(bone_stl, 'FaceColor', '#fad390', 'FaceAlpha', 0.5, 'LineStyle', 'none');
lighting gouraud;
camlight('headlight');
zlim([-700, -625]);
% zlim([-800, -950]);
% zlim([-950, -850]);
% ylim([-80 0]);
% xlim([-50 25])