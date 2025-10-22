function probe_current_rigidbody = estimateProbeRigidBody(probe, markerstick, method)
% The purpose of this function is to estimate the rigid body transformation
% for the probe. Cylinder fitting (registration) is the most logical way to
% do it, but it is quite tricky since the result of the segmentation is not
% perfectly a cylinder. The front side somehow has a buldge, the rear side
% somehow concave. So, when i wrote this function, first i avoid
% registration method, and i indeed have another idea to estimate the rigid
% body, by using PCA.

if strcmp(method, 'pca')
    % This method is using the PCA. This is because of the assumption
    % that the probe mesh is almost tubular, so that the main component
    % of the PCA is good enough to determine to direction of the probe
    probe_current_PCAbasisvector = pca(probe.Points);
    probe_current_PCAcentroid    = mean(probe.Points, 1);
    
    % Prepare the variables needed to calculate the intersection 
    % between the probe's centroid to the probe surfaces in the
    % direction of the main principal component
    origins    = probe_current_PCAcentroid;
    direction  = probe_current_PCAbasisvector(:,1);
    % if the direction of the main principal component the same as
    % markerstick, negate it 
    if(dot(direction, markerstick.basisvector(:,3)) > 0)
        direction = -direction;
    end
    
    % 3 vertices that makes the faces. This variable is required by 
    % TriangleRayIntersection() function.
    vert1 = probe.Points(probe.ConnectivityList(:,1),:);
    vert2 = probe.Points(probe.ConnectivityList(:,2),:);
    vert3 = probe.Points(probe.ConnectivityList(:,3),:);
    % with PCA method, the origin of the point is assumed to be the
    % intersection between the probe's centroid to the probe surface in 
    % the direction of the main principal component.
    [lineplane_intersect_logic, ~, ~, ~, xcoor] = TriangleRayIntersection(origins, direction, vert1, vert2, vert3);
    % get the intersection point from the logical variable (that is
    % returned by the function), then compute the distance
    lineplane_intersectpoints = xcoor(lineplane_intersect_logic, :);
    
    % construct the rigid body for the probes
    probe_current_centroid    = lineplane_intersectpoints;
    % Here, the main principal is the direction of the transducer,
    % since i want that direction to be the z-component of the basis
    % vector, so put the main principal to be the last column of the
    % matrix. Notice also i put negative sign, because before, i force
    % direction of the main principal is where the direction of
    % ultrasound beam. Now, we agree that, by convention, positive
    % z-direction (for the marker and the probe) is away from the bone.
    probe_current_basisvector = [probe_current_PCAbasisvector(:,2), probe_current_PCAbasisvector(:,3), -direction];
    probe_current_rigidbody   = [probe_current_basisvector, probe_current_centroid'; 0 0 0 1];

elseif strcmp(method, 'reg')

    % implementation for registration method here...

end

end

