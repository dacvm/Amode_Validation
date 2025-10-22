function [plane_center_globalframe, plane_basisvector, marker_points_esterror] = estimateMarkerRigidBody(marker_points)
% Estimate the plane (the normal)
% To obtain the rigid body transformation of the markers, we need to
% estimate the plane that is crossing the four points. Of course, there
% will be no exact solution, since there is exist noise, but we can
% estimate the plane. This can be done using the function provided here
% alled affine_fit(), however, even though it does estimate the plane
% (by giving us plane normal, plane basis vector, origin point), but
% they are not as exactly what we want.
%
% The plane normal is fine, we can use that. But, the other two results
% are not what we want. First, the Origin Point is calculated by the
% mean of the points, which is not the case as we don't want the mean
% of the markers as the origin. If you look again closely to our marker
% arrangement as the figure below:
% 1        4
%  o-------o
%  |        \
%  |         \
%  |    v x   \
%  |           \
%  |            \
%  o-------------o
% 2              3
% the mean point is represented by (v) and the actual center we want is
% the (x). Why? because that was the arrangement when we do the
% experiemnt, that is the point where the screw is, the center of the
% "platform". (v) is more towards to the left, as the center of mass of
% these points are in the left area.
%
% Second problem is the Plane Basis Vector. affine_fit() function just
% calculate arbitraty two vectors that is perpendicular to the normal,
% so it does not represent the geometry of our markers. 
% 
% So, this is what we want. We want the origin of the plane in the point
% x in the figure above, and we want (2-3) as the direction of X and 
% (2-1) as the direction of Y.


% 1) Calculate plane vector X, Y, and Z -------------------------------

% fit the plane, and obtain only the plane normal
[plane_normal, ~, ~] = affine_fit(marker_points);
% Calculate vector X and vector Y
plane_vector_x  = (marker_points(3,:) - marker_points(2,:)) / norm((marker_points(3,:) - marker_points(2,:)));
plane_vector_y  = (marker_points(1,:) - marker_points(2,:)) / norm((marker_points(1,:) - marker_points(2,:)));
% Calculate vector Z by cross product between the two vector. 
plane_vector_z  = cross(plane_vector_x, plane_vector_y);

% Okay, you might ask, why the hell you do this, you already get the 
% normal vector from affine_fit() function? The problem is, the direction 
% of the normal from that funtion is also arbitrary, it can be upward or
% downward (relative to the plane). No, we don't want that. We want
% certainty, we will follow left-hand rules. If you consider the
% diagram above, normal vector should be pointing out of the screen.
%
% However, i am believing more towards normal vector that is estimated
% by affine_fit() then by the cross product. Since the function also
% accomodate the fourth point. So, i will take the direction but i will
% change the sign according to the vector from cross product
sign_corrector  = sign(dot(plane_vector_z, plane_normal'));
plane_basisvector = [plane_vector_x', plane_vector_y', sign_corrector*plane_normal];

% 2) Calculate the center of the plane -------------------------------

% It is more difficult to calculate it in 3D, so what i will do, i will
% consider the markers in their local reference frame, in 2D. We can do 
% this by transforming our markers with the invers of its rigid body 
% transformation matrix. But we only have the basis vector (the
% rotation matrix), we still dont have the translation matrix? We can
% use Point 2 as our temporary origin point, aka translation vector.
marker_points_localframe = homogeneous2cartesian( inverseHMat([plane_basisvector, marker_points(2,:)'; 0 0 0 1]) * cartesian2homogeneous(marker_points) );
% We only consider the x and y column. There will be small number of Z
% as the plane is not perfectly fitting through all the points, but it
% is small enough that we can ignore it.
marker_points_2d         = marker_points_localframe(:,1:2);
% calculate the center of the "square" by only using Point 1-2-3
plane_center_2d          = [(marker_points_2d(2,1)+marker_points_2d(3,1))/2, (marker_points_2d(1,2)+marker_points_2d(2,2))/2];
% return the center back to 3D
plane_center_globalframe = homogeneous2cartesian( [plane_basisvector, marker_points(2,:)'; 0 0 0 1] * cartesian2homogeneous([plane_center_2d, 0]) );

% 3) Calculate marker placement estimation errors ---------------------

% This part is an addition. I notice that the shape of our markers are 
% not exactly like right-trapezoids, it skewed. This can be caused by:
% 1) 3D printing error
% 2) Human error (when placing the markers)
% 3) Segmentation error (from the segmentation software)
% For the sake of clarity, i will also calculate the placement error 
% (relative to Point 2). Since i know exactly how it should looks like, 
% so i put the numbers here. As you can see, Point 2 is (0,0) since we 
% are talking about relative placement error.
marker_points_2dGT       = [0 10; 0 0; 10 0; 5 10];
% calculate the mean square error.
marker_points_esterror   = mean ( sqrt( sum( (marker_points_2dGT-marker_points_2d).^2, 2 ) ) );
end

