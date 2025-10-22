function marker_order = rearrangeMarker(marker_points)
% Rearrange our markers so that we can do next operations easily.
% Please take a not. The code i write here is not the only solution of
% this problem, and obviously not the optimal solution, i just somehow 
% got this path.
%
% The fiducials (marker) is arranged like diagram below:
% o----o
% |     \
% |      \
% o-------o
% The goal is to make arrangement: (1) top-left, (2) botom-left, (3)
% bottom-right, (4) top-right. However, this is not an easy task, since
% we only have relational-distance information, nothing more. So we
% need to be smart looking the features.
% 
% In the design, (2-1) and (2-3) supposed to be the same (10mm) 
% however, it is not the case here as you can see from the distance 
% matrix below. It might caused by (a) human error when placing, (b)
% detection error from segmentation. If you feel this is not good, use
% other alternative, probe point cloud registration. If you feel, "yeah
% it is what it is" let's go on.
%
% When i checked the distance matrix below, looking for ~10mm distance
% between two points (to look for (2-1) and (2-3)) is not reliable.
% Then i realized there are 2 important feature.
% 1) distance between (1-4) will be always the smallest.
% 2) distance between (1-3) will be always the longest.
% 3) point (4-2) and (4-3) is always giving similar distance,
%    as (4-2-3) forms isoceles triangle. 
% 
% It is a good start, but using information 1 and 2 we cant determine
% which point is which, because the information we got is only the 
% relational distance. So, i start with information 3. We could 
% determine Point 4 because it is unique, it is connected to two other 
% point that has similar distance, and connected to a point with
% smallest distance.
%
% After we could determine Point 4, we can determine point 1, by
% looking the shortest connection. After we could determine Point 1, we
% could determine Point 3 by looking the longest link from Point 1.
% Then we can determine the last point is Point 2.

% 1) create distance matrix -------------------------------------------

dist_matrix = zeros(size(marker_points,1), size(marker_points,1));
for i=1:size(marker_points,1)
    for j=1:size(marker_points,1)
        dist_matrix(i,j) = sqrt(sum((marker_points(i,:) - marker_points(j,:)).^2));
    end
end

% 2) Now, i am looking for a link with a smallest distance ------------

% Here, because the diagonal of a dist_matrix is zero (a distance of a
% point to itself), i changed it to big number, 999, so that it will 
% become easy to look for the smallest distance
tmp = dist_matrix + eye(size(dist_matrix))*999;
% search for smallest distance
[min_val, min_idx] = min(tmp, [], 'all');
% determine the 2d matrix index of that smallest distance within the
% distance matrix
[potential_row, potential_col] = ind2sub(size(tmp), min_idx);

% Now, the information that we just got is ONLY a relational distance,
% it doesn't gave us which point is point number 4. Look at the diagram
% below, it can be left or it can be right.
% o----o        o----o
% | \               / \
% |   \     or     /   \ 
% |     \         /     \
% o       o     o        o
%
% This relation is represented in the distance matrix. Look at the
% example distance matrix below:
% dist_matrix =
%          0    4.8596   10.9908   14.5730
%     4.8596         0   12.0130   11.8411
%    10.9908   12.0130         0    9.8485
%    14.5730   11.8411    9.8485         0
% We can see that even though we know the smallest distance is at 
% dist_matrix(2,1), but we don't know if Point 4 is represented by 2 or
% by 1.
%
% Still focus on the smallest distance value, relating to the possible 
% diagram above, the other two distance links are represented in row-wise 
% and column-wise. We need to determine which one we should look for. 
% Luckily, we know that Point 4 has two similar distance connection. 
% So, that is what we look for now.

% 3) Determining the correct connection -------------------------------

% let's take two possible connection, row-wise and column-wise
possible1 = dist_matrix(potential_row, :);
possible2 = dist_matrix(:, potential_col);

% Possible1 and possible2 is a vector of distance between point
% (connection/link). To look for the "similar distance", we can take
% the absolute different between the two connection. If we find a value
% that is close to 0, it means we find it.
% 
% We can do thresholding with a certain number (for example 0.1), but i 
% tried to avoid threshold. So what i do, i just find the smallest
% distance between connections in possible1 and possible2. Then between 
% the two, i took the one that produces the smallest.

% Since there are 3 connections, there will be 3*2*1 combination 
% distance we need to check.
possible1_diff = [];
possible2_diff = [];
for i=1:length(possible1)
    for j=i+1:length(possible1)
        diff1 = abs(diff([possible1(i), possible1(j)]));
        diff2 = abs(diff([possible2(i), possible2(j)]));

        possible1_diff= [possible1_diff, diff1];
        possible2_diff= [possible2_diff, diff2];
    end
end

% take the minimum distance value between the connections both for
% possible1 and possible2
[min1_val, min1_idx] = min(possible1_diff, [], 'all');
[min2_val, min2_idx] = min(possible2_diff, [], 'all');

% if min1_val is lesser than min2_val, it means that possible1 has the
% "similar distance" connection.
if(min1_val<min2_val)

    % if possible1 is the answer, it means that we should look
    % row-wise. It means that Point 4 is determined by the row index
    % from the smallest distance we got within the distance matrix (see
    % the distance matrix example above)
    point4 = potential_row;
    point1 = potential_col;
else    

    % if possible2 is the answer, it means that we should look
    % column-wise. It means that Point 4 is determined by the column
    % index from the smallest distance we got within the distance
    % matrix (see the distance matrix example above)
    point4 = potential_col;
    point1 = potential_row;    
end

% 4) Determining Point 2 and 3 ----------------------------------------

% At this point, we know already which one is Point 1 and Point 4. Now,
% to determine the other point, we can start with Point 3, which has
% the longest distance to Point 1. So, let's look at the connections
% from Point 1. We can do this by taking the entire row of Point 1 (or,
% you can also choose the entire column of Point 1)
tmp = dist_matrix(point1,:);
% the index of the maximum value will be Point 3
[max_val, max_idx] = max(tmp, [], 'all');
point3 = max_idx;
% the second maximum value will be Point 4. Here, i just zeroed the
% last maximum value, so that if i, again, use max function, i got the
% next maximum value.
tmp(max_idx) = 0;
[max_val, max_idx] = max(tmp, [], 'all');
point2 = max_idx;

% Voila, we got the marker order
marker_order = [point1, point2, point3, point4];
end

