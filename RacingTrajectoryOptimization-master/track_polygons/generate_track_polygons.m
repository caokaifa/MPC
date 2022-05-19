function track = generate_track_polygons( checkpoints )

    track = struct;
    
    % Interlace left and right track boundary points
    track.vertices = nan(2, 2*length(checkpoints));
    track.vertices(:,1:2:end) = [checkpoints.left];
    track.vertices(:,2:2:end) = [checkpoints.right];

    % Generate initial set of polygons (quadrilaterals) that cover the track
    for i = 1:length(checkpoints)
        track.polygons(i).vertex_indices = mod1(2*i + (-1:2), size(track.vertices,2));
    end

    track = merge_polygons(track);%1
  
    left_points = [checkpoints(:).left];
    right_points = [checkpoints(:).right];
    forward_vectors = [checkpoints.forward_vector];
    figure
     % Draw track area
    pts=[right_points right_points(:,1) left_points left_points(:,1)];
    fill(pts(1,:), pts(2,:),[1 1 1]*.8,'EdgeAlpha',0)
    % Draw track outline with extra width for the vehicle
    width = .6;    
    normals = width*[0 -1;1 0]*forward_vectors;
    left_points = left_points + normals;
    right_points = right_points - normals;
    plot([left_points(1,:) left_points(1,1)],[left_points(2,:) left_points(2,1)],'k','LineWidth',1);
    hold on
    plot([right_points(1,:) right_points(1,1)],[right_points(2,:) right_points(2,1)],'k','LineWidth',1)
    
    hold on
    for i=1:length(track.polygons)
        dd=[]
           for j=1:length(track.polygons(i).vertex_indices)
                dd=[dd track.vertices(:,track.polygons(i).vertex_indices(j))];
             end
        plot(dd(1,:),dd(2,:),'k','LineWidth',1.5);
    end
    track = add_overlaps(track);%2
     hold on
    for i=1:length(track.polygons)
        dd=[]
           for j=1:length(track.polygons(i).vertex_indices)
                dd=[dd track.vertices(:,track.polygons(i).vertex_indices(j))];
             end
        plot(dd(1,:),dd(2,:),'r','LineWidth',1.5);
    end
    track = add_forward_directions(track, checkpoints);%3
      hold on
    for i=1:length(track.polygons)
        dd=[]
           for j=1:length(track.polygons(i).vertex_indices)
                dd=[dd track.vertices(:,track.polygons(i).vertex_indices(j))];
             end
        plot(dd(1,:),dd(2,:),'g','LineWidth',1.5);
    end
end