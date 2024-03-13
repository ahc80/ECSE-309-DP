%% Function computeWire
 % -----------------------------------------------------------------------
 % This function draws a piecewise linear wire in a 3D space by returning
 % an array of points on a linearly spaced grid
 %
 % Inputs: 
 %       1) V : coordinates of joints in the wire (x y z)        [Nx3]
 %       2) spacing : the spacing over which to return coordinates of the
 %       wires path in space
 %

function [wirePoints] = computeWire(V,spacing)

    wirePoints(1,:) = V(1,:);

    for k=2:length(V)
        xvec = (V(k,1) - V(k-1,1));
        yvec = (V(k,2) - V(k-1,2));
        zvec = (V(k,3) - V(k-1,3));
        mag = sqrt (xvec^2 + yvec^2 + zvec^2);

        a = xvec/mag; %x-slope
        b = yvec/mag; %y-slope
        c = zvec/mag; %z-slope
        
        points = [spacing:spacing:mag];

        xlin = points.*a + V(k-1,1);
        ylin = points.*b + V(k-1,2);
        zlin = points.*c + V(k-1,3);
        
        newpts = [xlin; ylin; zlin]';
        
        wirePoints = vertcat(wirePoints, newpts);

    end
end