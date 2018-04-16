% pedmdraw(x,P,rgbvec)
%   draw a realization in 2 or 3D with the graph edges

function pedmdraw(x,P,rgbvec)
  if (nargin < 3)
    rgbvec = [0,0,1];
  end
  [K,n] = size(x);
  labels = cellstr(num2str([1:n]'));
  clf
  hold on
  if (K == 2) 
    plot(x(1,:), x(2,:), 'o', 'color', rgbvec);
    text(x(1,:), x(2,:), labels, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    for i=1:n
      for j=1:n
        if (P(i,j) > 0) 
          line([x(1,i),x(1,j)],[x(2,i),x(2,j)], 'color', rgbvec); 
        end
      end
    end
  elseif (K == 3)  
    plot3(x(1,:), x(2,:), x(3,:), 'o', 'color', rgbvec);
    text(x(1,:), x(2,:), x(3,:), labels, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
    for i=1:n
      for j=1:n
        if (P(i,j) > 0) 
          line([x(1,i),x(1,j)],[x(2,i),x(2,j)],[x(3,i),x(3,j)], 'color', rgbvec); 
          %Cylinder(x(:,i), x(:,j), 0.1, 10, 'b', 1, 0);
        end
      end
    end
    view(3);
  else 
    printf('showrealization requires 2- or 3-dimensional realizations only\n');
  end
  grid on;
end