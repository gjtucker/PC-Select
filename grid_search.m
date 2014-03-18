function [val, v] = grid_search(func, grid)
v = zeros(length(grid), 1);
for i = 1:length(grid)
 %   try
        v(i) = func(grid(i));
  %  catch err
   %     err
    %    v(i) = NaN;
   % end
end
[~, i] = min(v);
val = grid(i);
end