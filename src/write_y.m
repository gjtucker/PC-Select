function write_y(filename, y)
% Write phenotype matrix
n = length(y);

f = fopen(filename, 'W');

fprintf(f, 'FID\tIID\tPheno\n');
for i = 1:n
    fprintf(f, '1\tIND%d\t%f\n', i, y(i));
end

fclose(f);

end
