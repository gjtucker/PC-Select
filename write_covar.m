function write_covar(filename, covar )
% Write covariates
f = fopen(filename, 'W');

n = size(covar, 1);
for i = 1:n
    fprintf(f, '1\tIND%d', i);
    
    for j = 1:size(covar, 2)
        fprintf(f, '\t%f', covar(i, j));
    end
    fprintf(f, '\n');
end
fclose(f);

end
