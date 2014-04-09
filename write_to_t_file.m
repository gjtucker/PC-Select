function write_to_t_file(base_filename, X, X_test, y, covar, write_sim )

if nargin < 6
    write_sim = false;
end

if write_sim
    write_sim_matrix([base_filename, '_sim'], X);
else
    write_genotype_matrix([base_filename, '_cov'], X);
end

write_genotype_matrix([base_filename, '_test'], X_test);

write_y([base_filename, '_pheno.txt'], y);
write_covar([base_filename, '_covar.txt'], covar);

end

function write_sim_matrix(base_filename, X)
n = size(X, 1);

f_sim = fopen([base_filename, '.txt'], 'W');
fprintf(f_sim, 'var');
for i = 1:n
    fprintf(f_sim, '\t1 IND%d', i);
end
fprintf(f_sim, '\n');

for i = 1:n
    fprintf(f_sim, '1 IND%d', i);
    fprintf(f_sim, '\t%g', X(i, :));
    fprintf(f_sim, '\n');
end

fclose(f_sim);

end

function write_genotype_matrix(base_filename, X)
[n, m] = size(X);

% Write tfam file
f_fam = fopen([base_filename, '.tfam'], 'W');
for i = 1:n
    fprintf(f_fam, '1\tIND%d\t0\t0\t0\t1\n', i);
end
fclose(f_fam);

% Write tped file
ped_matrix = zeros(m, 2*n);
ped_matrix(:, 1:2:end) = (X' > 0) + 1;
ped_matrix(:, 2:2:end) = (X' > 1) + 1;
ped_matrix = [ones(m, 1), (1:m)', zeros(m, 2), ped_matrix];

% tic;
% f_ped = fopen([base_filename, '.tped'], 'W');
% for i = 1:m
%     fprintf(f_ped, '%d\t', ped_matrix(i, :));
%     fprintf(f_ped, '\n');
% end
% fclose(f_ped);
% toc

fast_write([base_filename, '.tped'], ped_matrix');
%tic; dlmwrite([base_filename, '.tped'], ped_matrix, '\t'); toc

end