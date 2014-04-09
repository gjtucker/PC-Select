classdef (Abstract) Model < handle
    properties
        name
        use_pc
        supervised
    end

    methods (Abstract)
        run(X, X_test, y, worker_id)
    end

    methods (Static)
        function covar = form_covariates(X, y, use_pc, supervised)
            if use_pc   
                if supervised
                    covar = Model.supervised_PCA(X, y);
                else
                    W = normalize_genotypes(X);
                    %[~, scores] = princomp(W, 'econ');
                    [U, ~] = svd(W*W');
                    covar = U(:, 1:5);
                end
            else
                covar = [];
            end 
            
            % Always add the all ones vector
            covar = [ones(size(y)), covar];
        end
    
        function covar = supervised_PCA(X, y)
            W = normalize_genotypes(X);
            W_resid = W - bsxfun(@times, y'*W/norm(y), y);
            
            W_weighted = bsxfun(@times, W_resid, sqrt(y'*W));

            [~, scores] = princomp(W_weighted, 'econ');
            covar = scores(:, 1:5);
        end
    end
end



