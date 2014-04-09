classdef Fast_lmm_autoselect_model < Model  
    properties
        cutoffs
    end
    
    methods
        function m = Fast_lmm_autoselect_model(cutoffs, use_pc, supervised)
            m.cutoffs = cutoffs;
            m.use_pc = use_pc;
            if nargin < 3
                m.supervised = false;
            else
                m.supervised = supervised;
            end
        end

        function [ wald_stat, p_val, extra ] = ...
            run(m, X, X_test, y, worker_id)
        
            covar = Model.form_covariates(X, y, m.use_pc, m.supervised);

            % FLS doesn't use the all ones covariate
            covar = covar(:, 2:end);
        
            [ ~, wald_stat, p_val, mse, log_like, n_snps_chosen ] = ...
                fast_lmm_autoselect( X, X_test, y, covar, m.cutoffs, ...
                worker_id );
            
            % Extra results that are specific to this type of model
            extra = struct();
            extra.mse = mse;
            extra.log_like = log_like;
            extra.n_snps_chosen = n_snps_chosen;
        end
    end
end
    
