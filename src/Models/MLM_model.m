classdef MLM_model < Model  
    methods
        function m = MLM_model(use_pc, supervised)
            m.use_pc = use_pc;
            if nargin < 2
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

            [ ~, wald_stat, p_val ] = fast_lmm( X, X_test, y, covar, ...
                    worker_id );
            
            % Extra results that are specific to this type of model
            extra = struct();
        end
    end
end
