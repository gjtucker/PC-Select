classdef Lin_model < Model    
    methods
        function m = Lin_model(use_pc, supervised)
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
        
            [ ~, wald_stat, p_val ] = lin_reg(X_test, y, covar );
            
            % Extra results that are specific to this type of model
            extra = struct();
        end
    end
end