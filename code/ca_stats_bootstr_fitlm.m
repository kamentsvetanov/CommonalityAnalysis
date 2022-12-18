function [R2] = ca_stats_bootstr_fitlm(DV,PredV,doRobust)
% Function to bootstrap fitlm and returt desired parameters
% 
% Inputs (defined in cfg.):
%   Model    - model defined as per Wilkinson annotation
%   data     - data in table format with variable names defined in Model
%   doRobust - whether or not to run robust regression
%   varname  - cell array of variable names to output
    
% try doRobust    = cfg.doRobust; catch doRobust = 0; end 
% try Model       = cfg.Model;    catch error('Specify model.'); end
% try data        = cfg.data;     catch error('Specify data table.'); end
% try varname     = cfg.varname;  catch error('Specify at least one output variable.'); end
    
mlr = fitlm(PredV,DV,'RobustOpts',doRobust);

R2 = mlr.Rsquared.Ordinary;

end