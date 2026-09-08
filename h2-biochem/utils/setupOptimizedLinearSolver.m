function [nls, lsolve] = setupOptimizedLinearSolver(model, varargin)
% Setup optimized AMGCL CPR linear solver for compositional biochemical models
%
% SYNOPSIS:
%   [nls, lsolve] = setupOptimizedLinearSolver(model)
%   [nls, lsolve] = setupOptimizedLinearSolver(model, 'complexityLevel', 'high')
%
% PARAMETERS:
%   model           - BiochemistryModel instance
%   complexityLevel - 'low' (advection only), 'medium' (with diffusion),
%                     'high' (with both diffusion & dispersion) [default: 'medium']
%   solverTolerance - Tolerance for linear solver [default: 1e-4]
%   maxNonlinIter   - Max nonlinear iterations [default: 15]
%   cprDamp         - Retained for backward-compatible calls; AMGCL CPR
%                     does not expose this setting.
%
% RETURNS:
%   nls    - NonLinearSolver configured with linear solver
%   lsolve - CPRSolverAD linear solver with optimized AMGCL parameters

    opt = struct('complexityLevel', 'medium', ...
                 'solverTolerance', 1e-4, ...
                 'maxNonlinIter', 15, ...
                 'cprDamp', []);
    opt = merge_options(opt, varargin{:});

    % Select base linear solver
    lsolve = selectLinearSolverAD(model);

    % Configure AMGCL CPR parameters based on problem complexity
    amgclSettings = struct('max_levels', 20, 'verbose', false);

    switch lower(opt.complexityLevel)
        case 'low'
            % Simple advection-only: aggressive coarsening
            amgclSettings.aggr_eps_strong = 1e-1;
            amgclSettings.maxIterations = 50;

        case 'medium'
            % Diffusion included: balanced settings (default)
            amgclSettings.aggr_eps_strong = 1e-2;
            amgclSettings.maxIterations = 100;

        case 'high'
            % Both diffusion & dispersion: conservative coarsening
            amgclSettings.aggr_eps_strong = 5e-3;
            amgclSettings.maxIterations = 150;

        otherwise
            error('Unknown complexityLevel: %s', opt.complexityLevel);
    end

    % Apply AMGCL settings through the solver's actual configuration
    % interface. Small systems may use the backslash fallback instead.
    if isa(lsolve, 'AMGCLSolverAD')
        lsolve.setCoarsening('aggregation');
        lsolve.amgcl_setup.max_levels = amgclSettings.max_levels;
        lsolve.amgcl_setup.verbose = amgclSettings.verbose;
        lsolve.amgcl_setup.aggr_eps_strong = amgclSettings.aggr_eps_strong;
        lsolve.maxIterations = amgclSettings.maxIterations;
    end
    lsolve.tolerance = opt.solverTolerance;

    % Configure nonlinear solver
    nls = NonLinearSolver();
    nls.LinearSolver = lsolve;
    nls.maxIterations = opt.maxNonlinIter;
    nls.verbose = false;
end

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
