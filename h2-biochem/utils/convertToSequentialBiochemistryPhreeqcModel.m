function seqModel = convertToSequentialBiochemistryPhreeqcModel(model, varargin)
%Convert an existing BiochemistryPhreeqcModel to a SequentialBiochemistryPhreeqcModel.
%
% SYNOPSIS:
%   seqModel = convertToSequentialBiochemistryPhreeqcModel(model, ...
%       'reactionSubstepMaxDt', 0.4*day)
%
% DESCRIPTION:
%   Converts an already fully configured BiochemistryPhreeqcModel -- for
%   example the model returned by setupH2StorageExampleWithSRB_benchmark
%   with phreeqcBackend='sequential-h2biochem-phreeqc' and
%   phreeqcTimestepCoupling=true -- into an equivalent
%   SequentialBiochemistryPhreeqcModel that runs the general
%   coarse-flow/local-reaction split (one global flow+transport solve per
%   control step, followed by local biochemical-reaction/PHREEQC
%   substeps) instead of simulateSequentialH2BiochemPhreeqc's outer
%   Picard iteration.
%
%   Every public property of MODEL -- including customization applied
%   after its own constructor ran, such as a replaced EOSModel or
%   bact_capProp/bact_maxProp -- is copied onto the new model; see
%   copyBiochemistryModelProperties. Model-defining options that
%   SequentialBiochemistryPhreeqcModel's own constructor validates
%   immediately (phreeqcBackend, phreeqcTimestepCoupling, carbonateBuffer,
%   ...) are also forwarded directly to that constructor so the model can
%   be built in one step; this leaves nothing for the later property copy
%   to disagree with.
%
% REQUIRED PARAMETERS:
%   model - A BiochemistryPhreeqcModel (or subclass, but not already a
%           SequentialBiochemistryPhreeqcModel) instance with
%           phreeqcBackend='sequential-h2biochem-phreeqc' and
%           phreeqcTimestepCoupling=true.
%
% OPTIONAL PARAMETERS:
%   reactionSubstepMaxDt - Maximum duration of one local reaction/PHREEQC
%                          substep. Required (no default): must be
%                          supplied as a property/value pair.
%   Any other property/value pair accepted by
%   SequentialBiochemistryPhreeqcModel (e.g. flowNonLinearSolver,
%   reactionNonLinearSolver), applied after the property copy so they take
%   precedence over MODEL's own configuration.
%
% RETURNS:
%   seqModel - SequentialBiochemistryPhreeqcModel instance, already
%              validated (its flow-stage/reaction-stage sub-models are
%              built and ready to simulate).
%
% SEE ALSO:
%   SequentialBiochemistryPhreeqcModel, BiochemistryPhreeqcModel,
%   copyBiochemistryModelProperties, setupH2StorageExampleWithSRB_benchmark

    assert(isa(model, 'BiochemistryPhreeqcModel'), ...
        ['convertToSequentialBiochemistryPhreeqcModel requires a ', ...
         'BiochemistryPhreeqcModel (or subclass) instance.']);
    assert(~isa(model, 'SequentialBiochemistryPhreeqcModel'), ...
        'model is already a SequentialBiochemistryPhreeqcModel.');
    assert(model.phreeqcTimestepCoupling && model.isSequentialH2BiochemPhreeqcBackend(), ...
        ['model must use phreeqcBackend=''sequential-h2biochem-phreeqc'' ', ...
         'with phreeqcTimestepCoupling=true for the coarse-flow/local-', ...
         'reaction split to be meaningful.']);

    maxDt = getRequiredOption(varargin, 'reactionSubstepMaxDt', ...
        ['convertToSequentialBiochemistryPhreeqcModel requires ', ...
         '''reactionSubstepMaxDt'' to be supplied as a property/value pair.']);

    % Model-defining options that the constructor validates immediately
    % (before any later property copy could supply them), forwarded
    % directly from the template model.
    requiredOptions = { ...
        'bacteriamodel', model.bacteriamodel, ...
        'carbonateBuffer', model.carbonateBuffer, ...
        'carbonateBufferPH', model.carbonateBufferPH, ...
        'carbonateBufferPka1', model.carbonateBufferPka1, ...
        'phreeqcTimestepCoupling', model.phreeqcTimestepCoupling, ...
        'phreeqcBackend', model.phreeqcBackend, ...
        'phreeqcDatabaseFile', model.phreeqcDatabaseFile, ...
        'phreeqcComProgId', model.phreeqcComProgId, ...
        'phreeqcCouplingOptions', model.phreeqcCouplingOptions, ...
        'bacterialDecayOrder', model.bacterialDecayOrder, ...
        'reactionSubstepMaxDt', maxDt};

    seqModel = SequentialBiochemistryPhreeqcModel(model.G, model.rock, model.fluid, ...
        model.compFluid, model.biochemFluid, model.water, model.AutoDiffBackend, ...
        requiredOptions{:});

    % Copy every remaining shared public property (EOSModel,
    % bact_capProp/bact_maxProp, molecularDiffusion, gammak, ...).
    seqModel = copyBiochemistryModelProperties(model, seqModel);

    % Apply any final caller overrides (including reactionSubstepMaxDt
    % again, harmlessly, plus e.g. flowNonLinearSolver/reactionNonLinearSolver).
    seqModel = merge_options(seqModel, varargin{:});

    % Rebuild the flow-stage/reaction-stage sub-models from the fully
    % copied and overridden configuration, and let validateModel run any
    % other standard model-consistency checks.
    seqModel = seqModel.validateModel();
end

function value = getRequiredOption(options, name, errmsg)
value = [];
found = false;
for i = 1:2:numel(options)
    if strcmpi(options{i}, name)
        value = options{i + 1};
        found = true;
    end
end
assert(found && ~isempty(value), errmsg);
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
