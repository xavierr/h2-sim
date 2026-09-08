function target = copyBiochemistryModelProperties(source, target)
%Copy every shared public property from SOURCE onto TARGET.
%
% SYNOPSIS:
%   target = copyBiochemistryModelProperties(source, target)
%
% DESCRIPTION:
%   Copies every publicly gettable/settable property of TARGET that also
%   exists (and is publicly readable) on SOURCE, using metaclass
%   reflection. This is used by SequentialBiochemistryPhreeqcModel to
%   build its internal flow-stage/reaction-stage sub-models so they share
%   an identical configuration with the fully configured meta-model they
%   were built from (including any options set after the base
%   BiochemistryPhreeqcModel constructor ran, e.g. a custom EOSModel or
%   bact_capProp/bact_maxProp), without maintaining a manually duplicated
%   list of constructor option names. It is also reused by
%   convertToSequentialBiochemistryPhreeqcModel to convert an existing,
%   fully configured BiochemistryPhreeqcModel benchmark model into a
%   SequentialBiochemistryPhreeqcModel.
%
%   Properties that are Constant, Transient, Dependent, or not publicly
%   accessible on either SOURCE or TARGET are skipped, as are properties
%   that only one of the two objects has (e.g. SequentialBiochemistry-
%   PhreeqcModel-only properties are simply not present on the plain
%   BiochemistryPhreeqcModel stage sub-models, so they are silently
%   skipped when copying onto those sub-models).
%
% PARAMETERS:
%   source - Fully configured model instance to copy FROM.
%   target - Model instance to copy properties ONTO.
%
% RETURNS:
%   target - TARGET with every shared public property set equal to the
%            corresponding property on SOURCE.
%
% SEE ALSO:
%   SequentialBiochemistryPhreeqcModel,
%   convertToSequentialBiochemistryPhreeqcModel

    assert(isobject(source) && isobject(target), ...
        'copyBiochemistryModelProperties requires two class instances.');

    % NOTE: for plain (non-handle) classdef objects, metaclass() returns a
    % matlab.metadata.Class instance whose findprop(name) method does not
    % reliably resolve statically defined properties by name (it can
    % return empty even for a property that is clearly present in
    % PropertyList), and findprop(instance, name) does not apply at all
    % to non-handle objects. PropertyList is therefore searched directly
    % (by name) on both source and target instead of using findprop.
    mcTarget = metaclass(target);
    mcSource = metaclass(source);
    for i = 1:numel(mcTarget.PropertyList)
        prop = mcTarget.PropertyList(i);
        name = prop.Name;
        if prop.Constant || prop.Transient || prop.Dependent
            continue
        end
        if ~isPublicAccess(prop.SetAccess) || ~isPublicAccess(prop.GetAccess)
            continue
        end
        if ~isprop(source, name)
            continue
        end
        srcProp = findPropertyByName(mcSource, name);
        if isempty(srcProp) || srcProp.Constant || srcProp.Transient || ...
                srcProp.Dependent || ~isPublicAccess(srcProp.GetAccess)
            continue
        end
        target.(name) = source.(name);
    end
end

function prop = findPropertyByName(mc, name)
% Search a matlab.metadata.Class/meta.class PropertyList by name. See the
% NOTE above copyBiochemistryModelProperties for why findprop is avoided.
    prop = [];
    for i = 1:numel(mc.PropertyList)
        if strcmp(mc.PropertyList(i).Name, name)
            prop = mc.PropertyList(i);
            return
        end
    end
end

function tf = isPublicAccess(access)
% Property access specifiers ('SetAccess'/'GetAccess') are usually the
% character vector 'public', but MATLAB returns a cell array of access
% specifiers (one per class in the hierarchy) for some properties; be
% conservative and require every listed specifier to be 'public'.
    if iscell(access)
        tf = ~isempty(access) && all(cellfun(@(a) ischar(a) && strcmp(a, 'public'), access));
    else
        tf = ischar(access) && strcmp(access, 'public');
    end
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
