function bc = addBacterialBC(bc, faces, nbactval)
% Attach Dirichlet bacterial concentrations to a boundary-condition struct.
%
% SYNOPSIS:
%   bc = addBacterialBC(bc, faces, nbactval)
%
% DESCRIPTION:
%   Sets the prescribed bacterial concentration `nbact` on a set of
%   boundary faces of an existing MRST boundary-condition struct, for use
%   with the bacterial diffusion Dirichlet conditions implemented in
%   `BiochemistryModel.addBacterialDiffusionBC`.
%
%   The values are stored in the extra field `bc.nbact`, aligned with
%   `bc.face`. Faces that are not given a value keep NaN, i.e. the natural
%   no-flux condition for the bacterial equation.
%
%   The target faces must already be present in `bc.face`. If a boundary is
%   meant to be impermeable to flow but still impose a bacterial
%   concentration, add a zero-flux condition there first, e.g.
%       bc = addBC([], faces, 'flux', 0, 'sat', sat);
%       bc = addBacterialBC(bc, faces, nbactVal);
%
% PARAMETERS:
%   bc       - Boundary-condition struct from `addBC` (must be non-empty).
%   faces    - Boundary face indices to assign a bacterial value to.
%   nbactval - Prescribed bacterial concentration. Either a scalar (applied
%              to all `faces`) or a vector matching numel(faces).
%
% RETURNS:
%   bc - Updated struct with field `bc.nbact` (numel(bc.face) x 1).
%
% SEE ALSO:
%   addBC, BiochemistryModel

    assert(~isempty(bc) && isfield(bc, 'face'), ...
        'A non-empty boundary-condition struct (from addBC) is required.');

    faces = faces(:);
    if isscalar(nbactval)
        nbactval = repmat(nbactval, numel(faces), 1);
    else
        nbactval = nbactval(:);
        assert(numel(nbactval) == numel(faces), ...
            'nbactval must be scalar or match numel(faces).');
    end

    % Initialise the storage on first use (NaN = no bacterial Dirichlet BC)
    if ~isfield(bc, 'nbact') || isempty(bc.nbact)
        bc.nbact = nan(numel(bc.face), 1);
    end

    [tf, loc] = ismember(faces, bc.face);
    assert(all(tf), ...
        ['All bacterial BC faces must already exist in bc.face. ', ...
         'Add a flux/pressure condition on these faces first.']);

    bc.nbact(loc) = nbactval;
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
