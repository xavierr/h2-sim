classdef PorosityFromRock < StateFunction
    % Simple porosity provider for cases without bio-clogging.
    % Returns static porosity from model.rock.poro.

    methods
        function sf = PorosityFromRock(model, varargin)
            sf@StateFunction(model, varargin{:});
            sf.label = '\phi';
        end

        function phi = evaluateOnDomain(sf, model, state) %#ok<INUSD>
            phi = model.rock.poro;
        end
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