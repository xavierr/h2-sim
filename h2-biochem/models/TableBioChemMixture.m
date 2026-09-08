classdef TableBioChemMixture
    % Biochemical reaction parameters for a mixture of microorganisms.
    %
    % SYNOPSIS:
    %   biochem = TableBioChemMixture({'MethanogenicArchae'});
    %
    % PARAMETERS:
    %   names - Cell array of valid reaction names. Use
    %           TableBioChemMixture.getFluidList() to see allowed names.
    %
    % RETURNS:
    %   biochem - Instance with populated property vectors (one entry per
    %             reaction). All vectors are row vectors of length N.
    %
    % SEE ALSO:
    %   bioChemFluidsStructs

    properties
        metabolicReaction   % Names of the metabolic reactions (1 x N string)
        bactnames           % name of bacteria
        rH2                 % Reactant name - H2 (1 x N string)
        rsub                % Reactant name - substrate (CO2, SO4, etc.)
        pH2O                % Product name - water
        p2                  % Product name - second product (e.g., CH4, H2S)
        gamrH2              % Stoichiometric coefficient for H2
        gamrsub             % Stoichiometric coefficient for substrate
        gampH2O             % Stoichiometric coefficient for water
        gamp2               % Stoichiometric coefficient for second product
        Y_H2                % Yield scale used by reaction source terms
        alphaH2             % Half-saturation constant for H2 (mol/mol)
        alphasub            % Half‑saturation constant for substrate (mol/mol)
        Psigrowthmax        % Maximum specific growth rate (1/s)
        bbact               % Decay rate constant (1/s)
        nbactMax            % Population scale used by reaction source terms
        xch_seuil           % Chemotaxis coefficient
        bactdiff            % microbial diffusion coefficient
    end

    properties (Dependent)
        nbioreact           % Number of metabolic reactions
    end

    methods
        function biochem = TableBioChemMixture(names,bactnames)
            % Constructor for TableBioChemMixture.
            if nargin == 0
                names = {};
            end
            if ischar(names) || isstring(names)
                names = {char(names)};
            end
            if nargin == 1
                bactnames = names;
            else
                if ischar(bactnames) || isstring(bactnames)
                    bactnames = {char(bactnames)};
                end
                assert(numel(names) == numel(bactnames));
            end
            % Load the database of known reactions
            allReactions = bioChemFluidsStructs();
            validNames   = {allReactions.metabolicReaction};

            % Validate input
            [isValid, idx] = ismember(lower(names), lower(validNames));
            if ~all(isValid)
                invalid = names(~isValid);
                error('TableBioChemMixture:UnknownReaction', ...
                    ['The following reaction names are not recognized: ', ...
                    strjoin(invalid, ', ')]);
            end

            N = numel(names);
            % Initialize property vectors
            biochem.metabolicReaction = strings(1, N);
            biochem.rH2      = strings(1, N);
            biochem.rsub     = strings(1, N);
            biochem.pH2O     = strings(1, N);
            biochem.p2       = strings(1, N);
            biochem.gamrH2   = zeros(1, N);
            biochem.gamrsub  = zeros(1, N);
            biochem.gampH2O  = zeros(1, N);
            biochem.gamp2    = zeros(1, N);
            biochem.Y_H2     = zeros(1, N);
            biochem.alphaH2  = zeros(1, N);
            biochem.alphasub = zeros(1, N);
            biochem.Psigrowthmax = zeros(1, N);
            biochem.bbact    = zeros(1, N);
            biochem.nbactMax = zeros(1, N);
            biochem.xch_seuil = zeros(1, N);
            biochem.bactdiff = zeros(1, N);

            % Fill from database
            for i = 1:numel(bactnames)
                cts = sum(strcmp(bactnames{i}, bactnames));
                if cts > 1
                    warning(['bacteria ', num2str(i), ': ', bactnames{i}, ' occurs multiple times.']);
                end
            end
            for i = 1:N
                entry = allReactions(idx(i));
                biochem.metabolicReaction(i) = entry.metabolicReaction;
                biochem.rH2(i)      = entry.rH2;
                biochem.rsub(i)     = entry.rsub;
                biochem.pH2O(i)     = entry.pH2O;
                biochem.p2(i)       = entry.p2;
                biochem.gamrH2(i)   = entry.gamrH2;
                biochem.gamrsub(i)  = entry.gamrsub;
                biochem.gampH2O(i)  = entry.gampH2O;
                biochem.gamp2(i)    = entry.gamp2;
                biochem.Y_H2(i)     = entry.Y_H2;
                biochem.alphaH2(i)  = entry.alphaH2;
                biochem.alphasub(i) = entry.alphasub;
                biochem.Psigrowthmax(i) = entry.Psigrowthmax;
                biochem.bbact(i)    = entry.bbact;
                biochem.nbactMax(i) = entry.nbactMax;
                biochem.xch_seuil(i) = entry.xch_seuil;
                biochem.bactdiff(i) = entry.bactdiff;
            end
            biochem.bactnames=bactnames;
        end

        function n = get.nbioreact(biochem)
            n = numel(biochem.metabolicReaction);
        end
    end

    methods (Static)
        function varargout = getFluidList()
            % Return list of available metabolic reaction names.
            allReactions = bioChemFluidsStructs();
            names = {allReactions.metabolicReaction};
            if nargout == 0
                disp('Available biochemical reactions:');
                fprintf('  %s\n', names{:});
            else
                varargout{1} = names(:)';   % row cell array
            end
        end
    end
end