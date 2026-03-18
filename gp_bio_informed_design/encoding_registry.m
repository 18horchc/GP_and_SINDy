function [configs, config] = encoding_registry(request, flags)
%ENCODING_REGISTRY Mix-and-match encoding configs for bio-informed GPs (acute transient).
%
%   configs = encoding_registry()  returns all predefined configs.
%   config  = encoding_registry('vanilla')  returns the vanilla config by ID.
%   config  = encoding_registry('custom', flags)  builds a config from boolean flags.
%
%   Flags (for 'custom' request): struct with optional fields (default false):
%     .informative_mean   - Fix mean at baseline b=0.1
%     .kernel_matern32    - Use Matérn 3/2 instead of SE
%     .late_pseudo_obs    - Add pseudo-observations at late times near baseline
%     .derivative_constraints - Add derivative sign constraints (rise-then-fall)
%
%   Predefined IDs: 'vanilla', 'mean_only', 'matern_only', 'pseudo_only',
%   'deriv_only', 'mean_matern', 'mean_pseudo', 'mean_deriv', 'matern_pseudo',
%   'matern_deriv', 'pseudo_deriv', 'mean_matern_pseudo', 'mean_matern_deriv',
%   'mean_pseudo_deriv', 'matern_pseudo_deriv', 'full'.
%
%   See also: run_bio_experiment

    if nargin < 1, request = 'all'; end
    if nargin < 2, flags = struct(); end

    % Build config from flags (mix and match)
    if strcmpi(request, 'custom')
        config = build_config_from_flags(flags);
        configs = config;
        return;
    end

    % Predefined configs (vanilla + 15 encoding combinations)
    all_configs = {
        encode_cfg('vanilla',          false, false, false, false),
        encode_cfg('mean_only',        true,  false, false, false),
        encode_cfg('matern_only',      false, true,  false, false),
        encode_cfg('pseudo_only',      false, false, true,  false),
        encode_cfg('deriv_only',       false, false, false, true),
        encode_cfg('mean_matern',      true,  true,  false, false),
        encode_cfg('mean_pseudo',      true,  false, true,  false),
        encode_cfg('mean_deriv',       true,  false, false, true),
        encode_cfg('matern_pseudo',     false, true,  true,  false),
        encode_cfg('matern_deriv',     false, true,  false, true),
        encode_cfg('pseudo_deriv',     false, false, true,  true),
        encode_cfg('mean_matern_pseudo', true,  true,  true,  false),
        encode_cfg('mean_matern_deriv',  true,  true,  false, true),
        encode_cfg('mean_pseudo_deriv',  true,  false, true,  true),
        encode_cfg('matern_pseudo_deriv', false, true,  true,  true),
        encode_cfg('full',             true,  true,  true,  true),
    };

    if strcmpi(request, 'all')
        configs = all_configs;
        config = [];
        return;
    end

    % Single config by ID
    id = lower(request);
    for i = 1:length(all_configs)
        if strcmpi(all_configs{i}.id, id)
            config = all_configs{i};
            configs = config;
            return;
        end
    end

    error('encoding_registry:unknown', 'Unknown config ID: %s', request);
end

function c = encode_cfg(id, use_mean, use_matern, use_pseudo, use_deriv)
    c = struct();
    c.id = id;
    c.label = build_label(use_mean, use_matern, use_pseudo, use_deriv);
    c.informative_mean = use_mean;
    c.kernel_matern32 = use_matern;
    c.late_pseudo_obs = use_pseudo;
    c.derivative_constraints = use_deriv;
end

function lbl = build_label(use_mean, use_matern, use_pseudo, use_deriv)
    parts = {};
    if use_mean,   parts{end+1} = 'mean'; end
    if use_matern, parts{end+1} = 'Matern32'; end
    if use_pseudo, parts{end+1} = 'pseudo'; end
    if use_deriv,  parts{end+1} = 'deriv'; end
    if isempty(parts)
        lbl = 'vanilla';
    else
        lbl = strjoin(parts, '+');
    end
end

function c = build_config_from_flags(flags)
    use_mean   = get_flag(flags, 'informative_mean', false);
    use_matern = get_flag(flags, 'kernel_matern32', false);
    use_pseudo = get_flag(flags, 'late_pseudo_obs', false);
    use_deriv  = get_flag(flags, 'derivative_constraints', false);
    lbl = build_label(use_mean, use_matern, use_pseudo, use_deriv);
    c = encode_cfg(['custom_' lbl], use_mean, use_matern, use_pseudo, use_deriv);
end

function v = get_flag(flags, name, default)
    if isempty(flags) || ~isfield(flags, name)
        v = default;
    else
        v = logical(flags.(name));
    end
end
