function varargout = toCPUArray(varargin)
    varargout = cell(size(varargin));
    for a=1: numel(varargin)
        varargout{a} = gather(varargin{a});
    end
end

