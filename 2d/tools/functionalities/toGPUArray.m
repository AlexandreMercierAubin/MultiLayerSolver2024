function varargout = toGPUArray(varargin)
    varargout = cell(size(varargin));
    for a=1: numel(varargin)
        varargout{a} = gpuArray(varargin{a});
    end
end

