function vout = SetVOutput( n, varargin )

%SetVOutput Manage variables in output.
%
%   vout = SetVOutput( n, A,B,C,... ) return a structure with fields
%   A,B,C,... if n=1, a cell array {A,B,C,...} otherwise
%

na = length(varargin);    

if na<2
    error('At least a variable is required.')
end

if n==1
                
    MyS = struct();
    for k=1:na
        MyS.(inputname(k+1)) = varargin{k};
    end
    vout = {MyS};
    
    return;
    
end

vout = varargin;