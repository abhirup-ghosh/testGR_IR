function options = EOBSet_mod(varargin)

%EOBSet Create/alter EOB OPTIONS structure.
%
%   OPTION = EOBSet('NAME1',VALUE1,'NAME2',VALUE2,...) creates a structure
%   OPTIONS in which the named properties have the specified values. Any
%   unspecified properties have default values. 
%
%   OPTIONS = EOBSet(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS. Any new properties
%   overwrite corresponding old properties. 
%
%   EOBSet with no input arguments displays all property names and the
%   possible values. 
%


% Properties names
Names = [
    'Dynamics                  ' % EOB dynamics flag
    'RadReac                   ' % Radiation Reaction flag
    'ddotrMethodDyn            ' % how to compute ddotr during dynamics
    'ddotrMethodWav            ' % how to compute ddotr in the waves
    'HorizonFlux               ' % Use horizon flux ?
    'RadialFlux                ' % Use radial flux ?
    'FIterms                   ' % Include Fujita-Iyer term inm rho_lm ?
    'rholmarg                  ' % Argument of rho_lm in DIN flux 
    'PNorder                   ' % PN order to be used
    'resumD                    ' % Resum metric poential D ?
    'resumf22                  ' % Resum f_{22} amplitude
    'Tidal                     ' % Use tidal corrections ?
    'PNTidal                   ' % PN order of tidal correction
    'NewtonWaves               '
    'DetermineNQCab            '
    'fac                       ' % Overall multiplicative factor for EOBFluxDIN_mod (square root multiplies waveform)
    'a1                        ' % Factor multiplying the modes that first contribute to the flux at 1PN (square root multiplies 2,1, 3,3, and 3,1 modes of waveform)
    'a2                        ' % Factor multiplying the modes that first contribute to the flux at 2PN (square root multiplies 3,2, 4,4, and 4,2 modes of waveform)
    'iQNMfac                   ' % Factor multiplying the imaginary part of the QNM
    'NQCFileNR                 '
    'QNMDataDir                '
    'textend                   '
    'ODEEventFunRmin           '
    'ODESolverAbsTol           ' % options for Matlab ODExx ...
    'ODESolverRelTol           '
    'ODESolverRefine           '    
    ];

[m,n]=size(Names);
names = lower(Names);


% Print out possible values of properties
if (nargin == 0) && (nargout == 0)
    
    % Description
    Description = {
        ' string                     {''eob''}                '
        ' string                     {''din''}                '
        ' string                     {''boot2''}              '
        ' string                     {''fd''}                 '
        ' string                     {''yes''}                '
        ' string                     {''yes''}                '
        ' string                     {''yes''}                '
        ' string                     {''v_phi'' }             '
        ' string                     {''5pnlog''}             '
        ' string                     {''pade03''}             '
        ' string                     {''no''}                 '
        ' string                     {''no''}                 '
        ' string                     {''nnlo''}               '
        ' string                     {''no''}                 '
        ' string                     {''no''}                 '
        ' positive scalar            {1.0}                    '
        ' positive scalar            {0.0}                    '
        ' positive scalar            {1.5}                    '
        ' positive scalar            {0.0}                    '
        ' positive scalar            {1.0}                    '
        ' string                     {''''}                   '
        ' string                     {''''}                   '
        ' positive scalar            {1.0}                    '
        ' positive scalar            {100}                    '
        ' positive scalar or vector  {1e-12}                  '        
        ' positive scalar            {1e-9}                   '
        ' positive integer           {4}                      '
        };

    Names = mat2cell(Names,ones(m,1),n);
    cellfun(@(x,y) fprintf(' %s [ %s ] \n',x,y),Names,Description);
    
    return;
    
end


% Combine all leading options structures o1, o2, ... in odeset(o1,o2,...)
options = [];
for j = 1:m
  options.(deblank(Names(j,:))) = [];
end
i = 1;
while i <= nargin
  arg = varargin{i};
  if ischar(arg)                        % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error('EOBSet:NoPropNameOrStruct',...
            ['Expected argument %d to be a string property name ' ...
                     'or an options structure\ncreated with EOBSet.'], i);
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
        val = arg.(deblank(Names(j,:)));
      else
        val = [];
      end
      if ~isempty(val)
        options.(deblank(Names(j,:))) = val;
      end
    end
  end
  i = i + 1;
end


% A finite state machine to parse name-value pairs
if rem(nargin-i+1,2) ~= 0
  error('EOBSet:ArgNameValueMismatch',...
        'Arguments must occur in name-value pairs.');
end
expectval = 0;                       % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};
    
  if ~expectval
    if ~ischar(arg)
      error('EOBSet:NoPropName',...
            'Expected argument %d to be a string property name.', i);
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error('EOBSet:InvalidPropName',...
            'Unrecognized property name ''%s''.', arg);
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
        msg = sprintf('Ambiguous property name ''%s'' ', arg);
        msg = [msg '(' deblank(Names(j(1),:))];
        for k = j(2:length(j))'
          msg = [msg ', ' deblank(Names(k,:))];
        end
        msg = sprintf('%s).', msg);
        error('EOBSet:AmbiguousPropName', msg);
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    options.(deblank(Names(j,:))) = arg;
    expectval = 0;
      
  end
  i = i + 1;
end

if expectval
  error('EOBSet:NoValueForProp',...
        'Expected value for property ''%s''.', arg);
end
