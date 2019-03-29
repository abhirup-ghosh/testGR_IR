function WriteASCII(fname,x,y,varargin)

%WriteASCII write ASCII data in a GNUPLOT file
%
%   WriteASCII(fname,x,Y);
%
%   WriteASCII(fname,x,Y, header);
%   Add a header specified by the string header
%
%   WriteASCII(fname,x,Y, [],info);
%   Print info (if info)
%


% default opt
header = [];
info   = 0;


% args in
optargs = {header info};
na = length(varargin);
if (na>2),  error('too many input args'); end
newvals = cellfun(@(x) ~isempty(x), varargin);
optargs(newvals) = varargin(newvals);
[header info] = optargs{:};


% check size
[rowsx,colsx] = size(x);
[rowsy,colsy] = size(y);

if( rowsx ~= rowsy )
  error('X and Y rows must be the same!')
end
if( colsx ~= 1 )
  error('Only one col in X is admittted!')
end


if info
    tstart = tic;
    fprintf(1,'====> WriteGPlot.m\n');
    fprintf(1,'writing X:(%d x %d) Y:(%d x %d) in %s ...',...
            rowsx,colsx,rowsy,colsy,fname);
end


% memory
d = zeros(colsy+1,rowsx);
d(1,:) = x;
d(2:end,:) = y';


% format
numformat = '%.16e ';
format = [repmat(numformat,1,colsy+1) '\n'];


% write
fid = fopen( fname, 'w' );
if ~isempty(header)
  fprintf( fid, '# %s\n', header );
end
fprintf( fid, format, d );
fprintf( fid, '\n' );
fclose(fid);


if info
  fprintf(1,' done\n');
  fprintf(1,'<==== (%4.3f s)\n',toc(tstart));
end
