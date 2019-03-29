function PrintSFields( s, fid )

%PrintSFields print fields of a structure and their content 
%
%   PrintSFields( s, fid )
%
    

Name = fieldnames(s);
for k=1:length(Name)
  pr(fid,Name{k},s.(Name{k}));    
end


function pr(fid,name,value)
if isa(value,'char')
  fprintf(fid,' %-20s = %s\n',name,value);  
elseif isa(value,'double')
  if isscalar(value)
    fprintf(fid,' %-20s = %.12e\n',name,value);
  else
    [n,m]= size(value);
    fprintf(fid,' %-20s = <%s> %dx%d\n',name,class(value),n,m);
  end
elseif isa(value,'function_handle')
  fprintf(fid,' %-20s = %s\n',name,func2str(value));
else
  fprintf(fid,' %-20s = <%s>\n',name,class(value));
end

