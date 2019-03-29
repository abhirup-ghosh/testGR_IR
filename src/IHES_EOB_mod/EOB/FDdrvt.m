function [d1f, d2f] = FDdrvt(f,t, varargin)

%FDdrvt Compute finite differences of the function f(t).
%
%   [d1f, d2f] = FDdrvt(f,t) compute finte differences of f(t) 
%
%   [d1f, d2f] = FDdrvt(f,t, order) specify order, a negative order uses
%   expressions for data non uniformly sampled (only 1st derivative
%   implemented at the moment)  
%

% Manage args in 
order   = 2;  % order
optargs = {order};
na = length(varargin);
if (na>1), error('too many input args'); end
newvals = cellfun(@(x) ~isempty(x), varargin);
optargs(newvals) = varargin(newvals);
[order] = optargs{:};

n     = length(t);

% check for uniform spacing
%FIXME: better way?
if ( all(diff(diff(t))) )
    order = - order;
else
    dt    = t(2)-t(1);
    oodt  = 1./dt;   
    oodt2 = oodt^2;        
end
 
t   = reshape(t,size(f));
df  = 0*f;
d2f = 0*f;

switch order
  
 case 2
  
  i = 2:n-1;
  d1f(i) =  (f(i+1)-f(i-1))*0.5*oodt;
  d2f(i) =  (f(i+1) - 2*f(i) + f(i-1))*oodt2;
  
  i = 1;
  d1f(i) = -(3*f(i)-4*f(i+1)+f(i+2))*0.5*oodt;
  d2f(i) =  (2*f(i)-5*f(i+1)+4*f(i+2)-f(i+3))*oodt2;
  
  i = n;
  d1f(i) =  (3*f(i)-4*f(i-1)+f(i-2))*0.5*oodt;
  d2f(i) =  (2*f(i)-5*f(i-1)+4*f(i-2)-f(i-3))*oodt2;
   
 case 4 
  
  c = 1/12;

  i = 3:n-2;
  d1f(i) = c*(8*(f(i+1)-f(i-1)) - f(i+2) + f(i-2))*oodt;
  d2f(i) = c*(-30*f(i)+16*(f(i+1)+f(i-1))-(f(i+2)+f(i-2)))*oodt2;
  
  i = 1;
  d1f(i) = c*(-25*f(i)+48*f(i+1)-36*f(i+2)+16*f(i+3)-3*f(i+4))*oodt;
  d2f(i) = c*(45*f(i)-154*f(i+1)+214*f(i+2)-156*f(i+3)+61*f(i+4)-10*f(i+5))*oodt2;
  
  i = 2;
  d1f(i) = c*(-3*f(i-1)-10*f(i)+18*f(i+1)-6*f(i+2)+f(i+3))*oodt;
  d2f(i) = c*(10*f(i-1)-15*f(i)-4*f(i+1)+14*f(i+2)-6*f(i+3)+f(i+4))*oodt2;
  
  i = n-1;
  d1f(i) = - c*(-3*f(i+1)-10*f(i)+18*f(i-1)-6*f(i-2)+f(i-3))*oodt;
  d2f(i) = c*(10*f(i+1)-15*f(i)-4*f(i-1)+14*f(i-2)-6*f(i-3)+f(i-4))*oodt2;
  
  i = n;
  d1f(i) = - c*(-25*f(i)+48*f(i-1)-36*f(i-2)+16*f(i-3)-3*f(i-4))*oodt;
  d2f(i) = c*(45*f(i)-154*f(i-1)+214*f(i-2)-156*f(i-3)+61*f(i-4)-10*f(i-5))*oodt2;
    
 case -2
     
  %FIXME: 2nd drvt   
  if nargout>1
    warning('vectors are non-uniformly spaced, 2nd derivate not implemented yet');
  end
  
  i = 2:n-1;      
  df(i) = (f(i+1)-f(i-1))/(t(i+1)-t(i-1));
    
  i = 1;
  d1f(i) = -(3*f(i)-4*f(i+1)+f(i+2))./(t(i+2)-t(i));
  
  i = n;
  d1f(i) =  (3*f(i)-4*f(i-1)+f(i-2))./(t(i)-t(i-2));  
  
 case -4    
     
  %FIXME: 2nd drvt      
  if nargout>1
    warning('vectors are non-uniformly spaced, 2nd derivate not implemented yet');
  end
  
  c = 1/3;   
     
  i = 3:n-2;
  d1f(i) = c*(8*(f(i+1)-f(i-1)) - f(i+2) + f(i-2))./(t(i+2)-t(i-2));
    
  %FIXME: the following expressions need fixing:
  i = 1; 
  d1f(i) = c*(-25*f(i)+48*f(i+1)-36*f(i+2)+16*f(i+3)-3*f(i+4))./(t(i+4)-t(i));
  
  i = 2;
  d1f(i) = c*(-3*f(i-1)-10*f(i)+18*f(i+1)-6*f(i+2)+f(i+3))./(t(i+3)-t(i-1));
  
  i = n-1;
  d1f(i) = - c*(-3*f(i+1)-10*f(i)+18*f(i-1)-6*f(i-2)+f(i-3))./(t(i+1)-t(i-3));
  
  i = n;
  d1f(i) = - c*(-25*f(i)+48*f(i-1)-36*f(i-2)+16*f(i-3)-3*f(i-4))./(t(i)-t(i-4));    
    
 otherwise
  error(' order not implemented')
end



%{
% test
x=0:0.01:2*pi;f=sin(x);
[df1,df2] = FDdrvt(f,x,2); plot(x,f,x,cos(x),x,df1,x,-df2)
[df1,df2] = FDdrvt(f,x,4); plot(x,f,x,cos(x),x,df1,x,-df2)
%}
