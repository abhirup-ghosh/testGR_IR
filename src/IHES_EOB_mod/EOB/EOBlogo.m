function EOBlogo()


EOBopt.a5 = +23.5;
EOBopt.a6 = -121;

EOBopt.Dynamics = 'eob';
EOBopt.PNorder  = '5pnlog';
EOBopt.resumD   = 'pade03';

EOBopt.Tidal   = '';
EOBopt.PNTidal = '';

nu = 1/4;
r  = 1.1:(25-1.1)/255:20;
th = -pi:2*pi/63:pi;

[A dA] = EOBMetric( nu, r, EOBopt );

figure
plot(r,A)

[rg,tg]=meshgrid(r,th);

xg = rg.*cos(tg);
yg = rg.*sin(tg);

A = repmat(A,length(th),1);

figure
surfc(xg,yg,A,...
    'FaceColor','interp',...
    'EdgeColor','none',...
    'FaceLighting','phong',...
    'FaceAlpha',0.8)

axis off

camlight left
set(gca,'color','none'); % trasparent bckgrd
set(gca, 'visible', 'off'); % invisible axis
set(gcf,'Color','none');

minx = 0.9*min(min(xg));
maxx = 0.9*max(max(xg));
mina = 1.01*min(min(A));
maxa = 1.01*max(max(A));

axis([minx maxx minx maxx mina maxa]) 

%daspect([1 1 1])
%view(-30,36)
