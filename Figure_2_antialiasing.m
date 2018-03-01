clear
close all
addpath('Files')
%
c = 343.1;
v = 200;


N_ssd = 120;                                 % Number of loudspeakers
fi = (0:2*pi/N_ssd:2*pi-2*pi/N_ssd)';
Rssd = 2;
x0 = Rssd*[ cos(fi) sin(fi) ];

Rref = 0.1;
xref = Rref*[ cos(fi) sin(fi) ];

n0 = -[ cos(fi) sin(fi) ];
v0 =  [ sin(fi) -cos(fi)];
%
dx = 0.025;
x_field = (-3.5:dx:2.5);
y_field = (-2.5:dx:3.5 );
[X,Y] = meshgrid(x_field,y_field);

%
x_a = [   -3   -3  -3  -2.5   -1.5    0  2.5;    %x_coordinates
    -2.5  -1   0   1.5    2.5    3  3 ];  %y_coordinates
[ p,xp,yp ] = make_path( x_a(1,:), x_a(2,:), 150 );

Tsim = p(end)/v;
fs = 44.1e3;
t = (0:1/fs:Tsim-1/fs)';
Nt = length(t);
w = 2*pi*fftshift( (-Nt/2:Nt/2-1)/(Nt)*fs );
s = repmat([ zeros(5,1);tukeywin(15,0.5);zeros(105,1)],100,1)';
s = s(1:length(t));

fir = fir1(100,[800,5e3]/fs*2);
s = filter(fir,1,s,[]);

xs = get_trajectory( p,xp,yp, t, v );
%% Calculate driving functions
Nt = length(t);
w = 2*pi*fftshift( (-Nt/2:Nt/2-1)/(Nt)*fs );
d_wfs = zeros(length(t),length(x0));
s_wfs = ifft( sqrt(1i*w/(c*2*pi)).*fft(s) );
%
dx0 = mean(sqrt(sum(diff(x0,1).^2,2)));
Tau = get_initial_position( v,c, x_a, x0 );
wb = waitbar(0,'Calculating driving functions');

wc = zeros(length(t),length(x0));
for n = 1 : length(t)
    waitbar(n/length(t),wb);
    xs_t = interp1( t,xs(:,1), t(n)-Tau, 'linear','extrap' );
    ys_t = interp1( t,xs(:,2), t(n)-Tau, 'linear','extrap' );
    vx = (xs_t - interp1( t,xs(:,1), t(n)-1/fs-Tau, 'linear','extrap' ))*fs;
    vy = (ys_t - interp1( t,xs(:,2), t(n)-1/fs-Tau, 'linear','extrap' ))*fs;
    Dvx = x0-[xs_t ys_t];
    R = sqrt( sum( Dvx.^2,2) );
    Vv = 1/c*sum([vx vy].*(Dvx),2);
    Delta = R - Vv;
    D = xs_t.*x0(:,2) - ys_t.*x0(:,1);
    xref = ( D.*Dvx(:,2) - sign(Dvx(:,2)).*Dvx(:,1).*sqrt(Rref^2.*R.^2-D.^2 ))./R.^2;
    yref = (-D.*Dvx(:,1) - abs(Dvx(:,2)).*sqrt(Rref^2.*R.^2-D.^2 ))./R.^2;
    dref = sqrt( sum( (x0-[xref yref]).^2, 2)  );
    Kn = sum(n0.*Dvx ,2);
    Kn =  Kn.*(Kn>0);
    d_wfs(n,:) = sqrt(dref./(dref+R)).*Kn.*interp1(t,s_wfs,t(n)-Tau,'linear')./Delta.^(3/2)*dx0;
    kt = sum(Dvx./R.*v0,2);    wc(n,:) = pi/dx0*c./abs(kt);
    
    Tau = Tau - 1/fs*Vv./Delta;
    
    
end
close(wb);
d_wfs(isnan(d_wfs)) = 0;
%% Antialising filtering
wlen = 2048/32;
hop = wlen/4;
nfft = wlen;
d_wfs_aa  = anti_aliasing_fd(fs,x0,d_wfs,t,wc,wlen,hop,nfft);

%%
field_synth    = zeros(size(X));
field_synth_aa = zeros(size(X));
i=1;
Rfield_full = reshape( sqrt((bsxfun( @minus, X(:), x0(:,1)' )).^2+ (bsxfun( @minus, Y(:), x0(:,2)' )).^2),...
    length(y_field), length(x_field) , length(x0));

%%
ftsize = 13;
f = figure('Units','points','Position',[200,100,280,550]);
set(f,'defaulttextinterpreter','latex')

fig_pos = [ 0.125   0.575  0.82 .42;
            0.125   0.0725 0.82 .42];

sp1 = axes('Units','normalized','Position',fig_pos(1,:));
p1 = pcolor(sp1, x_field,y_field,real(field_synth));
shading interp
axis equal tight
caxis([-1,1]*2.5e-2)
hold on
plot(sp1, xs(:,1),xs(:,2),'--k','LineWidth',2)
xlim([x_field(1),x_field(end)])
ylim([y_field(1),y_field(end)])
draw_ssd( sp1, x0(1:1:end,:), n0(1:1:end,:), 0.03 )
pos1 = plot( sp1, xs(1,1),xs(1,2),'ok','MarkerFaceColor','white');
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-2);

sp2 = axes('Units','normalized','Position',fig_pos(2,:));
p2 = pcolor(sp2, x_field,y_field,real(field_synth_aa));
shading interp
axis equal tight
caxis([-1,1]*2.5e-2)
hold on
plot(sp2, xs(:,1),xs(:,2),'--k','LineWidth',2)
xlim([x_field(1),x_field(end)])
ylim([y_field(1),y_field(end)])
draw_ssd( sp2, x0(1:1:end,:), n0(1:1:end,:), 0.03 )
pos2 = plot( sp2, xs(1,1),xs(1,2),'ok','MarkerFaceColor','white');
xlabel( '$x \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
ylabel( '$y \rightarrow [\mathrm{m}]$' , 'FontSize', ftsize );
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
b = get(gca,'XTickLabel');
set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize-2);
%%
%
[~,i0] = min(abs(t-0.0237));
for i = 1 : i0
    field_synth    = zeros(size(X));
    field_synth_aa = zeros(size(X));
    for n = 1 : length(x0)
        Rfield_synth = squeeze( Rfield_full(:,:,n) );
        field_synth    = field_synth    + 1/(4*pi)*interp1( t, d_wfs(:,n)   , t(i)-Rfield_synth/c,'linear','extrap'   )./Rfield_synth;
        field_synth_aa = field_synth_aa + 1/(4*pi)*interp1( t, d_wfs_aa(:,n), t(i)-Rfield_synth/c,'linear','extrap'   )./Rfield_synth;
    end
    %
    set(p1,'CData',real(field_synth));
    delete(pos1);
    pos1 = plot(sp1, xs(i,1),xs(i,2),'ok','MarkerFaceColor','white');
    
    set(p2,'CData',real(field_synth_aa));
    delete(pos2);
    pos2 = plot(sp2, xs(i,1),xs(i,2),'ok','MarkerFaceColor','white');
    drawnow
    
end
%%
set(gcf,'PaperPositionMode','auto');
print( '-r300', 'antialiased_synth' ,'-dpng')