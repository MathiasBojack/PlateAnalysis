Input.Geometry.a    = 12;
Input.Geometry.b    = 12;
Input.Geometry.h    = 0.15;
Input.Properties.Nu = 0.2;
Input.Properties.E  = 19.2e9;
Input.Load.Time     = 7200;
Input.Load.gamma    =  2500*9.8*Input.Geometry.h;
Input.SolParam.M    = 20;
Input.SolParam.N    = 20;

Input = dim_Prop_Cal( Input);

Solution = Four_ss_rectangular_Plate_VK( Input );
%%
f1 = figure(1);
s1= subplot(1,2,1);
colormap(jet);
H=surf(Solution.X,Solution.Y,Solution.W);
shading interp
set(gcf,'color','white')

xx = get(H, 'XData');
zz = get(H, 'Zdata');
set(H, 'XData', zz, 'ZData', xx);

set(H, 'CData', get(H, 'YData'))
set(H, 'CData', get(H, 'ZData'))
set(H, 'CData', get(H, 'XData')) 

set(H, 'CData', get(H, 'XData')) 
ax = gca;
ax.LineWidth = 1;
ax.DataAspectRatio = [0.1 1 1];
axis tight
xlabel('$W$ (/m)','Interpreter','LaTex','FontSize',12) 
ylabel('$X_2 $','Interpreter','LaTex','FontSize',12) 
zlabel('$X_1 $','Interpreter','LaTex','FontSize',12)
min(Solution.W(:))