%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          LOWER BOUND STATIC APPROACH FOR THICK SHELLS             %%%%
%%%%         Jeremy Bleyer, Laboratoire Navier (02/05/2013)            %%%%
%%%%                                                                   %%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all


%%
% Input.Geometry.a    = 12;
% Input.Geometry.b    = 10;
% Input.Geometry.h    = 0.15;
% Input.Properties.Nu = 0.2;
% Input.Properties.E  = 19.2e9;
% Input.Load.gamma    =  2500*9.8*Input.Geometry.h;

% a               = Input.Geometry.a;
% b               = Input.Geometry.b;
% thickness_list  = Input.Geometry.h;
a               = 14;
b               = 14;
thickness_list  = 0.15;
% Time_list       = [60 1200 :1200:14400];
Time_list       = [5400:600:8400];

% model           = 'VK';
model           = 'KL';

%% 
global dim

dim = 3;

% epaisseur
% critere de Mohr-Coulomb tronque en traction avec renforcements
crit.type = 'MCt-r';
% nombre de couches dans l'épaisseur pour approcher le critere 
crit.ngz = 8;
crit.ft = 0.;
crit.fc = 32;
% nombre d aciers 
crit.ac.nbac = 4;
% direction des aciers
crit.ac.dir = [1 0;0 1;0 1;1 0];
% effort membranaire maxi des aciers par unite de longueur
crit.ac.Ns = 0.1414;                         % eta = pi*phi^2/4/c/h (c : espacement), Ns = fs*eta*h
% position normalisee des aciers dans l'epaisseur [-1;1]
crit.ac.zs = [-0.56 -0.48 0.48 0.56];   % 2*2 lits orthogonaux (3 cm d'enrobage, phi = 6mm)


%% USER INTERFACE
%   TODO

% File name

fname = strcat('Geometry_a_',num2str(a),'_b_',num2str(b));

switch model
    case 'KL'
        Deformee_Path = '.\3.Thermoelastic\Kirchhoff-Love\DeformedShape';
    case 'VK' 
        Deformee_Path = '.\3.Thermoelastic\vonKarman\DeformedShape';
end


for i_thick = 1:length(thickness_list)
    thickness = thickness_list(i_thick);

%% CALLING GMSH FOR MESHING
mesh = read_input(fname);

%% BOUNDARY CONDITIONS
% 1 = fixed, 0 = free
% Condition Order = [Bottom Axis of symmetry Top Left]


% mesh.bc.unorm = [1 1 0 1];
% mesh.bc.utang = [1 1 0 1];
% mesh.bc.utrans = [1 1 0 1 ];
% mesh.bc.rottan = [1 1 0 1 ];
% mesh.bc.rotnorm = [0 1 0 1];



%------------------------------------------------------------------------------
%       simply supported boundary conditions on the four sides
% ----------------------------------------------------------------------------
% displacement normal to the boundary
mesh.bc.unorm = [1 1 0 0];
% displacement tangent to the boundary
mesh.bc.utang = [0 0 0 0];
% out-of-plan displacement
mesh.bc.utrans = [1 0 1 1];
mesh.bc.rottan = mesh.bc.utrans;
mesh.bc.rotnorm = [0 1 0 0];

% %------------------------------------------------------------------------------------
% %       simply supported on the top and bottom and free on the lateral
% %       sides
% % ------------------------------------------------------------------------------------
% % displacement normal to the boundary
% mesh.bc.unorm = [1 1 0 0];
% % displacement tangent to the boundary
% mesh.bc.utang = [0 0 0 0];
% % out-of-plan displacement
% mesh.bc.utrans = [1 0 1 0];
% mesh.bc.rottan = mesh.bc.utrans;
% mesh.bc.rotnorm = [0 1 0 0];



mesh = make_dof(mesh);
mesh.NED = mesh.act_edges.all.NED;



    


 %% PLOTTING MESH1
% figure(1000)
% clf
% plot_mesh(mesh)



for i_Time = 1:length(Time_list)  
    Time = Time_list(i_Time);

    switch model
        case 'KL'
            deformee_mat = strcat('Deformee_a_',num2str(a),'_b_',...
                num2str(b),'_h_',num2str(thickness),'_Time_',...
                num2str(Time),'_KL.mat');
        case 'VK' 
            deformee_mat = strcat('Deformee_a_',num2str(a),'_b_',...
                num2str(b),'_h_',num2str(thickness),'_Time_',...
                num2str(Time),'.mat');
    end
    deformee_name = fullfile(Deformee_Path,deformee_mat);
    load(deformee_name);
    X(:,1) = reshape(Solution.W,[length(Solution.W(:,1))^2,1]);
    X(:,2) = reshape(Solution.Y,[length(Solution.W(:,1))^2,1]);
    X(:,3) = reshape(Solution.X,[length(Solution.W(:,1))^2,1]);
    mesh.coor(:,1) = griddata(X(:,2),X(:,3),X(:,1),mesh.coor(:,2),mesh.coor(:,3));
    figure(1001)
    clf
    plot_mesh(mesh)
% %     load deformee_12x12_param.mat;
%     % points configuration deformee (Toan)
%     X = def{i};
%     % interpolation sur le maillage
%     mesh.coor(:,1) = griddata(X(:,2),X(:,3),X(:,1),mesh.coor(:,2),mesh.coor(:,3));
    
    Profil = read_Temperature( thickness );
    z = Profil.Position;
    T = Profil.Temperature(2:end,Profil.Temperature(1,:)==Time);

    

    % calcul des coefficients de reduction suivant EC2
    [kc,ka,kt] = fire_red_coeff(T);
    % ordonnees (normalisees) du milieu des courches 
    [az,wz] = quadrature_1D(crit.ngz,'uniform');
    % reinterpolation des coeffs de reduction en ces points
    crit.red_fact = interp1(z,kc,thickness/2.*az);    
    crit.red_fact_t = interp1(z,kt,thickness/2.*az);    
    crit.ac.red_fact = interp1(z,ka,thickness/2.*crit.ac.zs);



   mesh = average_normals(mesh);

   % changement de la normale sur les facettes situees au bord du mur pour
   % imposition des conditions aux limites
    nodes = mesh.edges(mesh.edges(:,3)==1,1:2);
    for j=1:size(nodes,1)
        n1 = nodes(j,1);n2=nodes(j,2);
        mesh.edg_norm(mesh.act_edges.all.node2edg(n1,n2),:) = [1 0 0];
    end
     nodes = mesh.edges(mesh.edges(:,3)==2,1:2);
    for j=1:size(nodes,1)
        n1 = nodes(j,1);n2=nodes(j,2);
        mesh.edg_norm(mesh.act_edges.all.node2edg(n1,n2),2) = 0;
    end
    nodes = mesh.edges(mesh.edges(:,3)==3,1:2);
    for j=1:size(nodes,1)
        n1 = nodes(j,1);n2=nodes(j,2);
        mesh.edg_norm(mesh.act_edges.all.node2edg(n1,n2),:) = [1 0 0];
    end
    nodes = mesh.edges(mesh.edges(:,3)==4,1:2);
    for j=1:size(nodes,1)
        n1 = nodes(j,1);n2=nodes(j,2);
        mesh.edg_norm(mesh.act_edges.all.node2edg(n1,n2),:) = [1 0 0];
    end
mesh.edg_norm = mesh.edg_norm./(sqrt(mesh.edg_norm(:,1).^2 + mesh.edg_norm(:,2).^2+mesh.edg_norm(:,3).^2)*[1,1,1]);
       
    
    loading = 'unif';

    switch loading
        case 'unif'
            p = sparse(mesh.NNE,3);
            p(:,3) = -2500*9.81*thickness*1e-6;
            Fu = assembl_unif_load(mesh,p,'global');
            Fl = assembl_lin_load([],[],mesh);
        case 'line'
            Fu = sparse(9*mesh.NNE,1);
            pos = find(mesh.edges3(:,3)== line_idx);
            line = mesh.edges3(pos,1:2);
            Fl = assembl_lin_load(F,line,mesh);
        otherwise
            error('TODO');
    end
    
    % construction des matrices d'equilibre et de continuite
    H = equilibrium(mesh);
    Cont = continuity(mesh);
    
   % uniquement pour les joints (physical line = 5)
    list_node = mesh.edges(mesh.edges(:,3)==5,1:2);
    id_list = [];
    for kk = 1:size(list_node,1)
        id_list = [id_list;mesh.act_edges.all.node2edg(list_node(kk,1),list_node(kk,2))];
    end
    Liais = liaison(mesh,id_list);
    Cont = [Cont;Liais];
    Fl = [Fl;sparse(size(Liais,1),1)];


    NNE = mesh.NNE;
    Nu = mesh.NEQ;


    % local abscissae of checking points
    check_loc = [1 0;0 1;0 0;1/2 1/2;0 1/2;1/2 0;1/3 1/3;2/3 1/6;1/6 2/3;1/6 1/6];
    NM = sparse(10,6);
    NNV = sparse(10,3);
    for j=1:10
        NM(j,:) = quadratic_shape_fct(check_loc(j,1),check_loc(j,2))';
        NNV(j,:) = [check_loc(j,:) 1-check_loc(j,1)-check_loc(j,2)];
    end
    % number of checking points for the bending criterion
    Ncheck = 10;
    list_N = sparse(9*NNE,1);
    list_M = sparse(18*NNE,1);
    list_V = sparse(6*NNE,1);
    for e=1:NNE
        list_N(9*(e-1)+(1:9)) = (33*(e-1)+(1:9)); 
        list_M(18*(e-1)+(1:18)) = (33*(e-1)+9+(1:18)); 
        list_V(6*(e-1)+(1:6)) = (33*(e-1)+9+(19:24)); 
    end

    normalize = thickness^2/4;
  %   normalize = 1;
    
    clear prob
    prob = building_mosek(crit,thickness,Ncheck,NM,NNV,Fu,Fl,H,Cont,mesh);
    prob.c = prob.c.*normalize;
    param.MSK_IPAR_PRESOLVE_USE = 'MSK_PRESOLVE_MODE_ON';
    param.MSK_IPAR_PRESOLVE_LINDEP_USE = 'MSK_OFF';
    
%     addpath('D:\Mosek\8\toolbox\r2014a');
    [r,res] = mosekopt('maximize',prob,param);
%     rmpath('D:\Mosek\8\toolbox\r2014a');

    lamb_stat(i_Time,i_thick) = res.sol.itr.pobjval;

% x = res.sol.itr.xx;
% N = x(1+list_N);
% Nxx = N(1:3:end);
% Nyy = N(2:3:end);
% Nxy = N(3:3:end);
% M = x(1+list_M);
% Mxx = M(1:3:end);
% Myy = M(2:3:end);
% Mxy = M(3:3:end);
% V = x(1+list_V);
% Vx = V(1:2:end);
% Vy = V(2:2:end);
% Mxx_n = [Mxx(1:6:end)';Mxx(2:6:end)';Mxx(3:6:end)'];
% Mxx_n = Mxx_n(:);
% Myy_n = [Myy(1:6:end)';Myy(2:6:end)';Myy(3:6:end)'];
% Myy_n = Myy_n(:);
% Mxy_n = [Mxy(1:6:end)';Mxy(2:6:end)';Mxy(3:6:end)'];
% Mxy_n = Mxy_n(:);
% 
% s11 = x(1+Nu+(1:9:9*crit.ngz*Ncheck*NNE));
% s11 = reshape(s11,crit.ngz,Ncheck*NNE)';
% s22 = x(1+Nu+(2:9:9*crit.ngz*Ncheck*NNE));
% s22 = reshape(s22,crit.ngz,Ncheck*NNE)';
% s12 = x(1+Nu+(3:9:9*crit.ngz*Ncheck*NNE));
% s12 = reshape(s12,crit.ngz,Ncheck*NNE)';
% 
% n = Nxx*thick/4;
% 
% scatter(s11(:),s12(:),'.')
% hold on

% n = mean(reshape(Nxx,3,mesh.NNE))'*thick/4;
% Mxx_n = mean(reshape(Mxx_n,3,mesh.NNE))';
% figure(k)
% % clf
% ch = convhull([n Mxx_n]);
% scatter(n,Mxx_n,'.')
% hold on
% lamb = compute_MCt_mosek(crit,4,time);
% alp=0:0.02:2*pi;
% D = [cos(alp)' sin(alp)'];
% V = [lamb'.*D(:,1)-0.1,lamb'.*D(:,2)];
% plot(V(:,1),V(:,2),'-+m','LineWidth',2)
% axis equal
% xlim([-1.2 0.2])
% ylim([-0.6 0.6])


 
NED_rt = mesh.act_edges.rottan.NED;
NED_rn = mesh.act_edges.rotnorm.NED;
NED_utr = mesh.act_edges.utrans.NED;
NED_utg = mesh.act_edges.utang.NED;
NED_un = mesh.act_edges.unorm.NED;

U = pseudo_mech_cont(res.sol.itr.y(size(H,1)+3*(NED_rt+NED_rn)+(1:2*(NED_utr+NED_utg+NED_un))),mesh);
[Utr,Utang] = pseudo_mech(res.sol.itr.y(1:size(H,1)),mesh);

cnt = (i_thick-1)*length(Time_list) +i_Time;
fig = figure(cnt);
clf
plot_mesh_def_glob(mesh,-U,0.1);
view(-56,14)
axis tight
switch model
    case 'KL'
        fig_name = strcat('.\7.CaseStudy\4sPlate\KLmodel\a_',...
            num2str(a),'_b_',num2str(b),'_h_',num2str(thickness),...
            '_Time_',num2str(Time),'_KL.tiff');
    case 'VK' 
        fig_name = strcat('.\7.CaseStudy\4sPlate\VKmodel\a_',...
            num2str(a),'_b_',num2str(b),'_h_',num2str(thickness),...
            '_Time_',num2str(Time),'_VK.tiff');
end

saveas(fig,fig_name)
close
% % 
% figure(1)
% clf
% plot_elem_V(Nxx,mesh)
% view(120,7)
% 
% figure(2)
% clf
% plot_elem_M(Myy,mesh)
% view(120,7)



end
end

switch model
    case 'KL'
        lambda_name =  strcat('.\7.CaseStudy\4sPlate\static\KLmodel\lambda_a_',...
            num2str(a),'_b_',num2str(b),'_h_',num2str(thickness_list),'_KL2.mat');
    case 'VK' 
        lambda_name =  strcat('.\7.CaseStudy\4sPlate\static\VKmodel\lambda_a_',...
            num2str(a),'_b_',num2str(b),'_h_',num2str(thickness_list),'_VK2.mat');
end

save(lambda_name,'lamb_stat','Time_list')
