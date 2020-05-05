%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%        UPPER BOUND KINEMATIC APPROACH FOR SHELLS                    %%%%
%%%%         Jeremy Bleyer, Laboratoire Navier (13/01/2014)                      %%%%
%%%%                                                                                              %%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% clc
%%
% a               = 12;
% b               = 12;
% thickness_list  = 0.15;
% Time_list       = [60 1200 :1200:14400];

% model           = 'VK';
% model           = 'KL';

%% GENERAL INPUTS

% Normalized time for fire condition simulation

%%% TYPE OF THE STRENGTH CRITERION AND PLATE THICKNESS


% Normalization with M0 = sig0*thick^2/4 = 1
% critere de Mohr-Coulomb tronque en traction avec renforcements
crit.type = 'MCt-r';

% nombre de couches dans l'¨¦paisseur pour approcher le critere 
crit.ngz = 8;
crit.ft = 0;
crit.fc = 32; 
% nombre d aciers 
crit.ac.nbac = 4;
% direction des aciers
crit.ac.dir = [1 0;0 1;0 1;1 0];
%crit.ac.dir = [1 1;-1 1;-1 1;1 1]./sqrt(2);
crit.ac.Ns = 0.1414;                         % eta = pi*phi^2/4/c/h (c : espacement), Ns = fs*eta*h  0.25*pi*36*500/0.1/10^6
crit.ac.zs = [-0.56 -0.48 0.48 0.56];   % 2*2 lits orthogonaux (3 cm d'enrobage, phi = 6mm)

%% INPUT FILE NAME


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
mesh.stabilization = 'drilling-off';
% type d'element fini
% mesh.el_type = 'discontinuous-only';
mesh.el_type = 'discontinuous';

% nombre de points de Gauss par arete pour integrer Prm_disc
mesh.nged = 5;


crit.strength = ones(mesh.NNE,1);


%% BOUNDARY CONDITIONS
% 1 = fixed, 0 = free
% possible conditions are :     unorm = normal in-plane displacement
%                                      utang = tangential in-plane displacement
%                                      utrans = transversal displacement
%                                      rottan = rotation along tangent
%                                      rotnorm = rotation along normal
% mesh.bc.unorm = [0 1 0 0];
% mesh.bc.utang = [0 1 0 1];
% mesh.bc.utrans = [0 1 0 1 ];
% mesh.bc.rottan = [0 1 0 1 ];
% mesh.bc.rotnorm = [0 0 0 0];

% % simple support
mesh.bc.unorm = [1 1 0 0];
mesh.bc.utang = [0 0 0 0];
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



% %------------------------------------------------------------------------------------
% %       Clambed boundaries
% % ------------------------------------------------------------------------------------
% % displacement normal to the boundary
% mesh.bc.unorm = [1 1 0 0];
% % displacement tangent to the boundary
% mesh.bc.utang = [0 0 0 0];
% % out-of-plan displacement
% mesh.bc.utrans = [1 0 1 0];
% mesh.bc.rottan = mesh.bc.utrans;
% mesh.bc.rotnorm = [0 1 0 0];



% 
% mesh.bc.unorm = [1 1 0 0 1];
% mesh.bc.utang = [0 0 0 0 1];
% mesh.bc.utrans = [1 0 1 1 1];
% mesh.bc.rottan = mesh.bc.utrans;
% mesh.bc.rotnorm = [0 1 0 0 1];

%% FORMING DOFS 
mesh = make_dof(mesh);
mesh.NED = sum(mesh.id_dof.nb_disc);      % total number of active edges (internal + boundary except free edges)
mesh.ndof = 27;
mesh.Nu = mesh.ndof*mesh.NNE;                       % total number of d.o.f (3x6 displacements + 3x3 rotations)




 %% PLOTTING MESH
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
    clear X;
    load(deformee_name);
    X(:,1) = reshape(Solution.W,[length(Solution.W(:,1))^2,1]);
    X(:,2) = reshape(Solution.Y,[length(Solution.W(:,1))^2,1]);
    X(:,3) = reshape(Solution.X,[length(Solution.W(:,1))^2,1]);
    
    mesh.coor(:,1) = griddata(X(:,2),X(:,3),X(:,1),mesh.coor(:,2),mesh.coor(:,3));
    
    figure(1001)
    clf
    plot_mesh(mesh)
    
    Profil = read_Temperature( thickness );
    z = Profil.Position;
    T = Profil.Temperature(2:end,Profil.Temperature(1,:)==Time);


    
    [kc,ka,kt] = fire_red_coeff(T);
    [az,wz] = quadrature_1D(crit.ngz,'uniform');
    crit.red_fact = interp1(z,kc,thickness/2.*az);    
    crit.red_fact_t = interp1(z,kt,thickness/2.*az);    
    crit.ac.red_fact = interp1(z,ka,thickness/2.*crit.ac.zs);
   
    mesh = average_normals(mesh);


    nodes = mesh.edges(mesh.edges(:,3)==1,1:2);
    for j=1:size(nodes,1)
        n1 = nodes(j,1);n2=nodes(j,2);
        mesh.edg_norm(mesh.act_edges.all.node2edg(n1,n2),:) = [1 0 0];
    end
     nodes = mesh.edges(mesh.edges(:,3)==2,1:2);
    for j=1:size(nodes,1)
        n1 = nodes(j,1);n2=nodes(j,2);
        mesh.edg_norm(mesh.act_edges.all.node2edg(n1,n2),2) = 0;
%        mesh.edg_norm(mesh.act_edges.all.node2edg(n1,n2),:) = [1 0 0];
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
    
%% COMPUTING THE GLOBAL LOAD VECTOR
loading = 'unif';
%loading = 'hydrostatic';
% loading = 'traction';
%loading = 'periodic';

switch loading
    case 'unif'
        p0 = 2500*9.81*thickness*1e-6;
        p = sparse(mesh.NNE,1);
        p(mesh.el_tag==1,3) = -p0;
        Fu = assembl_unif_pressure(mesh,p,'global');
    otherwise
        error('TODO');
end




%% COMPUTING THE GLOBAL EQUILIBRIUM MATRIX
B = strain_matrix(mesh);


%% COMPUTING THE GLOBAL DISCONTINUITY MATRIX
[Disc,edg_lgth,normal_loc] = discontinuity(mesh);



%% COMPUTING THE INDEX LIST FOR BCs
BC = boundary_conditions(mesh);

% uniquement pour les joints (physical_line = 5)
list_node = mesh.edges(mesh.edges(:,3)==5,1:2);
suppl_var = [];
for kk = 1:size(list_node,1)
    id = mesh.act_edges.prm.node2edg(list_node(kk,1),list_node(kk,2));
    suppl_var = [suppl_var;id];
end

    % Fonction Mosek sans joints
    prob = building_mosek(crit,thickness,Fu,BC,B,Disc,edg_lgth,normal_loc,mesh);
    
    % Fonction Mosek en presence de joints
%     prob = building_mosek_charn(crit,thick,Fu,BC,B,Disc,edg_lgth,normal_loc,mesh,suppl_var);
    
    dual = primal2dual(prob);
   
    param.MSK_IPAR_PRESOLVE_USE = 'MSK_PRESOLVE_MODE_OFF';
    param.MSK_IPAR_PRESOLVE_LINDEP_USE = 'MSK_OFF';
    
%     addpath('D:\Mosek\8\toolbox\r2014a');
    [r,res] = mosekopt('maximize',dual,param);
%     rmpath('D:\Mosek\8\toolbox\r2014a');

    x= res.sol.itr.y;    
    
   % interpolation quadratique par element
   Ufine = reshape(x(1:mesh.Nu),[27 mesh.NNE]);
   % interpolation lineaire par element
   U = Ufine(1:9,:);
   U = U(:);
   Theta = Ufine(19:27,:);
   Theta = Theta(:);
   Ufine = Ufine(1:18,:);
   Ufine = Ufine(:);
   
   % Calcul des courbures chi et des deformations membranaires eps
   d = B*x(1:mesh.Nu);
   chixx = reshape(repmat(d(1:18:18*mesh.NNE)',3,1),3*mesh.NNE,1);
   chiyy = reshape(repmat(d(2:18:18*mesh.NNE)',3,1),3*mesh.NNE,1);
   chixy = reshape(repmat(d(3:18:18*mesh.NNE)',3,1),3*mesh.NNE,1);
   list_eps = repmat(18*(0:mesh.NNE-1),9,1) + repmat((4:12)',1,mesh.NNE);
   epsxx = d(list_eps(1:3:end));
   epsyy = d(list_eps(2:3:end));
   epsxy = d(list_eps(3:3:end));
   
   ll = setdiff((1:15*mesh.NED)',BC);
   dd = Disc(ll,:)*x(1:mesh.Nu);

   % Post-process de la solution qui recalcule avec une meilleure precision
   % Prm et donc lambda
    ngz2 = 41;
    [az,wz] = quadrature_1D(ngz2,'uniform');
    crit.red_fact = interp1(z,kc,thickness/2.*az);  
    crit.red_fact_t = interp1(z,kt,thickness/2.*az);    
    crit.ac.red_fact = interp1(z,ka,thickness/2.*crit.ac.zs);
   [Prm_strain,Prm_disc] = prm_post_process(d,dd,BC,crit,thickness,edg_lgth,mesh,ngz2,41,normal_loc);
   lamb_kin(i_Time,i_thick) = Prm_strain + Prm_disc;
   disp('lamb_kin=');   lamb_kin
   % affiche et exporte le mecanisme de ruine
    cnt = (i_thick-1)*length(Time_list) +i_Time;
    fig = figure(cnt);
%    for alp = 0:0.3./max(abs(Ufine(1:3:end))):3./max(abs(Ufine(1:3:end)));
    alp = 3./max(abs(Ufine(1:3:end)));
    clf
    hold on
    plot_field(Ufine(1:3:end)./max(-Ufine(:)),Ufine(:),alp,mesh) 
    mesh2 = mesh;
    mesh2.coor(:,2) = 2*max(mesh.coor(:,2))-mesh.coor(:,2);
    Ufine2 = Ufine(:);
    Ufine2(2:3:end) = -Ufine(2:3:end);
    plot_field(Ufine(1:3:end)./max(-Ufine(:)),Ufine2(:),alp,mesh2) 
    view(-56,14)
    axis tight
    ylim([0 b])
    z_initial = zlim;
    zlim([z_initial(1) z_initial(1)+a])
%     xlim([-4 4])
%     ylim([-1 13])
%     zlim([0 1.1.*max(mesh.coor(:,3))])
    
%     title('Mur 12 m x 12 m','interpreter','latex','FontSize',12)
        
%     set(gca,'YTick',0:2:12)
%     set(gca,'ZTick',0:2.4:max(mesh.coor(:,3)))
%     set(gcf, 'PaperUnits', 'inches');
%     set(gcf, 'PaperSize', [5.2 5.8]);
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf, 'PaperPosition', [0 0 5.2 5.8]);
    switch model
        case 'KL'
            fig_name = strcat('.\7.CaseStudy\4sPlate\kinematic\KLmodel\a_',...
                num2str(a),'_b_',num2str(b),'_h_',num2str(thickness),...
                '_Time_',num2str(Time),'_KL_Kine.tiff');
        case 'VK' 
            fig_name = strcat('.\7.CaseStudy\4sPlate\kinematic\VKmodel\a_',...
                num2str(a),'_b_',num2str(b),'_h_',num2str(thickness),...
                '_Time_',num2str(Time),'_VK_Kine.tiff');
    end
    saveas(fig,fig_name)
    close
%     print('-dpng','-r200',horzcat('mur_12_12_',int2str(i),'.png'))
end
end



switch model
    case 'KL'
        lambda_name =  strcat('.\7.CaseStudy\4sPlate\kinematic\KLmodel\lambda_a_',...
            num2str(a),'_b_',num2str(b),'_h_',num2str(thickness_list),'_KL.mat');
    case 'VK' 
        lambda_name =  strcat('.\7.CaseStudy\4sPlate\kinematic\VKmodel\lambda_a_',...
            num2str(a),'_b_',num2str(b),'_h_',num2str(thickness_list),'_VK.mat');
end

save(lambda_name,'lamb_kin','Time_list')