function adim_Prop = adim_Prop_Cal( Time,h,nu )
%--------------------------------------------------------------------------
%
% This function caculates the adimensional section constants A0,B0,D0,NT0,MT0 
% for a plate with a thickness of h and temperature distribution caculated by
% Safir in advance
%
%--------------------------------------------------------------------------
% Temperature profile are registered every 60s from 60s to 7200s (2h)
Profil  = read_Temperature( h );

Safir_z = Profil.Position;    % vertical abscisse from -h/2:h/2 

index   = ismember(Profil.Temperature(1,:),Time);  
Safir_T = Profil.Temperature(2:end,index);            % temperature 

A0      = zeros(1,size(Safir_T,2));
B0      = zeros(1,size(Safir_T,2));
D0      = zeros(1,size(Safir_T,2));
NT0     = zeros(1,size(Safir_T,2));
MT0     = zeros(1,size(Safir_T,2));
% thermal curvature for simply supported sides.
CHI_SS  = zeros(1,size(Safir_T,2));
% thermal curvature for totally free boundary
CHI_FB  = zeros(1,size(Safir_T,2));


% Curve fitting Operation
% Fitting Method: 'smoothingspline','poly1' ect.
for i = 1:size(Safir_T,2)

        fitobject1  = fit(Safir_z,Safir_T(:,i),'smoothingspline');
        % Curve fitting object could be treated as a function.
        % Caculation of section constant
        dz          = 0.0001;
        z           = (-h/2:dz:h/2)';
        Temperature = fitobject1(z);
        kEc         = reduction_Factor(Temperature);
        epsilon     = thermal_Dilatation(Temperature);
        
        t       = z/h;
        dt      = dz/h;
        A0(i)   = sum(kEc*dt);
        B0(i)   = sum(kEc.*t*dt);
        D0(i)   = sum(kEc.*(t.^2)*dt);
        
        NT0(i)  = sum(kEc.*epsilon*dt);
        MT0(i)  = sum(kEc.*epsilon.*t*dt);
        % thermal curvature for simply supported boundary
        CHI_SS(i) = -(1+nu)/h*(A0(i)*MT0(i)-B0(i)*NT0(i))/(A0(i)*D0(i)-B0(i)^2);  
        % thermal curvature for totally free boundary
        CHI_FB(i) = -1/h*(A0(i)*MT0(i)-B0(i)*NT0(i))/(A0(i)*D0(i)-B0(i)^2);
end

adim_Prop.A0 = A0;
adim_Prop.B0 = B0;
adim_Prop.D0 = D0;
adim_Prop.NT0 = NT0;
adim_Prop.MT0 = MT0;
adim_Prop.CHI_SS = CHI_SS;
adim_Prop.CHI_FB = CHI_FB;
%-----------------------------------------------------------------------------------------------
%%% Dimensional properties

% Chi = (B0*NT0-A0*MT0)/(A0*D0-B0^2)/h*(1+nu)=Chi_T
% De = (A0*D0-B0^2)/A0*E*h^3;


